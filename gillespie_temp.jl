
using SpecialFunctions, Revise, Random, OffsetArrays, BenchmarkTools, Statistics, ProgressMeter, Profile
Random.seed!(1234)

include("C:\\Users\\franc\\OneDrive\\Desktop\\TESI_PELLIZZOLA\\SIMULAZIONI IN JULIA\\oggetti.jl")

#########################################################################
    #AUSILIARIE
#########################################################################
    function checktol(densità, total_times; tol=0.001)
        diff = norm.(datacoeph[1:end-1] .- datacoeph[2:end], Inf) # differenza massima in valore assoluto tra 
        flag = false
        for i in 1:length(diff)
            if diff[i] < tol
                println("Il sistema ha raggiunto la stazionarietà dopo $(total_times[i]) unità di tempo.")
                flag = true
            end

        end
        if !flag
            println("Il sistema non ha raggiunto la stazionarietà entro il tempo totale simulato.")
        end
    end


#########################################################################
    # INERENTI ALLA F GILLESPIE PRINCIPALE
##########################################################################
    function Tasep2Sys(tasep, ρ0)
        # estendo ν aggiungendo i serbatoi in 0 e L+1
        ν_data = BitVector(rand() < ρ0 for n in 0:tasep.L+1)
        ν = OffsetArray(ν_data, 0:tasep.L+1)
        # in 0 c'è sempre una particella (1) e in L+1 è sempre vuoto (0)
        ν[0] = true
        ν[tasep.L+1] = false
        vec = OffsetArray([ν[n] * (1 - ν[n+1]) for n in 0:tasep.L], 0:tasep.L)
        # estendo q aggiungendo i rate di ingresso e uscita α e β
        tasep.q[0] = tasep.α
        tasep.q[end] = tasep.β
        # rates
        W = tasep.q .* vec

        return ν, W, tasep.q
    end

    function onestep_gillespie(C::OffsetVector{Float64,Vector{Float64}})
        # restituisce un tempo e un indice {0,1,...,L} del vettore W 
        # C è cumsum(W) pre-allocato e aggiornato
        total = C[end]

        return (-log(rand()) / total), searchsortedfirst(C, rand() * total)
    end

    function updatesys!(ν,evento,L)
        if evento == 0  # l'evento 0 è quello che immette nel sistema una particella
            @inbounds ν[1] = true
        elseif evento == L      # l'evento L è quello che espelle dal sistema una particella
            @inbounds ν[L] = false
        else #casi evento ∈ {1...L-1}
            @inbounds ν[evento], ν[evento+1] = false, true  # sposto la particella dal sito n a n+1
        end
    end

    function nn(L)
        # calcola i gli indici dei processi da aggiornare
        [collect(max(0, e - 1):min(L, e + 1)) for e in 0:L]
    end

    update1rate!(i, ν, q) = q[i] * ν[i] * (1 - ν[i+1])

    function updaterates_nn!(W,processi_nn_selezionati, ν, q)
        @inbounds for i in processi_nn_selezionati   # +1 perche array Julia parte da 1
            W[i] = update1rate!(i, ν, q)
        end
    end

    function updatecumulative!(C, W, evento)
        @inbounds begin
            if evento == 0
                C[0] = W[0]   # base del cumulativo
                i0 = 1
            else
                i0 = max(0, evento - 1)
            end

            for i in i0:lastindex(W)
                C[i] = C[i-1] + W[i]
            end

        end
    end
##########################################################################################
         # GILLESPIE temporale
###############################################################################################
# 1 simulazione di gillespie - OTTIMIZZATA
    function Gillespie_temp(tasep, cond_init; termalization_time=10^5, n_measures=5*10^7, processi_nn=nn(tasep.L))
        # inizializzazione
        t = 0.0

        # inizializzazione sistema
        ν, W, q = Tasep2Sys(tasep, cond_init[1])

        # inizializzazione sistema
        sum_d = similar(ν) # per accumulare la somma delle densità
        fill!(sum_d, 0)

        current_accumulated = zeros(tasep.L) # per accumulare la somma delle correnti

        # PRE-ALLOCA cumsum una sola volta per simulazione
        # (non per ogni iterazione - questo è il grosso guadagno!)
        C_undef = Vector{Float64}(undef, tasep.L + 1)
        C = OffsetArray(C_undef, 0:tasep.L)

        # Calcola cumsum iniziale
        cumsum!(C, W)

        # scarto
        while t < termalization_time
            τ, evento = onestep_gillespie(C)
            t = t + τ
            # update sys
            updatesys!(ν, evento, tasep.L)
            # update rates
            updaterates_nn!(W, processi_nn[evento+1], ν, q)
            #cumsum!(C, W)  <-------vecchia versione
            updatecumulative!(C, W, evento) # aggiorno solo parte degli elementi della cumulativa, da evento in poi
        end

        # misuro
        iter = 0
        T = 0.0
        while iter < n_measures
            τ, evento = onestep_gillespie(C)
            T = T + τ
            iter += 1
            # update sys
            updatesys!(ν, evento, tasep.L)
            # update rates
            updaterates_nn!(W, processi_nn[evento+1], ν, q)
            updatecumulative!(C, W, evento) # aggiorno solo parte degli elementi della cumulativa, da evento in poi

            # somma delle DENSITA'
            sum_d = sum_d + τ .* ν

            # somma delle CORRENTI
            @inbounds @simd for i in 1:tasep.L
                current_accumulated[i] += tasep.q[i] * ν[i] * (1 - ν[i+1]) * τ
            end

        end

        d = sum_d ./ T
        c = current_accumulated ./ T

        return d[1:end-1], c[1:end-1]
    end
##########################################################################################
         # GILLESPIE P10
###########################################################################################
    function Gillespie_P10_temp(tasep, cond_init; termalization_time=2500, n_measures=1000, processi_nn=nn(tasep.L))
        # inizializzazione
        t = 0.0

        # inizializzazione sistema
        ν, W, q = Tasep2Sys(tasep, cond_init[1])

        # inizializzazione sistema
        sum_d = similar(ν) # per accumulare la somma delle densità
        fill!(sum_d, 0)

        P10_accumulated = zeros(tasep.L) # per accumulare la somma delle correnti

        # PRE-ALLOCA cumsum una sola volta per simulazione
        # (non per ogni iterazione - questo è il grosso guadagno!)
        C_undef = Vector{Float64}(undef, tasep.L + 1)
        C = OffsetArray(C_undef, 0:tasep.L)

        # Calcola cumsum iniziale
        cumsum!(C, W)

        # scarto
        while t < termalization_time
            τ, evento = onestep_gillespie(C)
            t = t + τ
            # update sys
            updatesys!(ν, evento, tasep.L)
            # update rates
            updaterates_nn!(W, processi_nn[evento+1], ν, q)
            #cumsum!(C, W)  <-------vecchia versione
            updatecumulative!(C, W, evento) # aggiorno solo parte degli elementi della cumulativa, da evento in poi
        end

        # misuro
        iter = 0
        T = 0.0
        while iter < n_measures
            τ, evento = onestep_gillespie(C)
            T = T + τ
            iter += 1
            # update sys
            updatesys!(ν, evento, tasep.L)
            # update rates
            updaterates_nn!(W, processi_nn[evento+1], ν, q)
            updatecumulative!(C, W, evento) # aggiorno solo parte degli elementi della cumulativa, da evento in poi

            # somma delle DENSITA'
            sum_d = sum_d + τ .* ν

            # somma delle CORRENTI
            @inbounds @simd for i in 1:tasep.L
                P10_accumulated[i] +=  ν[i] * (1 - ν[i+1]) * τ
            end

        end

        d = sum_d ./ T
        P10 = P10_accumulated ./ T

        return d[1:end-1], P10[1:end-1]
    end

