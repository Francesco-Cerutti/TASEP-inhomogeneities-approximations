using  OffsetArrays, Plots, LinearAlgebra, Statistics, Revise, ProgressMeter,Roots

################################################################################
                        # AUXILIARY DEFINITIONS
################################################################################
    function info(X)
        for (i,x) in enumerate(X)
            println( "$i has size = $(size(x)) and type = $(typeof(x))")
        end
    end

    norm0(u, v) = maximum(abs.(u .- v))


################################################################################
                        # TASEP DEFINITION
################################################################################
    struct Tasep
        L::Int
        α::Float64
        β::Float64
        q::OffsetVector{Float64, Vector{Float64}} 
    end

    function Tasep(L::Int, α::Float64, β::Float64, q_internal::Vector{Float64})
        @assert length(q_internal) == L-1 "The number of internal rates must be L-1"

        q = OffsetArray(vcat(1,q_internal,1), 0:L)
        return Tasep(L, α, β, q)
    end

    function lattice(L::Int64 # number of sites of the whole system
        ,q_defect1::Real      # hopping-rate value of the first defect
        ,l1::Int64             # length in bonds of the 1st defect
        ,from::Symbol;         # symbol selecting the case
        dbc::Int64=1          # position of the first defect from the reference boundary
                            # number of sites from the boundary to the defect
                            # if dbc = 1 the defect is a single site away from the boundary
        ,n_defects::Int64=1
        ,q_defect2::Real=1.0
        ,l2::Int64=0           # length in bonds of the 2nd defect
        ,dl::Int64=0)          # distance in sites between the 2 defects

    # define the vector of hopping rates
        q_internal = ones(L-1)
    # input checks
        l1 > 0 || throw(ArgumentError("l1 must be > 0"))
        l2 >= 0 || throw(ArgumentError("l2 must be >= 0"))
        dl >= 0 || throw(ArgumentError("dl must be >= 0"))
        n_defects in (1, 2) || throw(ArgumentError("n_defects must be 1 or 2"))
        from in (:center, :left, :right, :random) || throw(ArgumentError("invalid from"))

    # cases
        if (from == :center) && (n_defects == 1)
            center_idx = floor(Int, L / 2)
            q_internal[center_idx:center_idx+l1-1] .= q_defect1

        elseif  from == :left
            dbc >= 1 || throw(ArgumentError("Invalid configuration: dbc=$dbc must be >= 1 for from=:left"))
            dbc + l1 - 1 <= L - 1 || throw(ArgumentError("Invalid configuration: l1=$l1 is too large for dbc=$dbc. Reduce l1 or dbc"))
            q_internal[dbc:dbc+l1-1] .= q_defect1

        elseif from == :right
            # Defect on the right: ends at index (L-1) - dbc + 1, length l1
            end_idx = (L - 1) - dbc + 1
            start_idx = end_idx - l1 + 1
            start_idx >= 1 || throw(ArgumentError("Invalid configuration: l1=$l1 is too large for dbc=$dbc. Reduce l1 or dbc"))
            q_internal[start_idx:end_idx] .= q_defect1

        elseif (from == :center) && (n_defects == 2)
            l2 > 0 || throw(ArgumentError("if n_defects=2 then l2 must be > 0"))
            startl1, endl1, startl2, endl2 = centered_positions(q_internal, q_defect1, l1, q_defect2, l2, dl)

            q_internal[startl1:endl1] .= q_defect1
            q_internal[startl2:endl2] .= q_defect2

        elseif from == :random
            for i in 1:(L-1)
                q_internal[i] = rand()
            end
        end

        return q_internal
    end

    function centered_positions(q::AbstractVector{<:Real}, q_defect1::Real, l1::Int, q_defect2::Real, l2::Int, dl::Int)
        # dl = number of bonds between the two defective blocks
        L = length(q)
        total_length = l1 + dl + l2
        total_length > L && error("Invalid configuration: l1 + dl + l2 = $total_length > $L")

        startl1 = fld(L - total_length, 2) + 1
        endl1 = startl1 + l1 - 1
        startl2 = endl1 + dl + 1
        endl2 = startl2 + l2 - 1

        q[startl1:endl1] .= q_defect1
        q[startl2:endl2] .= q_defect2

        return startl1, endl1, startl2, endl2
    end

    function plot_lattice(q)
        L = length(q) + 1
        for i in 1:L-1
            val = "|"
            if q[i] != 1
                val = ""
            end
            println("   □: $i")
            println("$i: $(val)")
        end
        println("   □: $L")
    end




################################################################################
                        # PAIR APPROXIMATION DEFINITION
################################################################################
        
    function EuleroPA(tasep::Tasep,cond_init::Tuple{Float64,Float64};
                    Δt = 0.01,Nt = 1e6, tol = 1e-12, verbose = false, only_correlations = false)  

        ρ, P10 = set_cond_init_homo(tasep,cond_init[1],cond_init[2])

        ρ_next = copy(ρ)
        P10_next = copy(P10)
        density_variation = zeros(tasep.L)
        current_variation = zeros(tasep.L)

        for nt in 1:Nt
            # currents at boundaries + bulk
            for i in 0:tasep.L
                if i == 0 
                    P10_next[i] = tasep.α*(1-ρ[1])
                elseif i == tasep.L
                    P10_next[i] = tasep.β*ρ[tasep.L]
                else
                    current_variation[i] = (Gain(i,tasep.q,ρ,P10)-Loss(i,tasep.q,P10))
                    P10_next[i] = P10[i] + Δt*current_variation[i]
                end
            end
            
            # Euler update
            for i in 1:tasep.L
                density_variation[i] = (tasep.q[i-1]*P10[i-1] - tasep.q[i]*P10[i])
                ρ_next[i] = ρ[i] + Δt*density_variation[i]
            end
            # update everything
            P10 .= P10_next
            ρ .= ρ_next

            # tolerance checks
            actual_density_tol = maximum(abs.(density_variation))
            actual_current_tol = maximum(abs.(current_variation))

            if (actual_density_tol < tol) && (actual_current_tol < tol)
                verbose ? println("Num_it= $nt , density tol reached $(actual_density_tol), current tol reached $(actual_current_tol))") : nothing
                break
            end
        end
        J = tasep.q.*P10
        (only_correlations) ? (return P10) : (return ρ[1:end-1], J[1:end-1]) 
    end

    function EuleroPA_P10(tasep::Tasep,cond_init::Tuple{Float64,Float64};
                    Δt = 0.01,Nt = 1e6, tol = 1e-12, verbose = false)  

        ρ, P10 = set_cond_init_homo(tasep,cond_init[1],cond_init[2])

        ρ_next = copy(ρ)
        P10_next = copy(P10)
        density_variation = zeros(tasep.L)
        current_variation = zeros(tasep.L)

        for nt in 1:Nt
            # currents at boundaries + bulk
            for i in 0:tasep.L
                if i == 0 
                    P10_next[i] = tasep.α*(1-ρ[1])
                elseif i == tasep.L
                    P10_next[i] = tasep.β*ρ[tasep.L]
                else
                    current_variation[i] = (Gain(i,tasep.q,ρ,P10)-Loss(i,tasep.q,P10))
                    P10_next[i] = P10[i] + Δt*current_variation[i]
                end
            end
            
            # Euler update
            for i in 1:tasep.L
                density_variation[i] = (tasep.q[i-1]*P10[i-1] - tasep.q[i]*P10[i])
                ρ_next[i] = ρ[i] + Δt*density_variation[i]
            end
            # update everything
            P10 .= P10_next
            ρ .= ρ_next

            # tolerance checks
            actual_density_tol = maximum(abs.(density_variation))
            actual_current_tol = maximum(abs.(current_variation))

            if (actual_density_tol < tol) && (actual_current_tol < tol)
                verbose ? println("Num_it= $nt , density tol reached $(actual_density_tol), current tol reached $(actual_current_tol))") : nothing
                break
            end
        end

        return ρ, P10
    end

    #=   OLD VERSION
    function EuleroPA(tasep::Tasep,cond_init::Tuple{Float64,Float64},Nt,T) 

        ρ, P10 = set_cond_init_homo(tasep,cond_init[1],cond_init[2])

        Δt = T/Nt

        ρ_next = copy(ρ)
        P10_next = copy(P10)


        @showprogress for nt in 1:Nt
            # currents at boundaries + bulk
            for i in 0:tasep.L
                if i == 0 
                    P10_next[i] = tasep.α*(1-ρ[1])
                elseif i == tasep.L
                    P10_next[i] = tasep.β*ρ[tasep.L]
                else
                    P10_next[i] = P10[i] + Δt*(Gain(i,tasep.q,ρ,P10)-Loss(i,tasep.q,P10))
                end
            end
            
            # Euler update
            for i in 1:tasep.L
                ρ_next[i] = ρ[i] + Δt*(P10[i-1] - P10[i])
            end
            # update everything
            P10 .= P10_next
            ρ .= ρ_next
        end

        return ρ, P10
    end
    =#


    Gain(i,q,ρ,P10) = G1(q[i-1], ρ[i], ρ[i+1], P10[i-1], P10[i]) + G2(q[i+1], ρ[i], ρ[i+1], P10[i], P10[i+1])
    G1(qm1,ρ,ρp1,P10m1,P10) = qm1*(P10m1*(1-ρp1-P10))/(1-ρ) 
    G2(qp1,ρ,ρp1,P10,P10p1) = qp1*(P10p1*(ρ-P10))/ρp1
    Loss(i,q,P10) =  q[i]*P10[i]


    #=
    function EuleroPA(tasep::Tasep,cond_init::Tuple{Float64,Float64},Nt,T) 

        ρ, P10 = set_cond_init_homo(tasep,cond_init[1],cond_init[2])
        Δt = T/Nt
        q = tasep.q

        ρ_next = copy(ρ)
        P10_next = copy(P10)

        @showprogress for nt in 1:Nt
            # currents at boundaries + bulk
            for i in 0:tasep.L
                if i == 0 
                    P10_next[i] = tasep.α*(1-ρ_next[1])
                elseif i == tasep.L
                    P10_next[i] = tasep.β*ρ_next[tasep.L]
                else
                    P10_next[i] = P10[i] + Δt*(Gain(q[i-1],q[i+1],ρ_next[i],ρ_next[i+1],P10[i-1],P10[i],P10[i+1])-Loss(q[i],P10[i]))
                end

            # update P10 pointwise since it will be reused
            P10[i] = P10_next[i]
            end


            # Euler update
            for i in 1:tasep.L
                ρ_next[i] = ρ[i] + Δt*(P10[i-1] - P10[i])
            end
            # update everything
            ρ .= ρ_next
        end

        return ρ, P10
    end
    =#

    function set_cond_init_homo(tasep::Tasep, ρ0::S, P10_0::S) where S <:Real

        ρ_homo = OffsetArray(vcat(tasep.α, fill(ρ0, tasep.L), 1-tasep.β), 0:tasep.L+1)   
        P10_homo = OffsetArray(vcat(tasep.α*(1-ρ_homo[1]), fill(P10_0, tasep.L-1), tasep.β*ρ_homo[tasep.L]), 0:tasep.L)   
        
        return ρ_homo,P10_homo
    end



################################################################################
                        # TRIPLET APPROXIMATION DEFINITION
################################################################################
    function set_cond_init_homo(tasep::Tasep, ρ0::S, P10_0::S, P100_0::S, P110_0::S) where S<:Real
        # sets an initial value for all elements of the vector

        ρ_homo =  fill(ρ0, tasep.L)
        P10_homo = fill(P10_0, tasep.L - 1)
        P100_homo = fill(P100_0, tasep.L - 2)
        P110_homo =fill(P110_0, tasep.L - 2)

        return ρ_homo, P10_homo, P100_homo, P110_homo

    end

    # derived unknowns
        function P000_(P000, L, ρ, P10, P100)
            @inbounds @simd for i in 1:(L-2)
                P000[i] = (1.0 - ρ[i+2]) - P10[i+1] - P100[i]
            end
            return P000
        end

        function P101_(P101, L, P10, P100)
            @inbounds @simd for i in 1:(L-2)
                P101[i] = P10[i] - P100[i]
            end
            return P101
        end

        function P01_(P01, L, ρ, P10)
            @inbounds @simd for i in 1:(L-1)
                P01[i] = (1.0 - ρ[i]) - (1.0 - ρ[i+1]) + P10[i]   # = -ρ[i] + ρ[i+1] + P10[i]
            end
            return P01
        end

        function P00_(P00, L, ρ, P10)
            @inbounds @simd for i in 1:(L-1)
                P00[i] = (1.0 - ρ[i+1]) - P10[i]
            end
            return P00
        end

        function P010_(P010, L, P10, P110)
            @inbounds @simd for i in 1:L-2
                P010[i] = P10[i+1] - P110[i]
            end
            return P010
        end

        function P111_(P111, L, ρ, P10, P110)
            @inbounds @simd for i in 1:L-2
                P111[i] = ρ[i] - P10[i] - P110[i]
            end
            return P111
        end

        function P11_(P11, L, ρ, P10)
            @inbounds @simd for i in 1:(L-1)
                P11[i] = ρ[i] - P10[i]
            end
            return P11
        end



    # incognite
        function P100_(P100_next, P100_variation, L, Δt, α, β, q, P00, P01, P000, P010, P100, P101, P110) # ordered from smallest to largest, P000 > P10 since it has more digits

            # edge cases
            P100_variation[1] = α * P000[1] + q[3] * P101[1] * P010[2] / P01[2] - q[1] * P100[1]
            P100_next[1] = P100[1] + Δt * P100_variation[1]
            # general case 
            @inbounds @simd for i in 2:L-3
                P100_variation[i] = q[i-1] * P100[i-1] * P000[i] / P00[i] + q[i+2] * P101[i] * P010[i+1] / P01[i+1] - q[i] * P100[i]
                P100_next[i] = P100[i] + Δt * P100_variation[i]
            end
            P100_variation[L-2] = q[L-3] * P100[L-3] * P000[L-2] / P00[L-2] + β * P101[L-2] - q[L-2] * P100[L-2]
            P100_next[L-2] = P100[L-2] + Δt * P100_variation[L-2]

            return P100_next, P100_variation
        end

        function P110_(P110_next, P110_variation, L, Δt, α, β, q, P01, P11, P010, P101, P110, P111) # ordered from smallest to largest, P000 > P10 since it has more digits

            P110_variation[1] = α * P010[1] + q[3] * P111[1] * P110[2] / P11[2] - q[2] * P110[1]
            P110_next[1] = P110[1] + Δt * P110_variation[1]

            # general case 
            @inbounds @simd for i in 2:L-3
                P110_variation[i] = q[i-1] * P101[i-1] * P010[i] / P01[i] + q[i+2] * P111[i] * P110[i+1] / P11[i+1] - q[i+1] * P110[i]
                P110_next[i] = P110[i] + Δt * P110_variation[i]
            end

            P110_variation[L-2] = q[L-3] * P101[L-3] * P010[L-2] / P01[L-2] + β * P111[L-2] - q[L-1] * P110[L-2]
            P110_next[L-2] = P110[L-2] + Δt * P110_variation[L-2]

            return P110_next, P110_variation
        end

        function P10_(P10_next, P10_variation, L, Δt, α, β, q, P00, P10, P11, P100, P110) # ordered from smallest to largest, P000 > P10 since it has more digits

            P10_variation[1] = α * P00[1] + q[2] * P110[1] - q[1] * P10[1]
            P10_next[1] = P10[1] + Δt * P10_variation[1]
            # general case 
            @inbounds @simd  for i in 2:L-2
                P10_variation[i] = q[i-1] * P100[i-1] + q[i+1] * P110[i] - q[i] * P10[i]
                P10_next[i] = P10[i] + Δt * P10_variation[i]
            end
            P10_variation[L-1] = q[L-2] * P100[L-2] + β * P11[L-1] - q[L-1] * P10[L-1]
            P10_next[L-1] = P10[L-1] + Δt * P10_variation[L-1]

            return P10_next, P10_variation

        end

        function ρ_(ρ_next, density_variation, L, Δt, α, β, q, ρ, P10)

            density_variation[1] = α * (1 - ρ[1]) - q[1] * P10[1]
            ρ_next[1] = ρ[1] + Δt * density_variation[1]

            @inbounds @simd for i in 2:L-1
                density_variation[i] = q[i-1] * P10[i-1] - q[i] * P10[i]
                ρ_next[i] = ρ[i] + Δt * density_variation[i]
            end

            density_variation[L] = q[L-1] * P10[L-1] - β * ρ[L]
            ρ_next[L] = ρ[L] + Δt * density_variation[L]

            return ρ_next, density_variation
        end

    # triplet function
        function EuleroTRI(tasep::Tasep,cond_init::Tuple{Float64,Float64,Float64,Float64};
                    Δt = 0.01,Nt = 1e6, tol = 1e-12, verbose = false)  

            ρ, P10, P100, P110 = set_cond_init_homo(tasep,cond_init[1],cond_init[2],cond_init[3],cond_init[4])

            # init
                α = tasep.α
                β = tasep.β
                q = tasep.q
                L = tasep.L
                ρ_next = copy(ρ)
                P10_next = copy(P10)
                P100_next = copy(P100)
                P110_next = copy(P110)
                
                P000 = zeros(L-2)
                P101 = zeros(L-2)
                P00 = zeros(L-1)
                P01 = zeros(L-1)
                P11 = zeros(Float64, L-1)
                P111 = zeros(Float64, L-2)
                P010 = zeros(Float64, L-2)

                ρ_variation = zeros(tasep.L)
                P10_variation = zeros(tasep.L)
                P100_variation = zeros(tasep.L)
                P110_variation = zeros(tasep.L)

            @inbounds for nt in 1:Nt
                    
                    P000 = P000_(P000, L, ρ, P10, P100)
                    P101 = P101_(P101, L, P10, P100)
                    P01  = P01_(P01, L, ρ, P10)
                    P00  = P00_(P00, L, ρ, P10)
                    P010 = P010_(P010, L, P10, P110)
                    P111 = P111_(P111, L, ρ, P10, P110)
                    P11  = P11_(P11, L, ρ, P10)

                    P100_next, P100_variation = P100_(P100_next, P100_variation ,L, Δt,α, β,q, P00, P01, P000, P010, P100, P101, P110)
                
                    P110_next, P110_variation = P110_(P110_next, P110_variation,L, Δt,α, β,q, P01, P11, P010, P101, P110, P111)
            
                    P10_next, P10_variation= P10_(P10_next, P10_variation,L, Δt,α, β, q, P00, P10, P11, P100, P110)
                    
                    ρ_next, ρ_variation= ρ_(ρ_next, ρ_variation, L, Δt, α, β, q, ρ, P10)
                    
                    # update everything
                    ρ .= ρ_next
                    P10 .= P10_next
                    P100 .= P100_next
                    P110 .= P110_next

                    # tolerance checks
                    actual_ρ_tol = maximum(abs.(ρ_variation))
                    actual_P10_tol = maximum(abs.(P10_variation))
                    actual_P100_tol = maximum(abs.(P100_variation))
                    actual_P110_tol = maximum(abs.(P110_variation))


                    if (actual_ρ_tol < tol) && (actual_P10_tol < tol) && (actual_P100_tol < tol) && (actual_P110_tol < tol)
                        verbose ? println("Num_it= $nt , tol reached density:  $(actual_density_tol), current: $(actual_current_tol),P100: $(actual_P100_tol),P110: $(actual_P110_tol)") : nothing
                        break
                    end
                end

                return ρ, q[1:L-1] .* P10, P100, P110
        end

        function EuleroTRI_P10(tasep::Tasep,cond_init::Tuple{Float64,Float64,Float64,Float64};
                    Δt = 0.01,Nt = 1e6, tol = 1e-12, verbose = false)  

            ρ, P10, P100, P110 = set_cond_init_homo(tasep,cond_init[1],cond_init[2],cond_init[3],cond_init[4])

            # init
                α = tasep.α
                β = tasep.β
                q = tasep.q
                L = tasep.L
                ρ_next = copy(ρ)
                P10_next = copy(P10)
                P100_next = copy(P100)
                P110_next = copy(P110)
                
                P000 = zeros(L-2)
                P101 = zeros(L-2)
                P00 = zeros(L-1)
                P01 = zeros(L-1)
                P11 = zeros(Float64, L-1)
                P111 = zeros(Float64, L-2)
                P010 = zeros(Float64, L-2)

                ρ_variation = zeros(tasep.L)
                P10_variation = zeros(tasep.L)
                P100_variation = zeros(tasep.L)
                P110_variation = zeros(tasep.L)

            @inbounds for nt in 1:Nt
                    
                    P000 = P000_(P000, L, ρ, P10, P100)
                    P101 = P101_(P101, L, P10, P100)
                    P01  = P01_(P01, L, ρ, P10)
                    P00  = P00_(P00, L, ρ, P10)
                    P010 = P010_(P010, L, P10, P110)
                    P111 = P111_(P111, L, ρ, P10, P110)
                    P11  = P11_(P11, L, ρ, P10)

                    P100_next, P100_variation = P100_(P100_next, P100_variation ,L, Δt,α, β,q, P00, P01, P000, P010, P100, P101, P110)
                
                    P110_next, P110_variation = P110_(P110_next, P110_variation,L, Δt,α, β,q, P01, P11, P010, P101, P110, P111)
            
                    P10_next, P10_variation= P10_(P10_next, P10_variation,L, Δt,α, β, q, P00, P10, P11, P100, P110)
                    
                    ρ_next, ρ_variation= ρ_(ρ_next, ρ_variation, L, Δt, α, β, q, ρ, P10)
                    
                    # update everything
                    ρ .= ρ_next
                    P10 .= P10_next
                    P100 .= P100_next
                    P110 .= P110_next

                    # tolerance checks
                    actual_ρ_tol = maximum(abs.(ρ_variation))
                    actual_P10_tol = maximum(abs.(P10_variation))
                    actual_P100_tol = maximum(abs.(P100_variation))
                    actual_P110_tol = maximum(abs.(P110_variation))


                    if (actual_ρ_tol < tol) && (actual_P10_tol < tol) && (actual_P100_tol < tol) && (actual_P110_tol < tol)
                        verbose ? println("Num_it= $nt , tol reached density:  $(actual_density_tol), current: $(actual_current_tol),P100: $(actual_P100_tol),P110: $(actual_P110_tol)") : nothing
                        break
                    end
                end

                return ρ, P10, P100, P110
        end

################################################################################
                        # MEAN FIELD APPROXIMATION DEFINITION
################################################################################
   

    function EuleroMF(tasep::Tasep,cond_init::Tuple{Float64,Float64};
                    Δt = 0.01,Nt = 1_000_000, tol = 1e-12, verbose = false) 

        ρ, P10 = set_cond_init_homo(tasep,cond_init[1],cond_init[2])
        ρ_next = copy(ρ)
        variation = zeros(tasep.L)

        for nt in 1:Nt

        # currents
            for i in 0:tasep.L
                P10[i] = ρ[i]*(1 - ρ[i+1])
            end
        # density
            for i in 1:tasep.L
                variation[i] = tasep.q[i-1]*P10[i-1] - tasep.q[i]*P10[i]
                ρ_next[i] = ρ[i] + Δt*variation[i]
            end
        # update
            ρ .= ρ_next

        # tolerance checks
            actual_tol = maximum(abs.(variation))

            if actual_tol < tol
                verbose ? println("Convergence reached in $nt iterations = $(actual_tol))") : nothing
                break
            end

        end
        J = tasep.q.*P10
        return ρ[1:end-1], J[1:end-1]
    end


    #= OLD VERSION
        function EuleroMF(tasep::Tasep,cond_init::Tuple{Float64,Float64},Nt,T) 

        ρ, P10 = set_cond_init_homo(tasep,cond_init[1],cond_init[2])
        Δt = T/Nt

        ρ_next = copy(ρ)

        @showprogress for nt in 1:Nt
            
            for i in 0:tasep.L
                P10[i] = ρ[i]*(1 - ρ[i+1])
            end
        
        
            for i in 1:tasep.L
                ρ_next[i] = ρ[i] + Δt*(tasep.q[i-1]*P10[i-1] - tasep.q[i]*P10[i])
            end
            


            ρ .= ρ_next
        end

        return ρ, tasep.q.*P10
        end
    =#



    function set_cond_init(tasep::Tasep, ρ0::S, P10_0::S) where S <:Real

        ρ_homo = OffsetArray(vcat(tasep.α, fill(ρ0, tasep.L), 1-tasep.β), 0:tasep.L+1)   
        P10_homo = OffsetArray(vcat(tasep.α*(1-ρ_homo[1]), fill(P10_0, tasep.L-1), tasep.β*ρ_homo[tasep.L]), 0:tasep.L)   
        
        return ρ_homo,P10_homo
    end

    function stationaryMF(tasep::Tasep,cond_init::Tuple{Float64,Float64};
        N_iter = 10_000, tol = 1e-12, verbose = false) 
        
        ρ, P10 = set_cond_init(tasep, cond_init[1], cond_init[2]) # initialize ρ and P10 with homogeneous conditions
        J = copy(P10)

        for i in 1:N_iter
            ρ_previous = copy(ρ)
            for l in 1:tasep.L
                ρ[l] = (tasep.q[l-1] * ρ[l-1])/(tasep.q[l-1] * ρ[l-1] + tasep.q[l] * (1 - ρ[l+1]))
            end

            actual_tol = maximum(abs.(ρ - ρ_previous))
            if actual_tol < tol
                verbose ? println("Convergence reached in $i iterations = $(actual_tol))") : nothing
                break
            end
        end

        for i in 0:tasep.L
            J[i] = tasep.q[i]*ρ[i]*(1 - ρ[i+1])
        end
    
        return ρ, J

    end
    #= OLD VERSION
    function stationaryMF(tasep::Tasep,cond_init::Tuple{Float64,Float64};
        N_iter = 10_000, tol = 1e-12) 
        
        ρ, P10 = set_cond_init(tasep, cond_init[1], cond_init[2]) # initialize ρ and P10 with homogeneous conditions
        J = copy(P10)

        @showprogress for i in 1:N_iter
            ρ_previous = copy(ρ)
            for l in 1:tasep.L
                ρ[l] = (tasep.q[l-1] * ρ[l-1])/(tasep.q[l-1] * ρ[l-1] + tasep.q[l] * (1 - ρ[l+1]))
            end

            actual_tol = maximum(abs.(ρ - ρ_previous))
            if actual_tol < tol
                println("Convergence reached in $i iterations = $(actual_tol))")
                break
            end
        end

        for i in 0:tasep.L
            J[i] = tasep.q[i]*ρ[i]*(1 - ρ[i+1])
        end
    
        return ρ, J

    end
    =#

################################################################################
                        # ISA DEFINITION
################################################################################
    
 
    # single-defect system
        α_star(q) = (2+3*q)/4 - sqrt(((2+3*q)^2)/16 - q);

        function Jisa1(tasep::Tasep)
            
            q = only(Float64.(tasep.q[floor(Int, tasep.L/2)])) # assume there is a single defect at the central site

            return α_star(q)*(1-α_star(q))
        end

        Jisa1(q) = α_star(q)*(1-α_star(q))

    # system with l defects

    function Z(α, β, N; tol=1e-2) # N number of sites of the slow subsystem

        total = 0.0

        if abs(α-β)<tol
            for j in 1:N
                total += j*(factorial(2*N-1-j))/(factorial(N-j)*factorial(N))*(j+1)*(1/α)^j
            end
        else
            for j in 1:N
                total += j*(factorial(2*N-1-j))/(factorial(N-j)*factorial(N))*(((1/β)^(j+1) - (1/α)^(j+1))/(1/β - 1/α))
            end
        end

        return total
    end

    JisaN(α, β, N) = Z(α, β, N-1)/Z(α, β, N)


    function α_star(q, l) # l number of defective bonds

        f(α) = α*(1-α) - q*JisaN((1-α)/q, (1-α)/q, l+1) # l+1 because the slow subsystem has l+1 sites

        return find_zero(f, 2/5)
    end
    JisaN(q, l) = α_star(q, l)*(1-α_star(q, l))



################################################################################
                        # CORRELATIONS DEFINITION
################################################################################
    function corr_Gillespie(tasep, cond_init; n_simu=n_simu, total_time=total_time, processi_nn=nn(tasep.L), terms=false)

        ρ, P10 = Gillespie_P10(tasep, cond_init; n_simu=n_simu, total_time=total_time, processi_nn=nn(tasep.L))
        # product of mean values
        prod = [ρ[i] * (1 - ρ[i+1]) for i in 1:length(ρ)-1]

        # do I want the two terms of the correlation?
        terms ? (return ρ, P10, prod, P10 .- prod) : (return P10 .- prod)
    end

    function corr_PA(tasep, cond_init; terms=false, Δt=1, Nt=1e6)
        ρ, P10 = EuleroPA_P10(tasep, cond_init; Δt=Δt, Nt=Nt)
        prod = [ρ[i] * (1 - ρ[i+1]) for i in firstindex(ρ):(lastindex(ρ)-1)]


        # remove the first and last element, the 2 reservoirs
        P10_short = P10[firstindex(P10)+1:lastindex(P10)-1]
        prod_short = prod[firstindex(prod)+1:lastindex(prod)-1]
        ρ_short = ρ[firstindex(ρ)+1:lastindex(ρ)-1]

        # do I want the two terms of the correlation?
        terms ? (return ρ_short, P10_short, prod_short, P10_short .- prod_short) : (return P10_short .- prod_short)
    end

################################################################################
                        # CURRENT AS A FUNCTION OF q DEFINITION
################################################################################
    function current_q(qmin, qmax, nq;
        methods=["MF", "PA", "stat", "ISA", "GILL", "TRI"],
        L=200, α=0.5, β=0.5,
        Nt=1e6, Δt=0.1, tol=1e-10,
        defect1_length=1,
        from::Symbol=:center,
        distance_from_boundary::Int64=1, n_defects::Int64=1,
        q_defect2::Real=1.0, defect2_length::Int64=0,
        distance_between_defects::Int64=1,
        termalization_time=1e3, n_measures=1e7
        )

        q_values = range(qmin, qmax, length=nq)
        currents = Dict{String,Vector{Float64}}()

        for method in methods
            currents[method] = Float64[]
        end

        begin
            @inbounds @showprogress for q_defect in q_values

                # define the lattice 
                q = lattice(L, q_defect, defect1_length, from; dbc=distance_from_boundary, n_defects=n_defects, q_defect2=q_defect2, l2=defect2_length, dl=distance_between_defects)
                tasep = Tasep(L, α, β, q)  # example tasep with variable q
                cond_init = (0.5, 0.0)  # example initial conditions

                if "MF" in methods
                    ρ, J = EuleroMF(tasep, cond_init; Nt=Nt, Δt=Δt, tol=tol)
                    push!(currents["MF"], mean(J))
                end

                if "PA" in methods
                    ρ, J = EuleroPA(tasep, cond_init; Nt=Nt, Δt=Δt, tol=tol)
                    push!(currents["PA"], mean(J))
                end

                if "stat" in methods
                    ρ, J = stationaryMF(tasep, cond_init; verbose=false)
                    push!(currents["stat"], mean(J))
                end

                if "ISA" in methods
                    J = JisaN(q_defect, defect1_length)
                    push!(currents["ISA"], J)
                end
                if "GILL" in methods
                    _, J = Gillespie_temp(tasep, cond_init; n_measures=n_measures, termalization_time=termalization_time, processi_nn=nn(tasep.L))
                    push!(currents["GILL"], mean(J))
                end
                if "TRI" in methods
                    ρ, J, _, _ = EuleroTRI(tasep, (0.5, 0.25, 0.25, 0.25); Nt=Nt, Δt=Δt, tol=tol)
                    push!(currents["TRI"], mean(J))
                end

            end
        end
        return q_values, currents
    end



nothing