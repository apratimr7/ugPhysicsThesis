# PERFECTION 
using Distributed
@everywhere using DifferentialEquations
@everywhere using LinearAlgebra
using SharedArrays
using Plots
using Printf
using LaTeXStrings
using Colors


# Add worker processes if not already added
if nworkers() == 1
    addprocs(CPU_THREADS - 1)
end

@everywhere begin
     # Define parameters
     θ = 1.1
     α = 2.7
     γ = 0.55
     β = 1.4
     s = 10
     g₁ = 4
     g₂ = 4
     αₙₛ = 0.02
     λ = 0.6
     ω = 1.2
 
     # Spatial grid parameters
     Nx, Ny = 10,10
     Lx, Ly = 255, 255
     dx, dy = Lx/Nx, Ly/Ny
     x = range(-Lx, Lx, length=Nx)
     y = range(-Ly, Ly, length=Ny)
 
     potential_height = 44
     potential_period = 12.998

     V(x,y)= @. potential_height * sinc(potential_period * hypot(x,y))
     name = "Sinc_tcmsg_profile"

    # V(x,y) = @. sin(5* x*y/10^2) * exp(-hypot(x,y)^2/100)
    # name = "TCMSG"
   
    # using SpecialFunctions
    # V(x,y) = @. 5 * besselj(3, 1* hypot(x,y))*exp(-0.35*hypot(x,y)^2)
    # name = "Bessel-Gaussian"
    
    # Gaussian initial condition

    function gaussian_initial(x, y; x0=0, y0=0, σx=1, σy=1, A=1)
        return A * exp(-((x-x0)^2 / (σx^2) + (y-y0)^2 / (σy^2)))
    end
    
     # TCMSG initial condition
    function tcmsg_initial(x, y; x0=0, y0=0, σx=1, σy=1, A=1)
        return @.A*sin(5*(x-x0)*(y-y0)/σx^2)*(exp(-(x-x0)^2/100-(y-y0)^2/σy^2))
    end
     # 2D Laplacian using finite differences
     function laplacian(E)
         ∇²E = zeros(ComplexF64, Nx, Ny)
         for i in 2:Nx-1, j in 2:Ny-1
             ∇²E[i,j] = (E[i+1,j] + E[i-1,j] + E[i,j+1] + E[i,j-1] - 4E[i,j]) / (dx^2)
         end
         ∇²E[1,:], ∇²E[Nx,:] = ∇²E[Nx-1,:], ∇²E[2,:]
         ∇²E[:,1], ∇²E[:,Ny] = ∇²E[:,Ny-1], ∇²E[:,2]
         return ∇²E
     end
 
     # Define the complex Ginzburg-Landau equation
     function CGLE!(dE, E, p, t)
         μ, σ = p
         a = (σ * λ^2) / (λ^2 + ω^2)
         b = (σ * λ * ω) / (λ^2 + ω^2)
         ∇²E = laplacian(E)
         for i in 1:Nx, j in 1:Ny
             term1 = -(1 - im*θ) + (μ * (1 - im * α)) / (1 + g₁ * abs2(E[i,j]))
             term2 = -(γ * (1 - im * β)) / (1 + s* g₂ * abs2(E[i,j]))
             term3 = αₙₛ + a - im * b
             dE[i,j] = (term1 + term2 + term3) * E[i,j] + im * ∇²E[i,j] + im * V(x[i], y[j]) * E[i,j]
         end
     end
 

    # Function to compute steady state and return the full solution
    function simulate_system(p)
        E₀ = [complex(gaussian_initial(x[i], y[j], σx=10, σy=10), 0.25) for i in 1:Nx, j in 1:Ny]
        E₀ = [complex(tcmsg_initial(x[i], y[j], σx=10, σy=10), 0.25) for i in 1:Nx, j in 1:Ny]
        tspan = (0.0, 100.0)
        prob = ODEProblem(CGLE!, E₀, tspan, p)
        sol = solve(prob, save_everystep=true)
        return sol
    end

    # Compute the largest Lyapunov exponent
    function largest_lyapunov_exponent(p, num_iterations=1000)
        sol = simulate_system(p)
        E = sol[end]
        δ = 1e-10
        λ_sum = 0.0
        
        for _ in 1:num_iterations
            E_perturbed = E + δ * randn(ComplexF64, size(E))
            
            prob = ODEProblem(CGLE!, E, (0.0, 0.1), p)
            sol1 = solve(prob, save_everystep=false)
            
            prob_perturbed = ODEProblem(CGLE!, E_perturbed, (0.0, 0.1), p)
            sol2 = solve(prob_perturbed, save_everystep=false)
            
            d = norm(sol1[end] - sol2[end])
            λ_sum += log(d / δ) / 0.1
            
            E = sol1[end]
        end
        
        return λ_sum / num_iterations
    end
end

# Stability analysis
μ_range = 0.0:0.2:2
σ_range = 0.0:0.2:2
stability_matrix = SharedArray{Float64}(length(μ_range), length(σ_range))

@sync @distributed for i in 1:length(μ_range)
    for j in 1:length(σ_range)
        μ, σ = μ_range[i], σ_range[j]
        p = (μ, σ)
        lyap_exp = largest_lyapunov_exponent(p)
        stability_matrix[i, j] = lyap_exp
        @printf("μ = %.2f, σ = %.2f, Lyapunov Exponent = %.4f\n", μ, σ, lyap_exp)
    end
end

# Convert SharedArray to regular Array for plotting
stability_matrix = Array(stability_matrix)

# Plotting
# Use PyPlot backend for better quality

# 1. Heatmap of Lyapunov exponents
p1 = heatmap(μ_range, σ_range, stability_matrix', 
             xlabel=L"\mu", ylabel=L"\sigma", 
             title="Stability Analysis: Largest Lyapunov Exponent",
             c=:turbo, colorbar_title="Lyapunov Exponent",
             clims=(-1, 1), # Adjust color range
             aspect_ratio=:equal)
contour!(p1, μ_range, σ_range, stability_matrix', 
         levels=[0], color=:white, linewidth=2, alpha=0.7)

# 2. 3D surface plot of Lyapunov exponents
p2 = surface(μ_range, σ_range, stability_matrix', xlims=(0,2),ylims=(0,2),
             xlabel=L"\mu", ylabel=L"\sigma", zlabel="Lyapunov Exponent",
             title="3D Stability Landscape",
             c=:turbo, alpha=0.8)

# 3. Contour plot with filled contours
p3 = contour(μ_range, σ_range, stability_matrix', ylims=(0,10),
              levels=40, # Increase for smoother contours
              xlabel=L"\mu", ylabel=L"\sigma",
              title="Stability Regions",
              c=:coolwarm)
contour!(p3, μ_range, σ_range, stability_matrix', 
         levels=[0], color=:black, linewidth=2, alpha=0.7)

# 4. Slice plots
μ_slice_index = length(μ_range) ÷ 2
σ_slice_index = length(σ_range) ÷ 2

println(μ_slice_index," ",μ_range[μ_slice_index])
println(σ_slice_index," ",σ_range[σ_slice_index])

p4 = plot(σ_range, stability_matrix[μ_slice_index, :], 
          xlabel=L"\sigma", ylabel="Lyapunov Exponent",
          title=L"Slice at constant \mu = " * string(round(μ_range[μ_slice_index], digits=2)),
          legend=false, lw=2)
hline!(p4, [0], color=:red, linestyle=:dash)

p5 = plot(μ_range, stability_matrix[:, σ_slice_index], 
          xlabel=L"\mu", ylabel="Lyapunov Exponent",
          title=L"Slice at constant \sigma = " * string(round(σ_range[σ_slice_index], digits=2)),
          legend=false, lw=2)
hline!(p5, [0], color=:red, linestyle=:dash)


function plot_phase_trajectory(sol, i, j, title)
    trajectory = [sol[t][i, j] for t in 1:length(sol.t)]
    real_traj = real.(trajectory)
    imag_traj = imag.(trajectory)
    
    p = plot(real_traj, imag_traj, 
             xlabel="Re(E)", ylabel="Im(E)", 
             title=title, 
             legend=false, 
             aspect_ratio=:equal,
             linewidth=1.5,
             alpha=0.7)
    
    # Add arrows
    num_arrows = 10
    arrow_indices = round.(Int, range(1, length(trajectory), length=num_arrows))
    for k in arrow_indices[1:end-1]
        arrow_start = (real_traj[k], imag_traj[k])
        arrow_end = (real_traj[k+1], imag_traj[k+1])
        arrow_vec = arrow_end .- arrow_start
        quiver!([arrow_start[1]], [arrow_start[2]], 
                quiver=([arrow_vec[1]], [arrow_vec[2]]), 
                color=:red, 
                arrow=arrow(:closed, :head, 0.05, 0.03))
    end
    
    return p
end

# Select specific (μ, σ) values for phase trajectories
μ_values = [0.5, 1, 1.5,2]
σ_values = [0.5, 1, 1.5,2]

# Plot phase trajectories
phase_plots = []
for (idx, (μ, σ)) in enumerate(zip(μ_values, σ_values))
    sol = simulate_system((μ, σ))
    p = plot_phase_trajectory(sol, Nx÷2, Ny÷2, "Phase Trajectory (μ=$μ, σ=$σ)")
    push!(phase_plots, p)
end

# Combine all plots
plot(p1, p2, p3, p4, p5, phase_plots..., layout=(3,3), size=(1800,1800), link=:none)

# Save the combined plot
savefig("./imgs/CGLE With "*name*" potential Stability analysis.png")

# Save individual plots for detailed viewing
savefig(p1, "./imgs/"*name*"_lyapunov_heatmap.png")
savefig(p2, "./imgs/"*name*"_lyapunov_3d_surface.png")
savefig(p3, "./imgs/"*name*"_stability_regions_contour.png")
savefig(p4, "./imgs/"*name*"_mu_slice_plot.png")
savefig(p5, "./imgs/"*name*"_sigma_slice_plot.png")
for (idx, p) in enumerate(phase_plots)
    savefig(p, "./imgs/"*name*"_phase_trajectory_$idx.png")
end