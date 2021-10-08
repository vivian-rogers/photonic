push!(LOAD_PATH, "./src/")


using Plots
using PlotStuff
using Constants
using LinearAlgebra
using UsefulFunctions
using Operators
using GrapheneParameters
using Bands

function H(k)
	#in basis |1> |2> for the 2 atom unit cell
	H₀ = [ε t
	      t ε] #hamiltonian of the unit cell
	H₂₁ = [0 0
	       t 0] #so <2|1> has a hopping energy
	H₁₂ = [0 t
	       0 0] #overlap of <1|2>
	H = H₀ .+ exp(im*k⋅r₁)*H₂₁ .+ exp(im*k⋅r₂)*H₂₁ .+ exp(-im*k⋅r₁)*H₁₂ .+ exp(-im*k⋅r₂)*H₁₂ 
	#H = H₀ .+ exp(im*k⋅a₁)*H₂₁ .+ exp(-im*k⋅a₁)*H₁₂
	#H = H₀ .+ exp(im*k⋅a₁)*H₂₁ .+ exp(im*k⋅a₂)*H₂₁
	#H = H⊗I(2) + 0.00001*I(2)⊗σ₃
	return H
end


klist = ["Γ", "M", "K", "Γ"]
nk = 1028
println("Getting eigenvalues of graphene between k = ")
show(klist)
println("...")
E, Estates = getBands(klist, nk, a, H)
#display(27.2*E)
println("Plotting...")
plotBands(klist,nk,E)
println("Done! Press ctrl+d to quit")

