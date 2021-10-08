push!(LOAD_PATH, "./src/")


using Plots
using PlotStuff
using Constants
using LinearAlgebra
using UsefulFunctions
using Operators
using GMaterialParameters
using ConstructHamiltonian
using nearestNeighbors
using Nanospheres

using Bands


function runBands(nk, H, proj::Bool=false,Q=I(λ))
	klist = ["γ", "μ", "κ", "γ","κ'"]
	kdict = kdictGen(λ)
	println("Getting eigenvalues of SLgraphene (λ=$λ*a₀) between k = ")
	show(klist)
	println("...")
	E, Estates = getBands(klist, kdict, nk, a, H)
	println("Plotting...")
	if(proj)
		projStates = expectedValue(Q,Estates)
		plotBands(klist,nk,E, projStates)
	else
		plotBands(klist,nk,E,2)
	end
end


#multiplicity of a_SL/a_graphene, stddev of NS, Δheight of NS profile
λ = 6; σ = 1*nm; Δh0 = 2*nm;
N = 2*λ^2 #total number of atoms in system
#projKet = I(Int(N/2))⊗[0;1] #project onto B site

#x, y, z = NSprofileGen(λ,Δh0,σ)
#surface(x,y,z)


println("Generating static hamiltonian...")
H = Hgen(λ)
#Q = closestPeak(Rvals(λ),peaksNS(λ),λ)
runBands(2^7,H,false)
#runBands(2^7,H,true,Q)

#=kdict = kdictGen(λ)

klist = ["γ", "μ", "κ", "γ","κ'"]
nk = 2^8

println("Getting eigenvalues of SLgraphene (λ=$λ*a₀) between k = ")
show(klist)
println("...")
E, Estates = getBands(klist, kdict, nk, a, H)
#display(27.2*E)
println("Plotting...")
if(true)
	projStates = expectedValue(closestPeak(Rvals(λ),peaksNS(λ),λ),Estates)
	plotBands(klist,nk,E, projStates)
else
	plotBands(klist,nk,E,2)
end
=#
println("Done! Press ctrl+d to quit")


