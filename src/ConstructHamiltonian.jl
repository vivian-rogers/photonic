

module ConstructHamiltonian
#push!(LOAD_PATH, "./src/")
using LinearAlgebra
using SparseArrays
using Operators
using nearestNeighbors
using GMaterialParameters
using UsefulFunctions
using Constants

export Hgen


function Hgen(λ)
	a₁ = λ*a*[1;0]; a₂ = λ*a*[cosd(60);sind(60)]; 
	R = Rvals(λ)
	N = 2*λ^2
	
	pureGrapheneNNs = nn1(λ,R,t) #returns 1st nn hopping structs
	#one could modify NN radii, hopping, etc here
	
	#thinking about the NS profile
	Δh0 = 2.0*nm; Δxy = 1*nm;
	NSprofile = NSprofileGen(λ)	
	
	H_hop, edgeNNs = nnHoppingMat(pureGrapheneNNs,λ)
	#show(edgeNNs)
	#pseudoB = 0.1*eV*map(r->sin(k⋅r),R)
	#H_pseudoB = 0.1*eV*I(λ^2)⊗σ₃
	#Hstatic = sparse(H_hop.+H_pseudoB) # add on Hpot, Hcharge, H+U, etc here
	Hstatic = sparse(H_hop) # add on Hpot, Hcharge, H+U, etc here
	function H(k)
		N=2*λ^2
		#H_edge = zeros(ComplexF64, N,N)
		H_edge = spzeros(ComplexF64, N,N)
		for NN in edgeNNs
			#println("\nedge NN = \n")
			H_edge[NN.a,NN.b] += NN.t*exp(im*k⋅(NN.N[1]*a₁+NN.N[2]*a₂)) 
			#H_edge[NN.b,NN.a] += NN.t*exp(im*k⋅NN.r) 
			#H_edge[NN.a,NN.b] += NN.t*exp(-im*k⋅NN.r)
			#H_edge[NN.b,NN.a] += NN.t*exp(im*k⋅NN.r)
		end
		return dropzeros(H_edge.+Hstatic)
	end
	return H
end

end
#function hoppingModification(
