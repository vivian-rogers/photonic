
module Operators

using LinearAlgebra
using Constants
using SparseArrays

⊗(A,B) = kron(A,B)
×(u,v) = cross(u,v)

σ₀ = [
      1 0
      0 1
     ]

σ₁ = [
      0  1 
      1  0
     ] 
σ₂ = [
      0 -im 
      im  0
     ] 
σ₃ = [
      1  0
      0 -1
     ]

τ₀ = [
      1 0
      0 1
     ]

τ₁ = [
      0  1 
      1  0
     ] 
τ₂ = [
      0 -im 
      im  0
     ] 
τ₃ = [
      1  0
      0 -1
     ]


S₁ = (1/2)*ħ*σ₁; S₂= (1/2)*ħ*σ₂; S₃ = (1/2)*ħ*σ₃; 
S = cat(S₁,S₂,S₃, dims = 3);
#σ = Array[σ₁;σ₂;σ₃];
x₊ = √(1/2)*[1;1]; x₋ = √(1/2)*[1;-1]; 
y₊ = √(1/2)*[1;im]; y₋ = √(1/2)*[1;-im]; 
z₊ = [1;0]; z₋ = [0;1];



function closestPeak(Rvals,RvalsNS,λ)
	N = 2*λ^2
	peakDist = zeros(N)
	for i in 1:N
		r1 = Rvals[i,:]
		#show(RvalsNS)
		Δx = min.(map(r2->norm(r1-r2),RvalsNS))
		peakDist[i] = Δx[1]
	end
	return sparse(Diagonal(peakDist))
end

Hₙ = Diagonal([-Ry*q^2/n^2 for n = 1:6])


function L₃(l) 
	return Diagonal([ħ*m for m = -l:l])
end


for n in names(@__MODULE__; all=true)
               if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
                   @eval export $n
               end
end

end
