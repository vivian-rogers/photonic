

#push!(LOAD_PATH, "./src/")
module Nanospheres

#using Plots
#using PyPlot
using Constants
using LinearAlgebra
using UsefulFunctions
using Operators
using GMaterialParameters
using Bands
using ForwardDiff
#using Interpolations

export peaksNS, NSprofileGen, strain, testhess

function peaksNS(λ,θ=0) #generates the positions of the nanospheres at each of the corners of the unit cell
	a₁ = λ*a*[1;0]; a₂ = λ*a*[cosd(60);sind(60)]; 
	#needs work to rotate
	#RvalsNS = Array[]
	#RvalsNS = Array{Array{}(1),}(1)
	#RvalsNS = Array{Array{Float64, 2}}(4)
	#RvalsNS = Array{::Array{Float64,1},1}(undef, 4)
	RvalsNS = zeros(4,2)
	for i = 0:1; 
		for j = 0:1
			x = i*a₁ + j*a₂
			RvalsNS[2*j+i+1,:] = x
			#push!(RvalsNS,x)
		end
	end
	return [RvalsNS[i, :] for i in 1:size(RvalsNS, 1)]
end

#function height(Δh0,σ,RvalsNS,x,y)
#	return sum([Δh0*Gaussian(σ,[x,y],RvalsNS[i,:]) for i = 1:4]);

function NSprofileGen(λ,Δh0,σ,θ=0)
	Δedge = 2*nm
	#Δx = 0.05*nm
	nGrid = 10*λ
	RvalsNS = peaksNS(λ)
	MaxR = maximum(RvalsNS[4][:])
	#xVals = [(0-Δedge):Δx:(MaxR[1]+Δedge)] # create vector of Xvals over grid
	#yVals = [(0-Δedge):Δx:(MaxR[2]+Δedge)] # create vector of Xvals over grid

	xVals = range(0 - Δedge, stop=MaxR + Δedge, length=nGrid);
	yVals = range(0 - Δedge, stop=MaxR + Δedge, length=nGrid);

	z(R) = sum([Δh0*Gaussian(σ,R,RvalsNS[i]) for i = 1:4]);
	zVals = [z([x,y]) for y in yVals, x in xVals];
	#surface(xVals,yVals,z)
	return xVals, yVals, zVals, z
end

function Gaussian(σ,R,Rcenter)
	C = 1
	#C = 1/(σ*√(2*π))
	#show(R)
	#show(Rcenter)
	return Float64(C*exp((-1/2)*norm(R - Rcenter)^2/σ^2))
end	

function testhess(z)
	dux_dx = R-> a*ForwardDiff.hessian(z,R)[1,1]
	return dux_dx
end


function strain(λ,Δh0,fᵤ,σ,z)
	#gridpts = CartProd(x,y)
	
	dh_dx = R -> ForwardDiff.gradient(z,R)[1]
	dh_dy = R -> ForwardDiff.gradient(z,R)[2]
	
	# generate in-plane displacement profile
	u = R -> -fᵤ*a*ForwardDiff.gradient(z,R)
	# we don't actually need this though, just the derivatives
	# little bit jank but it gets the job done
	dux_dx = R-> -fᵤ*a*ForwardDiff.hessian(z,R)[1,1]
	duy_dx = R-> -fᵤ*a*ForwardDiff.hessian(z,R)[1,2]
	dux_dy = R-> -fᵤ*a*ForwardDiff.hessian(z,R)[2,1]
	duy_dy = R-> -fᵤ*a*ForwardDiff.hessian(z,R)[2,2]

	#dux1_dx2 = d2u_d2.(Rvals)
	ε_xx = R -> dux_dx(R) + (1/2)*(dh_dx(R))^2
	ε_yy = R -> duy_dy(R) + (1/2)*(dh_dy(R))^2
	ε_xy = R -> (1/2)*(dux_dy(R)+duy_dx(R)) + (1/2)*(dh_dx(R)*dh_dy(R))
	ε_eff = R -> ε_xx(R) + ε_yy(R)
	κ = 1; β = -3.37; Φ₀ = h/q;
	C = κ*β*Φ₀/(a^2*4*π)
	Ax = R -> C*(ε_xx(R)-ε_yy(R))
	Ay = R -> -C*(2*ε_xy(R))
	pseudoB = R -> ForwardDiff.gradient(Ay,R)[1] - ForwardDiff.gradient(Ax,R)[2]

	# in-plane displacement profile
	#u = -a*fᵤ*ForwardDiff.gradient(z,[x,y])
	return pseudoB, dux_dx, ε_eff
end

end
