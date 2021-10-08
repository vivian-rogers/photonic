

#push!(LOAD_PATH, "./src/")
module nearestNeighbors

using PyPlot
using Constants
using LinearAlgebra
using UsefulFunctions
using Operators
using GMaterialParameters
using Bands

function Rvals(λ)
	N = (λ^2)
	#println("N = $N\n")
	#δ = (a/sqrt(3))*[cosd(30);sind(30)];
	Rvals = zeros(N,2)
	n = 1
	a₁ = a*[1;0]; a₂ = a*[0;1]; 
	for i = 1:λ
		for j = 1:λ
			x1 = (j-1)*a₁+(i-1)*a₂
			#x2 = x1 + δ
			Rvals[n,:] = x1
			#Rvals[n+1,:] = x2
			n += 2
		end
	end
	return Rvals
end

function Rlist(λ)
	#N = 2*λ^2
	RvalList = Rvals(λ)
	#Rarr = Array[]
	#for i = 1:N
	#	Rarr = Rs[i,:]
	#end
	Rs = [RvalList[i,:] for i in size(RvalList)[1]]
	return Rs 
end

function Rofg(g,λ)
	x = g[1]; y = g[2];
	δ = (a/sqrt(3))*[cosd(30);sind(30)];
	a₁ = λ*a*[1;0]; a₂ = λ*a*[cosd(60);sind(60)]; 
	R = floor((x-1)/2)*a₁+(y-1)*a₂
	if(mod(x,2) == 0)
		R += δ
	end
	return R
end


function g2i(x,y,λ)
	return x + (y-1)*(2*λ)
	# i = x + (y-1)*2λ; (i - x)/2λ + 1= y 
	# i - (y-1)*2λ = x
end


function i2g(i,λ)
	y = ceil(i/(2*λ))
	x = i - (y-1)*(2λ)
	return x,y
end



mutable struct Hopping
	a::Int #site 1
	b::Int #site 2
	r # radius from a to b
	gB # grid-value of point B
	t  # hopping parameter affiliated with <b|a>
	edge::Bool
	N # vector describing the [n₁;n₂]⋅[a₁;a₂] unit cell
end

function nnHoppingMat(NNs,λ)
	N = 2*λ^2	
	H = zeros(N,N)
	edgeNNs = Any[]
	for NN in NNs
		if(NN.edge == true)
			push!(edgeNNs,deepcopy(NN))
		else
			H[NN.a,NN.b] = NN.t
		end
	end
	return H, edgeNNs
end

function nn1(λ,Rvals,t₀)
	N = 2*λ^2	
	NNs = Any[]
	for y = 1:λ
		for x = 1:2:(2*λ) #site A
			i1 = g2i(x,y,λ);
			r1 = [x;y]
			#do NNs
		
			r2 = [x+1,y-1] #down one
			i2 = g2i(r2[1],r2[2],λ);
			NN = Hopping(i1,i2,radius(r1,r2,λ),r2,t₀,false,[0;0])
			push!(NNs,NN)

			r2 = [x-1,y] #left one
			i2 = g2i(r2[1],r2[2],λ);
			NN = Hopping(i1,i2,radius(r1,r2,λ),r2,t₀,false,[0;0])
			push!(NNs,NN)
			
			r2 = [x+1,y] #right one
			i2 = g2i(r2[1],r2[2],λ);
			NN = Hopping(i1,i2,radius(r1,r2,λ),r2,t₀,false,[0;0])
			push!(NNs,NN)
		end
		for x = 2:2:(2*λ) #site B
			i1 = g2i(x,y,λ);
			r1 = [x;y]
			#do NNs
			r2 = [x-1,y+1] #up one
			i2 = g2i(r2[1],r2[2],λ);
			NN = Hopping(i1,i2,radius(r1,r2,λ),r2,t₀,false,[0;0])
			push!(NNs,NN)
			
			r2 = [x-1,y] #left one
			i2 = g2i(r2[1],r2[2],λ);
			NN = Hopping(i1,i2,radius(r1,r2,λ),r2,t₀,false,[0;0])
			push!(NNs,NN)
			
			r2 = [x+1,y] #down one
			i2 = g2i(r2[1],r2[2],λ);
			NN = Hopping(i1,i2,radius(r1,r2,λ),r2,t₀,false,[0;0])
			push!(NNs,NN)
		end
	end
	for NN in NNs # let us fix the edge designation
		x = NN.gB[1]; y = NN.gB[2] # gives grid-value of the adjoint
		#show(NN)
		if(x<1 || x>(2*λ) || y<1 || y>λ)
			
			# generate basis with wraparound
			xnew = mod(x-1,2*λ)+1
			ynew = mod(y-1,λ)+1
			NN.b =g2i(xnew,ynew,λ)
			
			# fill in N vector
			
			if(y<1)
				ny = -1
			elseif(y>λ)
				ny = 1
			else
				ny = 0
			end
			
			if(x<1)
				nx = -1
			elseif(x>2*λ)
				nx = 1
			else
				nx = 0
			end
			NN.N = [nx;ny]
			NN.edge = true
		end
		#println("\n")
	end
	# do nearest neighbor unstrained hopping
	return NNs
end

function radius(g1,g2,λ)
	r1 = Rofg(g1,λ)
	r2 = Rofg(g2,λ) 
	return (r2 - r1)
end

for n in names(@__MODULE__; all=true)
               if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
                   @eval export $n
               end
end

end

