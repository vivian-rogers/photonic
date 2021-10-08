module MaterialParameters

using Constants

#a = 1
# for ScN


a = 2.46*Å
t = 2.8*eV
ε = 0*eV

#real-space lattice vectors
a₁ = a*[1;0]; a₂ = a*[sind(30);cosd(30)]; 
A = hcat(a₁, a₂);
B = transpose((2*π)*inv(A));


#reciprocal lattice
b₁ = B[:,1]; b₂ = B[:,2];

kdict = Dict(
	    "Γ" => B*[0;     0],
	    "K" => B*[1/3; 2/3],
	    "K'" => B*[-1/3; 2/3],
	    "M" => B*[1/2; 0],
	    )	    

function nn1_vects(a)
	vects = Any[]
	for i = 1:3; 
		for j = [-1,1]
			# r₁ vector is defined in src/HolmiumParameters, may warrant reformulation in terms of lattice vectors	
			r = zeros(3);
			r[i] = j*(a/2);
			push!(vects,r)
		end
	end
	return vects
end


function nn2_vects(a)
	vects = Any[]
	for i = 1:3; 
		for j = [-1,1]
			for k = [-1,1]
				r = ones(3);
				#just consider only the two relevant dimensions
				dims = [1,2,3]
				dims = filter!(dim->dim!=i,dims)
				r[i] = 0
				r[dims[1]] = j*a
				r[dims[2]] = k*a
				push!(vects,r)
			end
		end
	end
	return vects
end








#exports all units
for n in names(@__MODULE__; all=true)
               if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
                   @eval export $n
               end
end

end
