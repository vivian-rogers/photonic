push!(LOAD_PATH, "./src/")

#using Plots
using nearestNeighbors
using Nanospheres
using Constants
using UsefulFunctions
using PlotStuff

λ = 80; σ = 4*10^-9; Δh0 = 10^-9; fᵤ = 1

x,y,z,fz = NSprofileGen(λ,Δh0,σ)
R = MtoA(Rvals(λ))
#dux_dx = testhess(fz)
#psB, u = strain(λ,Δh0,fᵤ,σ,fz)
#Ax, Ay, pseudoB = strain(λ,Δh0,fᵤ,σ,Rvals,z)
#testh = fz.(R)
#show(R)
#testf = dux_dx.(R)
plotFunct(R,fz,"x position (nm)","y position (nm)","pseudoB",:inferno,(1/nm),(1/nm))

#pseudoB = 
#surf = surface(x,y,z)

#plotSurf((1/nm)*x,(1/nm)*y,(1/nm)*z,"x position (nm)", "y position (nm)","Height profile", :plasma)
#surf = plot(x,y,z,st=:surface)
#gui(surf)
#surface(x,y,z)

