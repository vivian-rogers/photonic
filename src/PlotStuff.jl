
module PlotStuff
using Plots
using PyPlot
using ColorSchemes
using Constants

export plotBands, plotVec, plotSurf, plotFunct

function plotVec(x,yvecs, title)
	nY = size(yvecs)[1]
	for i = 1:nY
		plot(x,yvecs[i])
	end
	title(title)
	gcf()
end

function plotBands(klist, nk, E, projStates)
	#Plots.pyplot()
	nSymPts = size(klist)[1]
	indices = LinRange(0,nSymPts-1, size(E)[1])
	nE = size(E)[2]
	#display(E[:,1])
	#display(plot!(indices,E[1,:]))
	#display(plot!(indices,E[:,2]))
	#Eplot = transpose(E)
	kSymPts = [i for i =0:(nSymPts-1)]
	for kTick in kSymPts
		plot([kTick,kTick],[-30,30],c="#666666",lw=0.5)
	end
	xticks(kSymPts,klist)
	xlim(0, nSymPts-1)
	maxE = maximum(E)
	minE = minimum(E)
	ylabel("E - Ef (eV)")
	ylim(minE-0.5,maxE+0.5)


	set_cmap("rainbow")
	for iE = 1:nE
		Evals = collect(E[:,iE])
		Projvals = collect(projStates[:,iE])
		#Projvals = collect(projStates[:,iE])
		#scatter(indices,Evals,c=Projvals, s=0.9)
		scatter(indices,Evals,c=Projvals, vmin=0, s=0.9)
		#scatter(indices,Evals,c=Projvals, vmin=0, vmax=1,s=0.9)
		#display(plot!(indices,Evals))
	end
	gcf()
end

function plotPoints(Rvals,z,xlab="",ylab="",zlab="",name="",cmap= :redgreensplit,save=false)

	dx = maximum(x)-minimum(x); dy = maximum(y)-minimum(y)
	C = 500
	width = C
	height = C*dy/dx
	surf = scatter(Rvals[:,1],Rvals[:,2],z, xlabel=xlab, ylabel=ylab, title=name, c = cmap, size=(width,height))
	gui(surf)
end

function plotSurf(x,y,z,xlab="",ylab="",name="",cmap= :redgreensplit,save=false)

	dx = maximum(x)-minimum(x); dy = maximum(y)-minimum(y)
	C = 500
	width = C
	height = C*dy/dx
	surf = heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, c = cmap, size=(width,height))
	#heatmap!
	gui(surf)
end

function plotScatter(R,z,xlab="",ylab="",name="",cmap= :redgreensplit,xyscale=(1/nm),zscale=1,save=false)
	dx = maximum.(R)-minimum.(R); dy = maximum.(R)-minimum.(R)
	C = 500
	width = C
	height = C*dy/dx
	fig, ax = PyPlot.subplots();
	#ax.plot(x,y)
	zplot = ax.scatter(Tuple([xyscale*r[1] for r in R]),Tuple([xyscale*r[2] for r in R]), c=Tuple([zscale*zi for zi in z]))
	PyPlot.xlabel(xlab);
	#fig.colorbar(zplot, ax=ax);
	PyPlot.colorbar(zplot, label=name);
	PyPlot.ylabel(ylab);
	show(surf)
end


function plotFunct(R,f::Function,xlab="",ylab="",name="",cmap= :redgreensplit,xyscale=(1/nm),zscale=1,save=false)
	dx = maximum.(R)-minimum.(R); dy = maximum.(R)-minimum.(R)
	C = 500
	width = C
	height = C*dy/dx
	z = f.(R)
	#show(z)
	#surf = Plots.scatter(R[:][1],R[:][2],z, xlabel=xlab, ylabel=ylab, title=name, c = cmap, size=(width,height))
	
	#surf = Plots.scatter(R[:,1],R[:,2], marker=:c, color=cmap, zcolor=z, size=(width,height))
	#show([r[1] for r in R])
	#surf = Plots.scatter([r[1] for r in R],[r[2] for r in R], marker=:c, color=cmap, zcolor=z, size=(width,height))
	fig, ax = PyPlot.subplots();
	#ax.plot(x,y)
	zplot = ax.scatter(Tuple([xyscale*r[1] for r in R]),Tuple([xyscale*r[2] for r in R]), c=zscale*z, s = 0.8)
	PyPlot.xlabel(xlab);
	#fig.colorbar(zplot, ax=ax);
	PyPlot.colorbar(zplot, label=name);
	PyPlot.ylabel(ylab);
	#PyPlot.title(name);
	#surf.title(name)
	#surf.xlabel(xlab)
	#surf.ylabel(ylab)
	#scatter(Tuple.(points))	
	#show(R[:])
	#surf = Plots.scatter(R[1][2],R[:][2], marker=:c, color=cmap, zcolor=z, size=(width,height))
	#surf = PyPlot.scatter(R[1],R[2], c=z)
	#heatmap!
	show(surf)
end


function plotBands(klist, nk, E)
	#Plots.pyplot()
	nSymPts = size(klist)[1]
	indices = LinRange(0,nSymPts-1, size(E)[1])
	nE = size(E)[2]
	#display(E[:,1])
	#display(plot!(indices,E[1,:]))
	#display(plot!(indices,E[:,2]))
	#Eplot = transpose(E)
	kSymPts = [i for i =0:(nSymPts-1)]
	for kTick in kSymPts
		plot([kTick,kTick],[-30,30],c="#666666",lw=0.5)
	end
	xticks(kSymPts,klist)
	xlim(0, nSymPts-1)
	maxE = maximum(E)
	minE = minimum(E)
	#maxE = Emax
	#minE = -Emax
	ylabel("E - Ef (eV)")
	ylim(minE,maxE)
	for iE = 1:nE
		Evals = collect(E[:,iE])
		#plot(indices,Evals)
		PyPlot.scatter(indices,Evals,s=1)
		#display(plot!(indices,Evals))
	end
	gcf()
end

end
