
using Plots
using LaTeXStrings
include("Vortice.jl")

Pkg.add("")



W1=Vortice.VORTEX(0.5,1.0,1,1,4);
W1[6]

plot(W1[1],W1[2],label="f(r)",xlabel=L"R",ls=:solid,lw=[3 1],legend=:bottomright,size=(600,300),title=L"tesis,\lambda=0.3,n=1",color=:blue)
prueba1=plot!(W1[1],W1[3],label="a(r)",ls=:dash,lw=[3 1],dpi=300,color=:blue,size=(600,300));

savefig("Prueba1.png")
