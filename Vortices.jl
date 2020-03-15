
using Plots
include("Vortice.jl")
using LaTeXStrings
using BoundaryValueDiffEq





W1=Vortice.VORTEX(0.5,1.0,1,1,4);
W1[6]

plot(W1[1],W1[2],label="f(r)",xlabel=L"R",ls=:solid,lw=[3 1],legend=:bottomright,size=(600,300),title=L"tesis,\lambda=0.3,n=1",color=:blue)
plot!(W1[1],W1[3],label="a(r)",ls=:dash,lw=[3 1],dpi=300,color=:blue,size=(600,300));


plot(W1[1],W1[4],label="B(r)",ls=:solid,lw=[3 1],dpi=300,color=:blue)
plot!(W1[1],W1[5],label="E(r)",ls=:dashdotdot,lw=[3 1],dpi=300,color=:blue);

W2=Vortice.VORTEX(0.5,14.0,2,1,4);

plot(W2[1],W2[2],label="f(r)",xlabel=L"R",ls=:solid,lw=[3 1],legend=:bottomright,size=(600,300),title=L"tesis,\lambda=0.3,n=2",color=:red)
plot!(W2[1],W2[3],label="a(r)",ls=:dash,lw=[3 1],dpi=300,color=:red);


plot(W2[1],W2[4],label="B(r)",ls=:dash,lw=[3 1],dpi=300,color=:red)
plot!(W2[1],W2[5],label="E(r)",ls=:dash,lw=[3 1],dpi=300,color=:red);


W3=Vortice.VORTEX(0.5,14.0,3,1,4);



plot(W3[1],W3[2],label="f(r)",xlabel=L"R",ls=:solid,lw=[3 1],legend=:bottomright,size=(600,300),title=L"tesis,\lambda=0.3,n=3",color=:green)
plot!(W3[1],W3[3],label="a(r)",ls=:dash,lw=[3 1],dpi=300,color=:green);


plot(W3[1],W3[4],label="B(r)",ls=:solid,lw=[3 1],dpi=300,color=:green)
plot!(W3[1],W3[5],label="E(r)",ls=:dashdotdot,lw=[3 1],dpi=300,color=:green);

W4=Vortice.VORTEX(0.5,14.0,4,1,4);



plot(W4[1],W4[2],label="f(r)",xlabel=L"R",ls=:solid,lw=[3 1],legend=:bottomright,size=(600,300),title=L"tesis,\lambda=0.3,n=4",color=:orange)
plot!(W4[1],W4[3],label="a(r)",ls=:dash,lw=[3 1],dpi=300,color=:orange);


plot(W4[1],W4[4],label="B(r)",ls=:dash,lw=[3 1],dpi=300,color=:orange)
plot!(W4[1],W4[5],label="E(r)",ls=:dash,lw=[3 1],dpi=300,color=:orange);


plot(W1[1],W1[2],label=L"f,n=1",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:blue)
plot!(W1[1],W1[3],label=L"a,n=1",ls=:dash,lw=[1 1],dpi=300,color=:blue)
#plot!(st.t,sol3tn1,label=L"B,n=1",ls=:dashdotdot,lw=[3 1],dpi=300,color=:blue)

plot!(W2[1],W2[2],label=L"f,n=2",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:red)
plot!(W2[1],W2[3],label=L"a,n=2",ls=:dash,lw=[1 1],dpi=300,color=:red)
#plot!(st.t,sol3tn2,label=L"B,n=2",ls=:dashdotdot,lw=[3 1],dpi=300,color=:red)

plot!(W3[1],W3[2],label=L"f,n=3",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:green)
plot!(W3[1],W3[3],label=L"a,n=3",ls=:dash,lw=[1 1],dpi=300,color=:green)
#plot!(st.t,sol3tn3,label=L"B,n=3",ls=:dashdotdot,lw=[3 1],dpi=300,color=:green)

plot!(W4[1],W4[2],label=L"f,n=4",xlabel=L"u",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=14.0,\delta=0.5,\kappa=1.44",color=:orange)
perfiles1=plot!(W4[1],W4[3],label=L"a,n=4",ls=:dash,lw=[1 1],dpi=300,color=:orange,legendfontsize=7)
#comparacion2=plot!(st.t,sol3tn4,label=L"B,n=4",ls=:dashdotdot,lw=[3 1],dpi=300,color=:orange)

savefig("perfiles1.png")


plot(W1[1],W1[4],label=L"B,n=1",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:blue)
#plot!(W1[1],W1[5],label=L"E,n=1",ls=:dash,lw=[3 1],dpi=300,color=:blue)
#plot!(st.t,sol3tn1,label=L"B,n=1",ls=:dashdotdot,lw=[3 1],dpi=300,color=:blue)

plot!(W2[1],W2[4],label=L"B,n=2",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:red)
#plot!(W2[1],W2[5],label=L"E,n=2",ls=:dash,lw=[3 1],dpi=300,color=:red)
#plot!(st.t,sol3tn2,label=L"B,n=2",ls=:dashdotdot,lw=[3 1],dpi=300,color=:red)

plot!(W3[1],W3[4],label=L"B,n=3",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:green)
#plot!(W3[1],W3[5],label=L"E,n=3",ls=:dash,lw=[3 1],dpi=300,color=:green)
#plot!(st.t,sol3tn3,label=L"B,n=3",ls=:dashdotdot,lw=[3 1],dpi=300,color=:green)

campos1=plot!(W4[1],W4[4],label=L"B,n=4",xlabel=L"u",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=14.0,\delta=0.5,\kappa=1.44",color=:orange,dpi=300)
#plot!(W4[1],W4[5],label=L"E,n=4",ls=:dash,lw=[3 1],dpi=300,color=:orange)
#comparacion2=plot!(st.t,sol3tn4,label=L"B,n=4",ls=:dashdotdot,lw=[3 1],dpi=300,color=:orange)

savefig("campos1.png")


#plot(W1[1],W1[4],label=L"B,n=1",xlabel=L"R",ls=:solid,lw=[3 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:blue)
plot(W1[1],W1[5],label=L"\mathcal{E},n=1",ls=:dash,lw=[1 1],dpi=300,color=:blue)
#plot!(st.t,sol3tn1,label=L"B,n=1",ls=:dashdotdot,lw=[3 1],dpi=300,color=:blue)

#plot!(W2[1],W2[4],label=L"B,n=2",xlabel=L"R",ls=:solid,lw=[3 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:red)
#plot!(W2[1],W2[5],label=L"\mathcal{E},n=2",ls=:dash,lw=[1 1],dpi=300,color=:red)
#plot!(st.t,sol3tn2,label=L"B,n=2",ls=:dashdotdot,lw=[3 1],dpi=300,color=:red)

#plot!(W3[1],W3[4],label=L"B,n=3",xlabel=L"R",ls=:solid,lw=[3 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:green)
#plot!(W3[1],W3[5],label=L"\mathcal{E},n=3",ls=:dash,lw=[1 1],dpi=300,color=:green)
#plot!(st.t,sol3tn3,label=L"B,n=3",ls=:dashdotdot,lw=[3 1],dpi=300,color=:green)

#plot!(W4[1],W4[4],label=L"B,n=4",xlabel=L"R",ls=:solid,lw=[3 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:orange,dpi=300)
#densidad1=plot!(W4[1],W4[5],label=L"\mathcal{E},n=4",xlabel=L"u",ls=:dash,lw=[1 1],dpi=300,color=:orange,title=L"\lambda=14.0,\delta=0.5,\kappa=1.44",size=(600,300))
#comparacion2=plot!(st.t,sol3tn4,label=L"B,n=4",ls=:dashdotdot,lw=[3 1],dpi=300,color=:orange)

savefig("densidad1.png")

w1=Vortice.VORTEX(0.5,4.0,1,1,4);
w2=Vortice.VORTEX(0.5,4.0,2,1,4);
w3=Vortice.VORTEX(0.5,4.0,3,1,4);
w4=Vortice.VORTEX(0.5,4.0,4,1,4);
w1[6]

plot(w1[1],w1[2],label=L"f,n=1",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:blue)
plot!(w1[1],w1[3],label=L"a,n=1",ls=:dash,lw=[1 1],dpi=300,color=:blue)
#plot!(st.t,sol3tn1,label=L"B,n=1",ls=:dashdotdot,lw=[3 1],dpi=300,color=:blue)

plot!(w2[1],w2[2],label=L"f,n=2",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:red)
plot!(w2[1],w2[3],label=L"a,n=2",ls=:dash,lw=[1 1],dpi=300,color=:red)
#plot!(st.t,sol3tn2,label=L"B,n=2",ls=:dashdotdot,lw=[3 1],dpi=300,color=:red)

plot!(w3[1],w3[2],label=L"f,n=3",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:green)
plot!(w3[1],w3[3],label=L"a,n=3",ls=:dash,lw=[1 1],dpi=300,color=:green)
#plot!(st.t,sol3tn3,label=L"B,n=3",ls=:dashdotdot,lw=[3 1],dpi=300,color=:green)

plot!(w4[1],w4[2],label=L"f,n=4",xlabel=L"u",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=4.0,\delta=0.5,\kappa=0.77",color=:orange)
perfiles2=plot!(w4[1],w4[3],label=L"a,n=4",ls=:dash,lw=[1 1],dpi=300,color=:orange,legendfontsize=7,size=(600,300))
#comparacion2=plot!(st.t,sol3tn4,label=L"B,n=4",ls=:dashdotdot,lw=[3 1],dpi=300,color=:orange)


savefig("perfiles2.png")

plot(w1[1],w1[4],label=L"B,n=1",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:blue)
#plot!(W1[1],W1[5],label=L"E,n=1",ls=:dash,lw=[3 1],dpi=300,color=:blue)
#plot!(st.t,sol3tn1,label=L"B,n=1",ls=:dashdotdot,lw=[3 1],dpi=300,color=:blue)

plot!(w2[1],w2[4],label=L"B,n=2",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:red)
#plot!(W2[1],W2[5],label=L"E,n=2",ls=:dash,lw=[3 1],dpi=300,color=:red)
#plot!(st.t,sol3tn2,label=L"B,n=2",ls=:dashdotdot,lw=[3 1],dpi=300,color=:red)

plot!(w3[1],w3[4],label=L"B,n=3",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:green)
#plot!(W3[1],W3[5],label=L"E,n=3",ls=:dash,lw=[3 1],dpi=300,color=:green)
#plot!(st.t,sol3tn3,label=L"B,n=3",ls=:dashdotdot,lw=[3 1],dpi=300,color=:green)

campos2=plot!(w4[1],w4[4],label=L"B,n=4",xlabel=L"u",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=4.0,\delta=0.5,\kappa=0.77",color=:orange,dpi=300)
#plot!(W4[1],W4[5],label=L"E,n=4",ls=:dash,lw=[3 1],dpi=300,color=:orange)
#comparacion2=plot!(st.t,sol3tn4,label=L"B,n=4",ls=:dashdotdot,lw=[3 1],dpi=300,color=:orange)

savefig("campos2.png")


#plot(W1[1],W1[4],label=L"B,n=1",xlabel=L"R",ls=:solid,lw=[3 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:blue)
plot(w1[1],w1[5],label=L"\mathcal{E},n=1",ls=:dash,lw=[1 1],dpi=300,color=:blue)
#plot!(st.t,sol3tn1,label=L"B,n=1",ls=:dashdotdot,lw=[3 1],dpi=300,color=:blue)

#plot!(W2[1],W2[4],label=L"B,n=2",xlabel=L"R",ls=:solid,lw=[3 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:red)
plot!(w2[1],w2[5],label=L"\mathcal{E},n=2",ls=:dash,lw=[1 1],dpi=300,color=:red)
#plot!(st.t,sol3tn2,label=L"B,n=2",ls=:dashdotdot,lw=[3 1],dpi=300,color=:red)

#plot!(W3[1],W3[4],label=L"B,n=3",xlabel=L"R",ls=:solid,lw=[3 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:green)
plot!(w3[1],w3[5],label=L"\mathcal{E},n=3",ls=:dash,lw=[1 1],dpi=300,color=:green)
#plot!(st.t,sol3tn3,label=L"B,n=3",ls=:dashdotdot,lw=[3 1],dpi=300,color=:green)

#plot!(W4[1],W4[4],label=L"B,n=4",xlabel=L"R",ls=:solid,lw=[3 1],legend=:bottomright,size=(600,300),title=L"\lambda=0.3",color=:orange,dpi=300)
densidad2=plot!(w4[1],w4[5],label=L"\mathcal{E},n=4",ls=:dash,lw=[1 1],dpi=300,color=:orange,xlabel=L"u",title=L"\lambda=4.0,\delta=0.5,\kappa=0.77",size=(600,300))
#comparacion2=plot!(st.t,sol3tn4,label=L"B,n=4",ls=:dashdotdot,lw=[3 1],dpi=300,color=:orange)

savefig("densidad2.png")

q1=Vortice.VORTEX(0.5,0.25,1,1,4);
#q2=VORTEX(0.5,3.0,1,1,4);
#q3=VORTEX(0.5,4.0,1,1,4);
#q4=VORTEX(0.5,10.0,1,1,4);

q1[6],q2[6],q3[6],q4[6]



#plot(q1[1],q1[2],label=L"f,\lambda=0.25,\kappa=0.27",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"n=1,\delta=0.5",color=:red)
#plot!(q1[1],q1[3],label=L"a,\lambda=0.25,\kappa=0.27",ls=:dash,lw=[1 1],dpi=300,color=:red)
#plot!(st.t,sol3tn1p,label=L"B,\lambda=0.5",ls=:dashdotdot,lw=[3 1],dpi=300,color=:red)
plot(q2[1],q2[2],label=L"f,\lambda=3.0,\kappa=0.66",xlabel=L"R",ls=:solid,lw=[.5 .5],legend=:bottomright,size=(600,300),title=L"n=1,\delta=0.5",color=:blue)
plot!(q2[1],q2[3],label=L"a,\lambda=3.0,\kappa=0.66",ls=:dash,lw=[.5 .5],dpi=300,color=:blue)
#plot!(st.t,sol3tn1,label=L"B,\lambda=1.0",ls=:dashdotdot,lw=[3 1],dpi=300,color=:blue)
plot!(q3[1],q3[2],label=L"f,\lambda=4.0,\kappa=0.77",xlabel=L"R",ls=:solid,lw=[.5 .5],legend=:bottomright,size=(600,300),title=L"n=1,\delta=0.5",color=:green)
plot!(q3[1],q3[3],label=L"a,\lambda=4.0,\kappa=0.77",ls=:dash,lw=[.5 .5],dpi=300,color=:green)
#plot!(st.t,sol3tn1p1,label=L"B,\lambda=1.5",ls=:dashdotdot,lw=[3 1],dpi=300,color=:green)
plot!(q4[1],q4[2],label=L"f,\lambda=10.0,\kappa=1.21",xlabel=L"u",ls=:solid,lw=[.5 .5],legend=:bottomright,size=(600,300),title=L"n=1,\delta=0.5",color=:orange)
perfiles3=plot!(q4[1],q4[3],label=L"a,\lambda=10.0,\kappa=1.21",ls=:dash,lw=[.5 .5],dpi=300,color=:orange,legendfontsize=5,size=(600,300))


savefig("perfiles3.png")

plot(q1[1],q1[4],label=L"B,\lambda=0.25,\kappa=0.19",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"n=1,\delta=0.5",color=:red)
#plot!(W1[1],W1[3],label=L"a,\lambda=0.5",ls=:dash,lw=[3 1],dpi=300,color=:red)
#plot!(st.t,sol3tn1p,label=L"B,\lambda=0.5",ls=:dashdotdot,lw=[3 1],dpi=300,color=:red)
#plot!(q2[1],q2[4],label=L"B,\lambda=1.0,\kappa=1.30",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"n=1,\delta=0.5",color=:blue)
#plot!(W2[1],W2[3],label=L"a,\lambda=1.0",ls=:dash,lw=[3 1],dpi=300,color=:blue)
#plot!(st.t,sol3tn1,label=L"B,\lambda=1.0",ls=:dashdotdot,lw=[3 1],dpi=300,color=:blue)
plot!(q3[1],q3[4],label=L"B,\lambda=4.0,\kappa=0.77",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"n=1,\delta=0.5",color=:green,dpi=300)
#plot!(W3[1],W3[3],label=L"a,\lambda=1.5",ls=:dash,lw=[3 1],dpi=300,color=:green)
#plot!(st.t,sol3tn1p1,label=L"B,\lambda=1.5",ls=:dashdotdot,lw=[3 1],dpi=300,color=:green)
campos3=plot!(q4[1],q4[4],label=L"B,\lambda=10.0,\kappa=1.21",xlabel=L"u",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"n=1,\delta=0.5",color=:orange,dpi=300)
#plot!(W4[1],W4[3],label=L"a,\lambda=4.0",ls=:dash,lw=[3 1],dpi=300,color=:orange)

savefig("campos3.png")

plot(q2[1],q2[5],label=L"\mathcal{E},\lambda=0.25,\kappa=0.27",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"n=1,\delta=0.5",color=:red)
#plot!(W1[1],W1[3],label=L"a,\lambda=0.5",ls=:dash,lw=[3 1],dpi=300,color=:red)
#plot!(st.t,sol3tn1p,label=L"B,\lambda=0.5",ls=:dashdotdot,lw=[3 1],dpi=300,color=:red)
#plot!(q2[1],q2[5],label=L"\mathcal{E},\lambda=1.0,\kappa=1.30",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"n=1,\delta=0.5",color=:blue)
#plot!(W2[1],W2[3],label=L"a,\lambda=1.0",ls=:dash,lw=[3 1],dpi=300,color=:blue)
#plot!(st.t,sol3tn1,label=L"B,\lambda=1.0",ls=:dashdotdot,lw=[3 1],dpi=300,color=:blue)
plot!(q3[1],q3[5],label=L"\mathcal{E},\lambda=2.0,\kappa=0.77",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"n=1,\delta=0.5",color=:green,dpi=300)
#plot!(W3[1],W3[3],label=L"a,\lambda=1.5",ls=:dash,lw=[3 1],dpi=300,color=:green)
#plot!(st.t,sol3tn1p1,label=L"B,\lambda=1.5",ls=:dashdotdot,lw=[3 1],dpi=300,color=:green)
densidad3=plot!(q4[1],q4[5],label=L"\mathcal{E},\lambda=7.0,\kappa=1.44",xlabel=L"u",ls=:solid,lw=[1 1],legend=:topright,size=(600,300),title=L"n=1,\delta=0.5",color=:orange,legendfontsize=8)
#plot!(W4[1],W4[3],label=L"a,\lambda=4.0",ls=:dash,lw=[3 1],dpi=300,color=:orange)

savefig("densidad3.png")

Q1=Vortice.VORTEX(0.001,4.0,1,1,4);
Q2=Vortice.VORTEX(0.5,4.0,1,1,4);
Q3=Vortice.VORTEX(9/8,4.0,1,1,4);


Q1[6],Q3[6]

gr()
plot(Q1[1],Q1[2],label=L"f,\delta=10^{-3},\kappa=0.75",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=4.0,n=1",color=:blue)
plot!(Q1[1],Q1[3],label=L"a,\delta=10^{-3},\kappa=0.75",ls=:dash,lw=[1 1],dpi=300,color=:blue)
#plot!(st.t,sol3tn1p1,label=L"a,k=0.5",ls=:dashdotdot,lw=[3 1],dpi=300,color=:blue)

#plot!(Q2[1],Q2[2],label=L"f,\delta=0.5,\kappa=1.30",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=1.0,n=1",color=:red)
#plot!(Q2[1],Q2[3],label=L"a,\delta=0.5,\kappa=1.30",ls=:dash,lw=[1 1],dpi=300,color=:red)
#plot!(st.t,sol2tn1p3,label=L"B,k=10^{-3}",ls=:dashdotdot,lw=[3 1],dpi=300,color=:red)

plot!(Q3[1],Q3[2],label=L"f,\delta=9/8,\kappa=0.79",xlabel=L"u",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=4.0,n=1",color=:red)
perfiles4=plot!(Q3[1],Q3[3],label=L"a,\delta=9/8,\kappa=0.79",ls=:dash,lw=[1 1],dpi=300,color=:red)
#comparación4=plot!(st.t,sol3tn1p3,label=L"B,k=10^{-3}",ls=:dashdotdot,lw=[3 1],dpi=300,color=:green)



savefig("perfiles4.png")


plot(Q1[1],Q1[4],label=L"B,\delta=10^{-3},\kappa=1.33",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=1.0,n=1",color=:blue)
#plot!(W1[1],W1[3],label=L"a,k=10^{-3}",ls=:dash,lw=[3 1],dpi=300,color=:blue)
#plot!(st.t,sol3tn1p1,label=L"a,k=0.5",ls=:dashdotdot,lw=[3 1],dpi=300,color=:blue)

#plot!(Q2[1],Q2[4],label=L"B,\delta=0.5,\kappa=1.30",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=1.0,n=1",color=:red)
#plot!(W2[1],W2[3],label=L"a,k=0.5",ls=:dash,lw=[3 1],dpi=300,color=:red)
#plot!(st.t,sol2tn1p3,label=L"B,k=10^{-3}",ls=:dashdotdot,lw=[3 1],dpi=300,color=:red)

campos4=plot!(Q2[1],Q3[4],label=L"B,\delta=9/8,\kappa=1.26",xlabel=L"u",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=2.0,n=1",color=:red,dpi=300)
#plot!(W2[1],W3[3],label=L"a,k=9/8",ls=:dash,lw=[3 1],dpi=300,color=:green)
#comparación4=plot!(st.t,sol3tn1p3,label=L"B,k=10^{-3}",ls=:dashdotdot,lw=[3 1],dpi=300,color=:green)

savefig("campos4.png")

plot(Q1[1],Q1[5],label=L"\mathcal{E},\delta=10^{-3},\kappa=1.33",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=1.0,n=1",color=:blue)
#plot!(W1[1],W1[3],label=L"a,k=10^{-3}",ls=:dash,lw=[3 1],dpi=300,color=:blue)
#plot!(st.t,sol3tn1p1,label=L"a,k=0.5",ls=:dashdotdot,lw=[3 1],dpi=300,color=:blue)

#plot!(Q2[1],Q2[5],label=L"\mathcal{E},\delta=0.5,\kappa=1.30",xlabel=L"R",ls=:solid,lw=[1 1],legend=:bottomright,size=(600,300),title=L"\lambda=1.0,n=1",color=:red)
#plot!(W2[1],W2[3],label=L"a,k=0.5",ls=:dash,lw=[3 1],dpi=300,color=:red)
#plot!(st.t,sol2tn1p3,label=L"B,k=10^{-3}",ls=:dashdotdot,lw=[3 1],dpi=300,color=:red)

densidad4=plot!(Q2[1],Q3[5],label=L"\mathcal{E},\delta=9/8,\kappa=1.26",xlabel=L"R",ls=:solid,lw=[1 1],legend=:topright,size=(600,300),title=L"\lambda=2.0,n=1",color=:red,dpi=300)
#plot!(W2[1],W3[3],label=L"a,k=9/8",ls=:dash,lw=[3 1],dpi=300,color=:green)
#comparación4=plot!(st.t,sol3tn1p3,label=L"B,k=10^{-3}",ls=:dashdotdot,lw=[3 1],dpi=300,color=:green)

savefig("densidad4.png")
