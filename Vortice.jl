__precompile__()
module Vortice
 using Pkg
 Pkg.add("BoundaryValueDiffEq")
 using BoundaryValueDiffEq

 function VORTEX(delta,lambda,m,u1,u2)
  tf=12.0
  v1=u1
  v2=u2
  λ=lambda/2
  k=delta
  n=m
  e=1.0
  ti=1e-5
  tspan = (ti,tf)
  vm=((3*v1+v2)-sqrt(-8*k+(v2-v1)^2))/4
  vM=((3*v1+v2)+sqrt(-8*k+(v2-v1)^2))/4
  masaA= sqrt(2)*e*v2
  masaH=sqrt((λ)*(k+(v2-v1)^2))
  kappa=masaH/masaA
  function vorticetesis!(dx,x,p,t)
     dx[1] = x[2]
     dx[2] = -x[2]/t+x[1]*(1-x[3])^2*(n^2/t^2)+(2*λ/e^2)*(x[1]-1)*(x[1]-(vM/v2))*(x[1]-(vm/v2))

     dx[3] = x[4]
     dx[4] = x[4]/t-(2.0*x[1]^2)*(1-x[3])
  end


  function bc1!(residual, u, p, t)
     residual[1] = u[1][1] - 0.0
     residual[2] = u[end][1] - (1-(1/sqrt(tf))*exp(-sqrt(2)*(masaA/masaH)*tf))
     residual[3] = u[1][3] - 0.0
     residual[4] = u[end][3] - (1-sqrt(tf)*exp(-sqrt(2)*tf))

  end
  bvp1 = BVProblem(vorticetesis!, bc1!, [0.01,1.0,0.0,0.005], tspan)

  st = solve(bvp1, GeneralMIRK4(), dt=0.1)
  CampoF=[]
  CampoA=[]
  CampoB=[]
  E=[]
  for i in 1:length(st.u)
     push!(CampoF,st.u[i][1])
     push!(CampoA,st.u[i][3])
     push!(CampoB,(n/e)*st.u[i][4]/st.t[i])
     energy=(v2^2)*((n^2/(2*st.t[i]^2))*(st.u[i][4])^2+st.u[i][2]^2+(n^2/st.t[i]^2)*((1-st.u[i][3])^2)*st.u[i][1]^2+(λ/e^2)*((st.u[i][1]-1)^2)*(delta+(st.u[i][1]-(v1/v2))^2))
     push!(E,energy)
  end
    return [st.t,CampoF,CampoA,CampoB,E,kappa]
 end
end
