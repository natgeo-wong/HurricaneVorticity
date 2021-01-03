using DrWatson
@quickactivate "HurricaneVorticity"
using Crayons.Box
using FourierFlows
import GeophysicalFlows.TwoDNavierStokes

using PyCall
using LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("ellipse.jl"))

### Modify the variables here

ngrid = 512; dt = 5; nstep = 86400; pstep = 72
ζ0 = 5e-3
ri = 20e3
r0 = 3
λ  = 0.5
ϵ  = 0.7

### Model Runs here

prob = TwoDNavierStokes.Problem(
    CPU();stepper="FilteredRK4",
	nx=ngrid,Lx=400e3,ny=ngrid,Ly=400e3,dt=dt,
	ν=50
)

sol,clock,vars,grid = prob.sol,prob.clock,prob.vars,prob.grid
x,y = grid.x,grid.y; xabs = maximum(abs.(x))/1e3; yabs = maximum(abs.(y))/1e3
ζ = generateζ(ζ0,x,y,ri,r0,λ,ϵ)
TwoDNavierStokes.set_ζ!(prob,ζ)

startwalltime = time()
lvl = vcat(0,0.05,0.1,0.2,0.5,1:5)
for j = 0 : round(Int, nstep / pstep)

	if iszero(mod(j*pstep*dt,3600))
		cfl = clock.dt * maximum([maximum(vars.u) / grid.dx, maximum(vars.v) / grid.dy])
	    @info "$(now()) - Running Elliptical Vortex Simulation - $(BOLD("Model Time:")) $(@sprintf("%05.1f",clock.t/3600)) | $(BOLD("CFL:")) $(@sprintf("%4.2f",cfl))"
	end

    pplt.close(); f,axs = pplt.subplots(axwidth=2)
    c = axs[1].contourf(x/1e3,y/1e3,(vars.ζ)'*1e3,cmap="Reds",levels=lvl)
    axs[1].format(ltitle="t = $(@sprintf("%05.1f",clock.t/3600))",rtitle=L"$\zeta$")
    axs[1].format(xlim=(-xabs/2,xabs/2),ylim=(-yabs/2,yabs/2),xlabel="x",ylabel="y")

    f.colorbar(c,loc="r")
    f.savefig(plotsdir("ellipse/$j.png"),transparent=false,dpi=200)

    stepforward!(prob,pstep)
    TwoDNavierStokes.updatevars!(prob)

end
