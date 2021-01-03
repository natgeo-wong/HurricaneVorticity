### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ b2cd00b6-4c8a-11eb-269a-df1124f1e325
begin
	using DrWatson

	using Pkg
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ b1dc99e6-4c8a-11eb-024e-55952c1da3d3
begin
	@quickactivate "HurricaneVorticity"
	Pkg.instantiate()
	
	using PlutoUI
	
	using FourierFlows
	import GeophysicalFlows.TwoDNavierStokes
	
	using ImageShow, QuartzImageIO
	using PyCall
	using LaTeXStrings
	pplt = pyimport("proplot")
	
md"Loading modules from the HurricaneVorticity project..."
end

# ╔═╡ 0eae2e8a-4c6a-11eb-1084-bfb4de5e213a
md"# 3. Elliptical Vortex Model

The initial elliptical vortex is constructed as follows (in polar coordinates):

$$\zeta(r,\phi,0) = \zeta_0\left\{\begin{array}{ll}
1, & 0\leq r \leq r_i\alpha(\phi) \\
1 - f_\lambda(r
'), & 0\leq r_i\alpha(\phi) \leq r \leq r_0\alpha(\phi) \\
1, & 0\leq r_0\alpha(\phi) \leq r
\end{array}\right.$$"

# ╔═╡ 93da4a92-4c79-11eb-2a63-1967c13ad36e
md"where $\alpha(\phi)$ augments the circle into an ellipse with eccentricity $\epsilon$

$$\alpha(\phi) = \sqrt{\frac{1-\epsilon^2}{1-\epsilon^2\cos^2\phi}},\quad
\epsilon=\sqrt{1-\frac{b^2}{a^2}}$$

And $f_\lambda(r')$ is the smoothing function between the vortex and the transition to the base $\zeta=0$ state of the surrounding atmosphere, given by

$$f_\lambda(r') = e^{-\lambda \exp(1/(r'-1))/r'},\quad
r'=\frac{r-r_i\alpha(\phi)}{r_0\alpha(\phi)-r_i\alpha(\phi)}$$"

# ╔═╡ b034ad44-4cc2-11eb-0db7-c18775ce1a6b
md"We set the dynamic viscosity to be 50 m$^2$ s$^{-1}$."

# ╔═╡ 8844bfae-4cc2-11eb-1e61-0b9e4507b5d8
md"These values are taken and modified from Hendricks & Schubert (2009), \"Transport and mixing in idealized barotropic hurricane-like vortices\", namely the grid coarsened and timestep increased such that it can be run on a normal computer without too much delay.  Should one wish to run these scripts on a cluster, I have provided pure julia scripts in the project under examples that one can use to run."

# ╔═╡ 25efd164-4cc1-11eb-1118-7db974e10d96
md"ngrid  : $(@bind ngrid PlutoUI.Slider((1:8)*128))"

# ╔═╡ 4e31f454-4cc1-11eb-1bc9-43011901ac87
md"dt     : $(@bind dt PlutoUI.Slider((1:3)*5))"

# ╔═╡ 94a34944-4cba-11eb-19bd-4963e38749b3
md"
### A. Setting up the Model

First, we define the problem in with the `TwoDNavierStokes` module provided in the package `GeophysicalFlows.jl`.  We use a 400km x 400km domain, with $(ngrid) x $(ngrid) equally spaced points.  The timestep was taken to be every $(dt) seconds.
"

# ╔═╡ a4947f96-4cbd-11eb-019c-452ffa00e634
prob = TwoDNavierStokes.Problem(
    CPU();stepper="FilteredRK4",
	nx=ngrid,Lx=400e3,ny=ngrid,Ly=400e3,dt=dt,
	ν=50
)

# ╔═╡ 8a20d63c-4c81-11eb-2401-150fb53d2062
md"
### B. Initial vorticity field

We explore different combinations of the initial vorticity field and plot them below.
"

# ╔═╡ 2ed757f0-4c8c-11eb-2017-bfa4824184f1
begin
	# calculate α using x,y cartesian coordinates
	fα(ϵ,x,y) = sqrt((1-ϵ^2)/(1-ϵ^2*(x/sqrt(x^2+y^2))^2))
	
	# convert r to r' based on ri(α) and r0(α)
	r2rp(r,riα,r0α) = (r-riα)/(r0α-riα)
	
	# calculate fλ using r, ri(α) and r0(α) to calculate r'
	function fλ(r,riα,r0α,λ)
		rp = r2rp(r,riα,r0α)
		return exp(-λ/rp*exp(1/(rp-1)))
	end
	
	# calculate vorticity ζ
	function generateζ(ζ0,x,y,ri,r0,λ,ϵ)
		
		nx = length(x); ny = length(y); ζ = ones(nx,ny)
		ygrid = reshape(y,1,:)
		r = sqrt.(x.^2 .+ ygrid.^2)
		α = fα.(ϵ,x,ygrid)
		riα = ri * α
		r0α = ri * α * r0
		
		for iy in 1 : ny, ix in 1 : nx
			
			rii = riα[ix,iy]
			r0i = r0α[ix,iy]
			r_i = r[ix,iy]
			
			if r_i > r0i
				ζ[ix,iy] = 0
			elseif r_i > rii
				ζ[ix,iy] = 1 - fλ(r_i,rii,r0i,λ)
			end
			
		end
		
		return ζ * ζ0
		
	end
	
md"We start by defining the functions for $\alpha$, $r'$, $f_λ$ and $\zeta$"
end

# ╔═╡ a9bb13be-4cb7-11eb-0af3-2ddd92bb4c85
md"We use sliders to vary the following parameters to see how they affect the initial vorticity distribution."

# ╔═╡ 21d4a34c-4c8b-11eb-0cf0-3932e702ac8c
md"ζ0 : $(@bind ζ0 PlutoUI.Slider((2:0.1:8)/1e3))"

# ╔═╡ 3152aae8-4cb4-11eb-05d9-23bf63c2b0ca
md"ri : $(@bind ri PlutoUI.Slider((20:40)*1e3))"

# ╔═╡ 744137da-4c8b-11eb-3ed6-1508719cf982
md"r0 : $(@bind r0 PlutoUI.Slider(1.1:0.1:2))"

# ╔═╡ 7c1636e0-4c8b-11eb-2683-817f32dc96d2
md"λ  : $(@bind λ PlutoUI.Slider(0.2:0.2:2))"

# ╔═╡ e5f09984-4c8b-11eb-302f-9b310bfc73be
md"ϵ  : $(@bind ϵ PlutoUI.Slider(0:0.1:0.9))"

# ╔═╡ 63a178c4-4cc0-11eb-3780-111fdfafd433
md"Next, we initialize the vorticity field $\zeta$ in our problem."

# ╔═╡ 08d6c2f0-4cc0-11eb-2316-ffc56e94ffa4
begin
	sol,clock,vars,grid = prob.sol,prob.clock,prob.vars,prob.grid
	x,y = grid.x,grid.y; xabs = maximum(abs.(x))/1e3; yabs = maximum(abs.(y))/1e3
	ζ = generateζ(ζ0,x,y,ri,r0,λ,ϵ)
	TwoDNavierStokes.set_ζ!(prob,ζ)
end

# ╔═╡ 47a317fa-4cc4-11eb-0ed3-3357245f9257
md"We plot the initial vorticity field $\zeta$ below for our reference and exploration as to how this field changes in the parameter space."

# ╔═╡ d3968802-4cb3-11eb-25e8-d950cc766df6
begin
	pplt.close(); fi,axsi = pplt.subplots(axwidth=2)
    ci = axsi[1].contourf(
		x/1e3,y/1e3,ζ'*1e3,
		cmap="Reds",levels=vcat(0,0.05,0.1,0.2,0.5,1:5),
		extend="max"
	)
    axsi[1].format(suptitle=L"Initial $\zeta$")
    axsi[1].format(xlim=(-xabs/2,xabs/2),ylim=(-yabs/2,yabs/2),xlabel="x",ylabel="y")

    fi.colorbar(ci,loc="r")
	
	if !isdir(plotsdir()); mkpath(plotsdir()) end
    fi.savefig(plotsdir("ellipse-initial.png"),transparent=false,dpi=200)
	load(plotsdir("ellipse-initial.png"))
end

# ╔═╡ 697aef4c-4cc4-11eb-249f-ab8134694953
md"Text here about how $\zeta$ varies in the parameter space ..."

# ╔═╡ dc0e2564-4cc0-11eb-024c-e5245ece9763
md"nstep : $(@bind nstep PlutoUI.Slider((1:10)*360))"

# ╔═╡ 2852a550-4cc4-11eb-12ea-e37bf2da5c69
md"### C. Running the Model

We set the model to run over $(nstep) timesteps, which corresponds to $(nstep*dt/3600) hours."

# ╔═╡ 77808efa-4d6c-11eb-232e-cf9ace21cb11
begin
	stepforward!(prob,nstep)
    TwoDNavierStokes.updatevars!(prob)
	
	pplt.close(); ff,axsf = pplt.subplots(axwidth=2)
    cf = axsf[1].contourf(
		x/1e3,y/1e3,(vars.ζ)'*1e3,
		cmap="Reds",levels=vcat(0,0.05,0.1,0.2,0.5,1:5),
		extend="max"
	)
    axsf[1].format(suptitle=L"Final $\zeta$",ultitle="t=$(nstep*dt/3600) hours")
    axsf[1].format(xlim=(-xabs/2,xabs/2),ylim=(-yabs/2,yabs/2),xlabel="x",ylabel="y")

    ff.colorbar(cf,loc="r")
    ff.savefig(plotsdir("ellipse-final.png"),transparent=false,dpi=200)
	load(plotsdir("ellipse-final.png"))
end

# ╔═╡ Cell order:
# ╟─0eae2e8a-4c6a-11eb-1084-bfb4de5e213a
# ╟─93da4a92-4c79-11eb-2a63-1967c13ad36e
# ╟─b2cd00b6-4c8a-11eb-269a-df1124f1e325
# ╠═b1dc99e6-4c8a-11eb-024e-55952c1da3d3
# ╟─94a34944-4cba-11eb-19bd-4963e38749b3
# ╟─b034ad44-4cc2-11eb-0db7-c18775ce1a6b
# ╟─8844bfae-4cc2-11eb-1e61-0b9e4507b5d8
# ╟─25efd164-4cc1-11eb-1118-7db974e10d96
# ╟─4e31f454-4cc1-11eb-1bc9-43011901ac87
# ╠═a4947f96-4cbd-11eb-019c-452ffa00e634
# ╟─8a20d63c-4c81-11eb-2401-150fb53d2062
# ╠═2ed757f0-4c8c-11eb-2017-bfa4824184f1
# ╟─a9bb13be-4cb7-11eb-0af3-2ddd92bb4c85
# ╟─21d4a34c-4c8b-11eb-0cf0-3932e702ac8c
# ╟─3152aae8-4cb4-11eb-05d9-23bf63c2b0ca
# ╟─744137da-4c8b-11eb-3ed6-1508719cf982
# ╟─7c1636e0-4c8b-11eb-2683-817f32dc96d2
# ╟─e5f09984-4c8b-11eb-302f-9b310bfc73be
# ╟─63a178c4-4cc0-11eb-3780-111fdfafd433
# ╠═08d6c2f0-4cc0-11eb-2316-ffc56e94ffa4
# ╟─47a317fa-4cc4-11eb-0ed3-3357245f9257
# ╠═d3968802-4cb3-11eb-25e8-d950cc766df6
# ╟─697aef4c-4cc4-11eb-249f-ab8134694953
# ╟─2852a550-4cc4-11eb-12ea-e37bf2da5c69
# ╟─dc0e2564-4cc0-11eb-024c-e5245ece9763
# ╟─77808efa-4d6c-11eb-232e-cf9ace21cb11
