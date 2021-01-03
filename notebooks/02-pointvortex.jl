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

# ╔═╡ 4fe6a84c-4c69-11eb-295f-b7052b606cc3
begin
	using DrWatson
	using Pkg
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 56ce9c32-4c69-11eb-2ebd-79c3f4f19990
begin
	@quickactivate "HurricaneVorticity"
	Pkg.instantiate()
	
	using PlutoUI
	
	using ImageShow, QuartzImageIO
	using PyCall
	using LaTeXStrings
	pplt = pyimport("proplot");
	
md"Loading modules from the HurricaneVorticity project..."
end

# ╔═╡ 0a3bc67c-4c69-11eb-1dbf-37273535c61c
md"# 2. Point Vortex Models

This model of a hurricane idealizes its structure such that its vorticity is concentrated into a point at its center with circulation strength $\Gamma$."

# ╔═╡ 72969256-4c69-11eb-0c43-fba850561b78
 md"Γ : $(@bind Γ PlutoUI.Slider(5:15))"

# ╔═╡ 7a1407e0-4c69-11eb-3d67-b56be81e8f83
md"resolution : $(@bind res PlutoUI.Slider(0.01:0.01:0.1))"

# ╔═╡ a26602e8-4c69-11eb-3bb9-21a107f594f0
begin
	x = collect(-2:res:2)
	y = reshape(collect(-2:res:2),1,:)
	v = Γ ./ (2*pi*sqrt.(x.^2 .+ y.^2))
	
	xw = collect(-1:0.2:1)
	yw = reshape(collect(-1:0.2:1),1,:)
	hw = Γ ./ (2*pi*sqrt.(xw.^2 .+ yw.^2))
	uw = hw .* -yw ./ sqrt.(xw.^2 .+ yw.^2)
	vw = hw .* xw ./ sqrt.(xw.^2 .+ yw.^2)
	uwn = uw ./ sqrt.(uw.^2 + vw.^2)
	vwn = vw ./ sqrt.(uw.^2 + vw.^2)
	
md"### A. Wind Field due to a Point Vortex

The speed of the wind $v$ at a distance $r$ from the point vortex is given by

$$v = \frac{\Gamma}{2\pi r}$$

In the plot below, I do not plot the vorticity $\zeta$, but the magnitude of the horizontal wind $v$.  This is because in the point vortex model, the vorticity is zero everywhere in the domain except at the point vortex itself.

The wind field of a vortex does not evolve with time, and therefore only the theoretical steady state is plotted here."
	
end

# ╔═╡ 83097010-4c69-11eb-2301-4b6124f1bd07
 md"plot limits : $(@bind lim PlutoUI.Slider(1:0.1:2))"

# ╔═╡ 8c189532-4c69-11eb-15b9-4b52cb5bee40
begin
	pplt.close(); f,axs = pplt.subplots(axwidth=2)
	
	cb = axs[1].contourf(x,y[:],v',levels=0:10,extend="max",cmap="viridis")
	axs[1].format(xlim=(-lim,lim),ylim=(-lim,lim),xlabel="x",ylabel="y")
	axs[1].format(suptitle="Magnitude of Wind")
	axs[1].quiver(xw,yw[:],uwn',vwn')
	
	f.colorbar(cb,loc="r")
	if !isdir(plotsdir()); mkpath(plotsdir()) end
	f.savefig(plotsdir("pointvortex-initial.png"),transparent=false,dpi=200)
	load(plotsdir("pointvortex-initial.png"))
end

# ╔═╡ f9ed505c-4cbb-11eb-1746-6bbedc5aa269
md"Note that in cases when there is a grid point at (0,0), the wind is undefined."

# ╔═╡ 984cba72-4c69-11eb-2337-57a2fba0b874
md"
### B. Advection of Point Vortices about each other

Of interest is the fact that when point vortices interact with one another, they advect each other, but do not impact the strength of other point vortices.  The wind field due to multiple interacting point vortices is simply the vector sum of the wind fields of the individual point vortices.  The position of the vortices evolve with the wind field.

Here, we investigate the stability of systems involving different combinations of point vortices (in both number of vortices and in their relative strengths)
"

# ╔═╡ Cell order:
# ╟─0a3bc67c-4c69-11eb-1dbf-37273535c61c
# ╟─4fe6a84c-4c69-11eb-295f-b7052b606cc3
# ╟─56ce9c32-4c69-11eb-2ebd-79c3f4f19990
# ╟─a26602e8-4c69-11eb-3bb9-21a107f594f0
# ╟─72969256-4c69-11eb-0c43-fba850561b78
# ╟─7a1407e0-4c69-11eb-3d67-b56be81e8f83
# ╟─83097010-4c69-11eb-2301-4b6124f1bd07
# ╟─8c189532-4c69-11eb-15b9-4b52cb5bee40
# ╟─f9ed505c-4cbb-11eb-1746-6bbedc5aa269
# ╟─984cba72-4c69-11eb-2337-57a2fba0b874
