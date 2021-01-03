### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ 8326524e-4bff-11eb-3bf7-f7bb0761da5a
md"# 1. Models of Hurricane Vorticity
		
In this project, I investigate idealized vortex models of hurricanes and their evolution using the Julia package `GeophysicalFlows.jl`, developed by Dr. Wagner and Dr. Constantinou.  This package uses Fourier-based pseudospectral methods to solve fluid dynamical problems in the atmosphere and ocean.  I use the `TwoDNavierStokes` module that is provided by this package, which simply solves the two-dimensional vorticity equation.
		
Many thanks to Argel Ramirez for inspiring me to try out using the package for hurricane-related projects.  This work is not original, but attempts to replicate work and results of Hendricks et al. (2009), \"Life Cycles of Hurricane-Like Vorticity Rings\", and Hendricks & Schubert. (2009), \"Transport and mixing in idealized barotropic hurricane-like vortices\".
		
The hurricane vorticity models I investigated are listed below:
* Point Vortex Model (Notebook 2)
* Elliptical Vortex Model (Notebook 3)
* Rankine Vortex Model (Notebook 4)
* Vortex Ring Model (Notebook 5+)"

# ╔═╡ Cell order:
# ╟─8326524e-4bff-11eb-3bf7-f7bb0761da5a
