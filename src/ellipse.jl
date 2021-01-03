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

	nx = length(x); ny = length(y); ζ = zeros(nx,ny)
	ygrid = reshape(y,1,:)
	r = sqrt.(x.^2 .+ ygrid.^2)
	α = fα.(ϵ,x,ygrid)
	riα = ri * α
	r0α = ri * α * r0

	for iy in 1 : ny, ix in 1 : nx

		rii = riα[ix,iy]
		r0i = r0α[ix,iy]
		r_i = r[ix,iy]

		if r_i < rii
			ζ[ix,iy] = 1
		elseif r_i < r0i
			ζ[ix,iy] = 1 - fλ(r_i,rii,r0i,λ)
		end

	end

	return ζ * ζ0

end
