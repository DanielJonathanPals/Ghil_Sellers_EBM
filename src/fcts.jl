using LinearAlgebra, Interpolations, DifferentialEquations

"""
    Grid{T}

1D discretization Grid 

# Initialization 

    Grid(x::AbstractRange{T})

* `x`: coordinate axis 

# Fields

* `N`: Number of grid points 
* `x`: Coordinate axis
* `Δx`: Spacing between grid points 
"""
struct Grid{T}
    N 
    x::AbstractVector{T}
    Δx::T
end 

function Grid(x::AbstractRange{T}) where T
    N = length(x)
    Δx = abs(x[2] - x[1])
    return Grid{T}(N, x, Δx)
end


"""
    NeumannFD{T}

1D - Finite Difference Scheme matrix with Neumann boundary conditions

# Initialization

    NeumannFD(T::DataType, n::Integer, Δx::Number=1)

* `T`: precision used (e.g. `Float64`)
* `n`: number of grid points 
* `Δx`: spacing between grid points

    NeumannFD(grid::Grid{T})

# Usage 

```julia
g = Grid(1:1:10)
∂ₓ = NeumannFD(g)
dfdx = ∂ₓ(f)
```
"""
struct NeumannFD{T} 
    M::AbstractArray{T,2}
end 

function NeumannFD(T::DataType, n::Integer, Δx::Number=1)
    M = diagm(-1=>(-1*ones(T, n-1)),1=>ones(T, n-1))
    M[1,2] = T(0)
    M[n,n-1] = T(0)
    M ./= T(2*Δx)
    NeumannFD(M)
end 

(FD::NeumannFD{T})(x::AbstractVector{T}) where T = FD.M * x

NeumannFD(grid::Grid{T}) where T = NeumannFD(T, grid.N, grid.Δx)


"""
    ContinousGhilSellersParameters{T}

Holds all parameters needed for the Ghil Sellers 1D EBM. The model is given in CGS units. 
Uses the tabled values from the Ghil Paper and interpolates on a new grid. 

# Initialization

    ContinousGhilSellersParameters(g::Grid; μ=1., m=0.5, order="4th")

* `g::Grid{T}` 
* `ϕ::AbstractVector{T}` latitude vector
* `∂ₓ::NeumannFD{T}` FD scheme
* `T_0::AbstractVector{T}` initial condition given in the Ghil paper
* `μ::T` modifier for the incoming solar irradiance, can be used to produce a bifuraction diagram
* `C::AbstractVector{T}` heat capacity
* `Q::AbstractVector{T}` solar irradiance
* `b::AbstractVector{T}` empirical albedo coefficient
* `z::AbstractVector{T}` average surface elevation 
* `k_1::AbstractVector{T}` sensible heat flux coefficient
* `k_2::AbstractVector{T}`coeffcient of latent heat flux in the atmosphere 
* `c_1::T = 0.009` empirical albedo coefficient
* `c_2::T = 0.0065` empirical albedo coefficient
* `c_3::T = 1.9e-15` empirical emissivity coeffcient
* `c_4::T = 6.105*0.75*exp(19.6)` latent heat flux coefficient 
* `c_5::T = 5350` latent heat flux coefficient 
* `σ::T = 1.356e-12` Stefan Boltzmann constant (in CGS system)
* `m::T = 0.5` Atmospheric attenuation coefficient (0.5 for present conditions)
* `T_m::T = 283.16` reference temperature 
* `α_max::T = 0.85` maximum albedo value, 0.65 in Bodai et al, 0.85 in Ghil and Sellers

* `order` is the order of the finite difference scheme, `2nd` and `4th` are currently supported
"""
Base.@kwdef struct ContinousGhilSellersParameters{T}
    g::Grid{T}
    ϕ::AbstractVector{T}
    ∂ₓ::NeumannFD{T}
    T_0::AbstractVector{T}
    μ::T
    C::AbstractVector{T}
    Q::AbstractVector{T}
    b::AbstractVector{T}
    z::AbstractVector{T}
    k_1::AbstractVector{T}
    k_2::AbstractVector{T}
    c_1::T = 0.009
    c_2::T = 0.0065
    c_3::T = 1.9e-15
    c_4::T = 6.105*0.75*exp(19.6)
    c_5::T = 5350. 
    σ::T = 1.356e-12 # Stefan Boltzmann constant (in CGS system)
    m::T = 0.5 # atmospheric attenuation coefficient (0.5 for present conditions)
    T_m::T = 283.16;
    α_max::T = 0.65;
end 

function ContinousGhilSellersParameters(g::Grid; μ=1., m=0.5)
    T_0, C, Q, b, z, k_1, k_2 = load_interpolated_parameters(g)  
    ϕ = (π.*g.x)./2
    ∂ₓ = NeumannFD(g)  
    ContinousGhilSellersParameters(g=g, ϕ=ϕ, ∂ₓ=∂ₓ, T_0=T_0, μ=μ, C=C, Q=Q, b=b, z=z, k_1=k_1, k_2=k_2, m=m)
end 

function load_interpolated_parameters(g::Grid) 
     # Input Parameters from Sellers Paper, estimated from observational data for one hemisphere
     colat_1 = (-90.:10.:90.)./90. # used for C,Q,T_0
     colat_2 = (-85.:10.:85.)./90. # used for b,z,k_1,k_2
     
     # An intial condition given in the paper
     T_0 = [247.3625, 252.0740, 262.5715, 271.2980, 278.9325, 285.7530, 291.4090, 296.0815, 298.7815, 299.3510];
     
     # Combined air-land-sea effective heat capacity [cal/cm^2/K] / on grid
     C = [500, 1000, 1500, 4725, 5625, 5812, 5813, 5625, 6000, 5625] 
 
     # High-frequency solar irradiance [cal/cm^2/sec]
     Q = [0.426, 0.440, 0.484, 0.579, 0.696, 0.804, 0.894, 0.961, 1.003, 1.017]*1e-2;
 
     # first value in the Ghil paper may be a typo, correction suggested by Tomas Bodai 
     b = [2.912, 2.96, 2.934, 2.914, 2.915, 2.868, 2.821, 2.804, 2.805]
 
     z = [1204.5, 820.0, 295.0, 150.5, 193.5, 301.0, 261.0, 133.5, 156.0] #[m]
 
     # Eddy diffusivity coefficients; k1*T'x: sensible heat flux,
     # k2*g(T)*T'x: latent heat flux
     k_1 = [0.47113, 0.61988, 1.19933, 1.50214, 1.51063, 1.69562, 2.02342, 3.20611, 4.80401]*1e-5; # [cal/K/cm^2/sec]
 
     # first value different from zero to prevent negative values at the pole during extrapolation
     k_2 = [0.3, 0.9314, 1.9772, 3.4348, 4.8316, 3.7359, 0.6903, -2.5401, -10.5975]*1e-2; # [cal/dyn/sec]
     
     # southern hemisphere gets the same constants 
     T_0 = [T_0; reverse(T_0[1:end-1])]
     C = [C; reverse(C[1:end-1])]
     Q = [Q; reverse(Q[1:end-1])]
     b = [b; reverse(b)]
     z = [z; reverse(z)]
     k_1 = [k_1; reverse(k_1)]
     k_2 = [k_2; reverse(k_2)]
     
     #flat boundary conditions correspond due to the Neumann BC
     T_int = CubicSplineInterpolation(colat_1, T_0, bc=Flat(OnGrid()))
     C_int = CubicSplineInterpolation(colat_1, C, bc=Flat(OnGrid()))
     Q_int = CubicSplineInterpolation(colat_1, Q, bc=Flat(OnGrid()))
     b_int = CubicSplineInterpolation(colat_2, b, bc=Flat(OnCell()))
     z_int = CubicSplineInterpolation(colat_2, z, bc=Flat(OnCell()))
     k_1_int = CubicSplineInterpolation(colat_2, k_1, bc=Flat(OnCell()))
     k_2_int = CubicSplineInterpolation(colat_2, k_2, bc=Flat(OnCell()))
     
     T_0 = T_int.(g.x)
     C = C_int.(g.x)
     Q = Q_int.(g.x)
     b = b_int.(g.x)
     z = z_int.(g.x)
     k_1 = k_1_int.(g.x)
     k_2 = k_2_int.(g.x)

     return T_0, C, Q, b, z, k_1, k_2
end 


"""
    R_i(T,p::ContinousGhilSellersParameters)

Incoming radiation at discretized coordinate `x` and temperature `T` 
"""
R_i(T,p::ContinousGhilSellersParameters) = p.Q .* (1 .- α(T,p))


"""
    α(T,p::ContinousGhilSellersParameters)

Albedo at discretized coordinate `x` and temperature `T`. Albedo is cutoff at minimum 0.25 and maximum α_max
"""
α(T,p::ContinousGhilSellersParameters) = cutoff_function.(p.b - p.c_1 .* (p.T_m .+ min.(0, T - p.c_2 .* p.z .- p.T_m)), p.α_max)


function cutoff_function(x, max_val) 
    if x <= 0.25
        return 0.25
    elseif x >= max_val
        return max_val
    else 
        return x
    end 
end


"""
    R_o(T,p::ContinousGhilSellersParameters)

Outgoing radiation according to Stefan-Boltzmann law with emissivity reduced to account for greenhouse gases
"""
R_o(T,p::ContinousGhilSellersParameters) = p.σ .* T.^4 .* c(T,p)  

"""
    c(T,p::ContinousGhilSellersParameters)

Emissivity coefficient accounting for greenhouse gases
"""
c(T,p::ContinousGhilSellersParameters) = 1 .- p.m .* tanh.(p.c_3 .* T.^6)


"""
    k(T,p::ContinousGhilSellersParameters)

Heat flux, first term sensible heat flux, second term latent heat flux
"""
k(T,p::ContinousGhilSellersParameters) = p.k_1 + p.k_2 .* g(T,p)
k(T,i::Integer,p::ContinousGhilSellersParameters) = p.k_1[i] + p.k_2[i] * g(T,p)

"""
    g(T,p::ContinousGhilSellersParameters)

Latent heat flux contribution
"""
g(T,p::ContinousGhilSellersParameters) = p.c_4 .* exp.(-p.c_5 ./ T) ./ (T.^2)


"""
    D(T,p::ContinousGhilSellersParameters)  

1D Diffusion with Neumann boundary conditions applied
"""
function D(T,p::ContinousGhilSellersParameters)  
    A = ((2/π)^2) .* (1 ./ cos.(p.ϕ))
    
    D = A .* p.∂ₓ( cos.(p.ϕ) .* k(T,p) .* p.∂ₓ(T))
    D[1] = 2 * A[1] * cos(p.ϕ[1]) * k(T[1],1,p) * (T[2] - T[1])/(p.g.Δx^2)
    D[end] = 2 * A[end] * cos(p.ϕ[end]) * k(T[end],p.g.N,p) * (T[end-1] - T[end])/(p.g.Δx^2)

    return D
end


"""
    ghilsellers_ebm!(du,u,p,t)   

RHS of the PDE, to be used with `DifferentialEquations.jl`
"""
function ghilsellers_ebm!(du,u,p,t)   
    du .= (R_i(u,p) - R_o(u,p) + D(u,p)) ./ p.C  
end