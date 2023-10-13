export cog
export rog
export com
export dihedral
export superimposition
export transform
export transform!
export rmsd
export cross_matrix
export axis_rotation
export axis_rotationd
export x_axis_rotation
export x_axis_rotationd
export y_axis_rotation
export y_axis_rotationd
export z_axis_rotation
export z_axis_rotationd
export random_rotation

struct Transformation{T}
    rotation::SMatrix{3,3,T,9}
    translation::SVector{3,T}
end

cog(r::Vector{SVector{3,T}}) where {T} = mean(r)

function rog(r::Vector{SVector{3,T}}) where {T}
    r .-= Ref(cog(r))
    s = sum(v -> dot(v,v), r)
    return sqrt(s / length(r))
end

function com(r::Vector{SVector{3,T}}, m::Vector{T}) where {T}
    acc = zero(SVector{3,T})
    for i in eachindex(r, m)
        acc += r[i] * m[i]
    end

    return acc / sum(m)
end

function Base.angle(v1::SVector{3,T}, v2::SVector{3,T}) where {T}
    n1 = norm(v1)
    n2 = norm(v2)
    return acosd(clamp(dot(v1,v2)/n1/n2, -1, 1))
end

function Base.angle(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where{T}
    v21 = v1 - v2
    v23 = v3 - v2
    return angle(v21, v23)
end

function dihedral(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}, v4::SVector{3,T}) where {T}
    v12 = v2 - v1
    v23 = v3 - v2
    v34 = v4 - v3

    c123 = cross(v12, v23)
    c234 = cross(v23, v34)
    dihe = angle(c123,c234)

    c1234 = cross(c123, c234)
    judge = dot(c1234, v23)

    return judge ≥ 0 ? dihe : -dihe
end

function superimposition(r1::Vector{SVector{3,T}}, r2::Vector{SVector{3,T}}) where {T}
    cog1 = cog(r1)
    cog2 = cog(r2)

    v1 = r1 .- Ref(cog1)
    v2 = r2 .- Ref(cog2)

    r = zero(SMatrix{3,3,T,9})
    for i in eachindex(v1, v2)
        r += v1[i] * v2[i]'
    end

    s = svd(r)
    d = det(s.V * s.U') < zero(T) ? -one(T) : one(T)
    m = SMatrix{3,3,T,9}(diagm([one(T), one(T), d]))

    rot = s.V * m * s.U'
    tra = cog2 - rot*cog1

    return Transformation(rot, tra)
end

function transform!(r::Vector{SVector{3,T}}, t::Transformation{T}) where {T}
    @inbounds @simd for i in eachindex(r)
        r[i] = t.rotation * r[i] + t.translation
    end
end

function transform(r::Vector{SVector{3,T}}, t::Transformation{T}) where {T}
    return Ref(t.rotation) .* r .+ Ref(t.translation)
end

function transform(r::SVector{3,T}, t::Transformation{T}) where {T}
    return t.rotation * r + t.translation
end

function rmsd(r1::Vector{SVector{3,T}}, r2::Vector{SVector{3,T}}) where {T}
    n = length(r1)
    t = superimposition(r1, r2)

    s = zero(T)
    @inbounds for i in eachindex(r1,r2)
        v = r2[i] - transform(r1[i], t)
        s += dot(v,v)
    end

    return sqrt(s/n)
end

function cross_matrix(v::SVector{3,T}) where {T}
    return SMatrix{3,3,T,9}([
        zero(T), -v[3], v[2],
        v[3], zero(T), -v[1],
        -v[2], v[1], zero(T),
    ])
end

function axis_rotation(axis::SVector{3,T}, θ::T) where {T}
    m1 = axis * axis'
    m2 = cross_matrix(axis)

    c = cos(θ)
    s = sin(θ)

    return c * SMatrix{3,3,T,9}(I) + (one(T) - c) * m1 + s * m2
end

axis_rotationd(axis::SVector{3,T}, θ::T) where {T} = axis_rotation(axis, deg2rad(θ))

function x_axis_rotation(θ::T) where {T}
    c = cos(θ)
    s = sin(θ)

    return SMatrix{3,3,T,9}([1 0 0; 0 c -s; 0 s c;])
end

x_axis_rotationd(θ::T) where {T} = x_axis_rotation(deg2rad(θ))

function y_axis_rotation(θ::T) where {T}
    c = cos(θ)
    s = sin(θ)

    return SMatrix{3,3,T,9}([c 0 s; 0 1 0; -s 0 c])
end

y_axis_rotationd(θ::T) where {T} = y_axis_rotation(deg2rad(θ))

function z_axis_rotation(θ::T) where {T}
    c = cos(θ)
    s = sin(θ)

    return SMatrix{3,3,T,9}([c -s 0; s c 0; 0 0 1])
end

z_axis_rotationd(θ::T) where {T} = z_axis_rotation(deg2rad(θ))

function random_rotation(T::Type{<:AbstractFloat}; rng=Random.default_rng())
    v = normalize(rand(rng, SVector{3,T}) .- T(0.5))
    θ = T(2π * rand(T))

    return axis_rotation(v, θ)
end

random_rotation(; rng=Random.default_rng()) = random_rotation(Float64; rng)