export cog
export rog
export dihedral
export superimposition
export transform
export transform!

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
    d = det(s.V * s.U') < 0.0 ? -1.0 : 1.0
    m = SMatrix{3,3,T,9}(diagm([1, 1, d]))

    rot = s.V * m * s.U'
    tra = cog2 - rot*cog1

    return Transformation(rot, tra)
end

function transform!(r::Vector{SVector{3,T}}, t::Transformation{T}) where {T}
    for i in eachindex(r)
        r[i] = t.rotation * r[i] + t.translation
    end
end

function transform(r::Vector{SVector{3,T}}, t::Transformation{T}) where {T}
    return Ref(t.rotation) .* r .+ Ref(t.translation)
end