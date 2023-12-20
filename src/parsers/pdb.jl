export PDBLine

struct PDBLine
    atom_serial::Int     # line[7:11]
    atom_name::String    # line[13:16]
    residue_name::String # line[18:21]
    chain_id::String     # line[22]
    residue_serial::Int  # line[23:26]
    coor_x::Float64      # line[31:38]
    coor_y::Float64      # line[39:46]
    coor_z::Float64      # line[47:54]
    occupancy::Float64   # line[55:60]
    tempfactor::Float64  # line[61:66]
    segment_id::String   # line[67:76]
    element_name::String # line[77:78]
    charge::Float64      # line[79:80]
end

function PDBLine(line::AbstractString)
    @views atom_serial    = tryparse(Int, line[7:11])
    @views atom_name      = strip(line[13:16])
    @views residue_name   = strip(line[18:21])
    @views chain_id       = line[22:22]
    @views residue_serial = tryparse(Int, line[23:26])
    @views coor_x         = tryparse(Float64, line[31:38])
    @views coor_y         = tryparse(Float64, line[39:46])
    @views coor_z         = tryparse(Float64, line[47:54])
    @views occupancy      = tryparse(Float64, line[55:60])
    @views tempfactor     = tryparse(Float64, line[61:66])
    @views segment_id     = strip(line[67:76])
    @views element_name   = strip(line[77:78])
    @views charge         = tryparse(Float64, line[79:80])

    isnothing(residue_serial) && throw(ArgumentError("Failed parsing 'residue serial' of PDB ATOM record: $line"))
    isnothing(coor_x)         && throw(ArgumentError("Failed parsing 'x' of PDB ATOM record: $line"))
    isnothing(coor_y)         && throw(ArgumentError("Failed parsing 'y' of PDB ATOM record: $line"))
    isnothing(coor_z)         && throw(ArgumentError("Failed parsing 'z' of PDB ATOM record: $line"))

    return PDBLine(
        if isnothing(atom_serial) 0 else atom_serial end,
        atom_name,
        residue_name,
        chain_id,
        residue_serial,
        coor_x,
        coor_y,
        coor_z,
        if isnothing(occupancy) 0.0 else occupancy end,
        if isnothing(tempfactor) 0.0 else tempfactor end,
        segment_id,
        element_name,
        if isnothing(charge) 0.0 else charge end,
    )
end

function Base.show(io::IO, pdbline::PDBLine)
    @printf(
        io,
        "ATOM  %5d %4s %4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%10s%2s%2d \n",
        pdbline.atom_serial,
        pdbline.atom_name,
        rpad(pdbline.residue_name, 4),
        pdbline.chain_id,
        pdbline.residue_serial,
        pdbline.coor_x,
        pdbline.coor_y,
        pdbline.coor_z,
        pdbline.occupancy,
        pdbline.tempfactor,
        pdbline.segment_id,
        pdbline.element_name,
        Int(round(pdbline.charge)),
    )
end