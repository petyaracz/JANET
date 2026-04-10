using CSV, DataFrames

# -- setup -- #

# Core data structure: a simple matrix representation
struct SegmentInventory
    segments::Vector{String}      # segment labels: ["a", "e", "i", ...]
    features::Vector{String}       # feature names: ["labial", "coronal", ...]
    matrix::Matrix{Union{Int, Missing}}  # binary feature matrix with possible missing values
end

# A natural class is just a set of feature specifications
struct FeatureSpec
    features::Dict{String, Int}    # e.g., {"labial" => 1, "coronal" => 0}
end

# A natural class picks out a subset of segments
struct NaturalClass
    spec::FeatureSpec
    members::Set{String}           # which segments belong to this class
end

# -- functions -- #


"""
Find which segments match a feature specification.
A segment matches if it has the specified value for all specified features.
Segments with missing values for a feature are treated as not matching that feature specification.
"""
function find_members(spec::FeatureSpec, inventory::SegmentInventory)
    members = Set{String}()
    
    for (i, seg) in enumerate(inventory.segments)
        matches = true
        for (feat, val) in spec.features
            feat_idx = findfirst(==(feat), inventory.features)
            seg_val = inventory.matrix[i, feat_idx]
            
            # If the segment has no value (missing) for this feature, it doesn't match
            if ismissing(seg_val) || seg_val != val
                matches = false
                break
            end
        end
        if matches
            push!(members, seg)
        end
    end
    
    return members
end

"""
Remove redundant natural classes.
A natural class is redundant if there exists another natural class with:
1. The same members
2. More feature specifications (more specific)
"""
function remove_redundant_classes(classes::Vector{NaturalClass})
    non_redundant = NaturalClass[]
    
    for nc in classes
        # Check if there's a more specific class with the same members
        is_redundant = false
        for other_nc in classes
            # If other class has same members but more features specified, nc is redundant
            if nc.members == other_nc.members && 
               length(other_nc.spec.features) > length(nc.spec.features)
                is_redundant = true
                break
            end
        end
        
        if !is_redundant
            push!(non_redundant, nc)
        end
    end
    
    return non_redundant
end

"""
Generate all natural classes from an inventory using intersection closure.

Instead of enumerating all 3^n feature specifications (which is infeasible for
large feature sets), this algorithm computes all natural-class member sets
directly from the segment inventory.

For each feature and value (0 or 1) we compute a "basic cut": the set of
segments that carry that feature value.  The natural classes are exactly all
non-empty intersections of basic cuts (the intersection closure).  The closure
is bounded by 2^|segments|, which is tractable even for 27+ features.

For each member set in the closure we derive the canonical FeatureSpec: the most
specific spec whose members equal that set (i.e. every feature on which the set
is uniformly non-missing valued).
"""
function generate_natural_classes(inventory::SegmentInventory)
    n_segs = length(inventory.segments)
    n_feats = length(inventory.features)

    # Build basic cuts: one BitVector per (feature, value) pair.
    basic_cuts = BitVector[]
    for f_idx in 1:n_feats
        for v in [0, 1]
            cut = falses(n_segs)
            for s_idx in 1:n_segs
                val = inventory.matrix[s_idx, f_idx]
                if !ismissing(val) && val == v
                    cut[s_idx] = true
                end
            end
            if any(cut)
                push!(basic_cuts, cut)
            end
        end
    end

    # Compute the intersection closure of all basic cuts.
    # The queue holds newly discovered sets that still need to be intersected
    # with every basic cut.  When the queue empties, no new sets can arise.
    closed = Set{BitVector}(basic_cuts)
    queue  = copy(basic_cuts)
    while !isempty(queue)
        C = pop!(queue)
        for B in basic_cuts
            I = C .& B
            if any(I) && I ∉ closed
                push!(closed, I)
                push!(queue, I)
            end
        end
    end

    # For each member set derive its canonical FeatureSpec: include a feature
    # only when every segment in the set carries the same non-missing value.
    classes = NaturalClass[]
    for member_bits in closed
        members = Set{String}(
            inventory.segments[i] for i in 1:n_segs if member_bits[i]
        )

        spec_dict = Dict{String, Int}()
        for f_idx in 1:n_feats
            vals = [inventory.matrix[s_idx, f_idx]
                    for s_idx in 1:n_segs if member_bits[s_idx]]
            non_missing = filter(!ismissing, vals)
            if length(non_missing) == length(vals)
                first_val = non_missing[1]
                if all(==(first_val), non_missing)
                    spec_dict[inventory.features[f_idx]] = first_val
                end
            end
        end

        # Skip classes with an empty spec (no distinguishing features).
        if !isempty(spec_dict)
            push!(classes, NaturalClass(FeatureSpec(spec_dict), members))
        end
    end

    return classes
end

"""
Calculate similarity between two segments using natural classes metric.
Similarity = |shared classes| / (|shared classes| + |non-shared classes|)
"""
function calculate_similarity(seg1::String, seg2::String, 
                             classes::Vector{NaturalClass})
    shared = 0
    non_shared = 0
    
    for nc in classes
        in_seg1 = seg1 ∈ nc.members
        in_seg2 = seg2 ∈ nc.members
        
        if in_seg1 && in_seg2
            shared += 1
        elseif in_seg1 || in_seg2  # exactly one, not both
            non_shared += 1
        end
        # if neither: don't count
    end
    
    if shared + non_shared == 0
        return 0.0  # non-homorganic (no shared features)
    end
    
    return shared / (shared + non_shared)
end

"""
Build a similarity matrix for all segment pairs.
"""
function similarity_matrix(inventory::SegmentInventory)
    classes = generate_natural_classes(inventory)
    n = length(inventory.segments)
    sim_matrix = zeros(Float64, n, n)
    
    for i in 1:n
        for j in 1:n
            seg1 = inventory.segments[i]
            seg2 = inventory.segments[j]
            sim_matrix[i, j] = calculate_similarity(seg1, seg2, classes)
        end
    end
    
    return sim_matrix, classes
end

# Pretty printing helpers:

function print_natural_class(nc::NaturalClass)
    spec_str = join(["$feat=$val" for (feat, val) in nc.spec.features], ", ")
    members_str = join(sort(collect(nc.members)), ", ")
    println("[$spec_str] → {$members_str}")
end

function print_similarity_table(inventory::SegmentInventory, sim_matrix::Matrix{Float64})
    # Create vectors to store the data
    seg1_col = String[]
    seg2_col = String[]
    similarity_col = Float64[]
    
    for (i, seg1) in enumerate(inventory.segments)
        for (j, seg2) in enumerate(inventory.segments)
            sim = round(sim_matrix[i, j], digits=3)
            push!(seg1_col, seg1)
            push!(seg2_col, seg2)
            push!(similarity_col, sim)
        end
    end
    
    # Create and return DataFrame
    return DataFrame(
        segment1 = seg1_col,
        segment2 = seg2_col,
        similarity = similarity_col
    )
end

"""
Add rows for "segment vs nothing" and "nothing vs segment" comparisons.
The similarity between any segment and nothing is always 1.0.
"""
function add_nothing_rows(sim_df::DataFrame)
    segments = unique(sim_df.segment1)
    
    # Create rows for segment vs nothing
    nothing_rows1 = DataFrame(
        segment1 = segments,
        segment2 = fill(" ", length(segments)),
        similarity = fill(1.0, length(segments))
    )
    
    # Create rows for nothing vs segment
    nothing_rows2 = DataFrame(
        segment1 = fill(" ", length(segments)),
        segment2 = segments,
        similarity = fill(1.0, length(segments))
    )
    
    # Also add nothing vs nothing
    nothing_nothing = DataFrame(
        segment1 = [" "],
        segment2 = [" "],
        similarity = [1.0]
    )
    
    # Append all to existing dataframe
    return vcat(sim_df, nothing_rows1, nothing_rows2, nothing_nothing)
end

# -- execute -- #

# Load CSV with custom parsing to handle empty cells as missing
# First, read the header to discover column names
header_df = CSV.read("source/features.tsv", DataFrame; delim='\t', limit=0)
column_names = names(header_df)

# Build types dictionary dynamically:
# - First column (segment) is String
# - All other columns (features) are Union{Int, Missing}
types_dict = Dict{Symbol, Type}()
for col_name in column_names
    col_symbol = Symbol(col_name)
    if col_symbol == :segment
        types_dict[col_symbol] = String
    else
        types_dict[col_symbol] = Union{Int, Missing}
    end
end

# Now read the full file with the dynamically constructed types
feature_matrix = CSV.read("source/features.tsv", DataFrame; 
    delim='\t',
    missingstring="",  # Treat empty strings as missing
    types=types_dict
)

# Extract segments, features, and matrix
segments = feature_matrix.segment
features = names(feature_matrix)[2:end]  # all columns except 'segment'
matrix = Matrix(feature_matrix[:, 2:end])

# Create inventory
inventory = SegmentInventory(segments, features, matrix)

# Generate all natural classes
classes = generate_natural_classes(inventory)

# Print them to see what we got
for nc in classes
    print_natural_class(nc)
end

# Calculate similarity
sim_matrix, _ = similarity_matrix(inventory)
sim_df = print_similarity_table(inventory, sim_matrix)

# Add nothing rows
sim_df = add_nothing_rows(sim_df)

CSV.write("out/segment_similarity.tsv", sim_df; delim='\t')