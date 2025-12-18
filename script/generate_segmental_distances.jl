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
Generate all possible feature specifications (non-empty subsets of features).
For each feature, we can specify it as 1, 0, or leave it unspecified.
"""
function generate_all_feature_specs(feature_names::Vector{String})
    # For each feature: include as +, include as -, or omit
    # This generates 3^n - 1 possibilities (excluding the empty spec)
    specs = FeatureSpec[]
    
    n_features = length(feature_names)
    
    # Generate all combinations: 3^n total (including empty)
    # We'll represent each as a base-3 number
    for i in 1:(3^n_features - 1)  # -1 to exclude empty spec (all omitted)
        spec_dict = Dict{String, Int}()
        
        # Convert i to base-3 representation
        num = i
        for (feat_idx, feat_name) in enumerate(feature_names)
            choice = num % 3  # 0 = omit, 1 = specify as 0, 2 = specify as 1
            num = div(num, 3)
            
            if choice == 1
                spec_dict[feat_name] = 0
            elseif choice == 2
                spec_dict[feat_name] = 1
            end
            # choice == 0: omit (don't add to dict)
        end
        
        # Only add non-empty specs
        if !isempty(spec_dict)
            push!(specs, FeatureSpec(spec_dict))
        end
    end
    
    return specs
end


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
Generate all natural classes from an inventory.
"""
function generate_natural_classes(inventory::SegmentInventory)
    specs = generate_all_feature_specs(inventory.features)
    classes = NaturalClass[]
    
    for spec in specs
        members = find_members(spec, inventory)
        if !isempty(members)  # only keep non-empty classes
            push!(classes, NaturalClass(spec, members))
        end
    end
    
    # Remove redundant classes (proper supersets with same members)
    classes = remove_redundant_classes(classes)
    
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
feature_matrix = CSV.read("source/features.tsv", DataFrame; 
    delim='\t',
    missingstring="",  # Treat empty strings as missing
    types=Dict(
        :segment => String,
        :cons => Union{Int, Missing},
        :son => Union{Int, Missing},
        :cont => Union{Int, Missing},
        :labial => Union{Int, Missing},
        :coronal => Union{Int, Missing},
        :anterior => Union{Int, Missing},
        :dorsal => Union{Int, Missing},
        :lateral => Union{Int, Missing},
        :voice => Union{Int, Missing},
        :delrel => Union{Int, Missing},
        :seg => Union{Int, Missing},
        :long => Union{Int, Missing},
        :open1 => Union{Int, Missing},
        :open2 => Union{Int, Missing}
    )
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