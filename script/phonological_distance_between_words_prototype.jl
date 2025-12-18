# Phonological Distance Between Words
# takes a list of words and a phonological distance dictionary
# and computes pairwise phonological distances using dynamic programming alignment
# note: it creates word1, word2, dist, but not word2, word1 dist also. beware.
# Author: Péter Rácz

# -- Setup -- #

using CSV, DataFrames

# -- Define functions -- #

# Align two words using dynamic programming based on phonological distances.
function align_words(s1::String, s2::String, distance_dict::Dict{Tuple{String,String}, Float64}, gap_penalty::Float64=1.0)
    # Convert strings to vectors of characters/segments
    seq1 = collect(s1)
    seq2 = collect(s2)
    
    n = length(seq1)
    m = length(seq2)
    
    # Initialize DP matrix
    M = zeros(Float64, n+1, m+1)
    
    # Fill first row and column (gaps from start)
    for i in 1:n
        M[i+1, 1] = M[i, 1] + gap_penalty
    end
    for j in 1:m
        M[1, j+1] = M[1, j] + gap_penalty
    end
    
    # Fill rest of matrix
    for i in 1:n
        for j in 1:m
            seg1 = string(seq1[i])
            seg2 = string(seq2[j])
            
            match_cost = M[i, j] + get(distance_dict, (seg1, seg2), 1.0)
            delete_cost = M[i, j+1] + gap_penalty
            insert_cost = M[i+1, j] + gap_penalty
            
            M[i+1, j+1] = min(match_cost, delete_cost, insert_cost)
        end
    end
    
    # Traceback to get alignment
    alignment1 = String[]
    alignment2 = String[]
    distances = Float64[]
    
    i, j = n, m
    while i > 0 || j > 0
        if i > 0 && j > 0
            seg1 = string(seq1[i])
            seg2 = string(seq2[j])
            match_cost = M[i, j] + get(distance_dict, (seg1, seg2), 1.0)
            
            # Use approximate comparison for floats
            if abs(M[i+1, j+1] - match_cost) < 1e-9
                # Match/substitution
                pushfirst!(alignment1, seg1)
                pushfirst!(alignment2, seg2)
                pushfirst!(distances, get(distance_dict, (seg1, seg2), 1.0))
                i -= 1
                j -= 1
            elseif abs(M[i+1, j+1] - (M[i, j+1] + gap_penalty)) < 1e-9
                # Deletion from seq1
                pushfirst!(alignment1, seg1)
                pushfirst!(alignment2, "-")
                pushfirst!(distances, gap_penalty)
                i -= 1
            else
                # Insertion to seq1
                pushfirst!(alignment1, "-")
                pushfirst!(alignment2, seg2)
                pushfirst!(distances, gap_penalty)
                j -= 1
            end
        elseif i > 0
            pushfirst!(alignment1, string(seq1[i]))
            pushfirst!(alignment2, "-")
            pushfirst!(distances, gap_penalty)
            i -= 1
        else
            pushfirst!(alignment1, "-")
            pushfirst!(alignment2, string(seq2[j]))
            pushfirst!(distances, gap_penalty)
            j -= 1
        end
    end
    
    total_dist = sum(distances)
    alignment_length = length(distances)
    
    return (
        segment1 = alignment1,
        segment2 = alignment2,
        dist = distances,
        phon_dist = total_dist,
        length = alignment_length
    )
end

# Align multiple word pairs and return results
function align_word_pairs(words::Vector{String}, distance_dict::Dict)
    results = []
    
    # Generate all unique pairs
    for i in 1:length(words)
        for j in (i+1):length(words)
            word1 = words[i]
            word2 = words[j]
            
            alignment = align_words(word1, word2, distance_dict)
            
            push!(results, (
                word1 = word1,
                word2 = word2,
                segment1 = alignment.segment1,
                segment2 = alignment.segment2,
                dist = alignment.dist,
                phon_dist = alignment.phon_dist,
                length = alignment.length
            ))
        end
    end
    
    return results
end


# -- Execute -- #

# load distance dictionary
distance_dict = CSV.read("out/segment_similarity.tsv", DataFrame; delim='\t')

# load df w/ words
words = CSV.read("out/forms.tsv", DataFrame; delim='\t')

# words is the unique values of the lemma column in words as String
words = unique(String.(words.lemma))

# convert similarity to distance (dist = 1 - similarity)
distance_dict.dist = 1.0 .- distance_dict.similarity
distance_dict = Dict((String(row.segment1), String(row.segment2)) => row.dist for row in eachrow(distance_dict))

results = align_word_pairs(words, distance_dict)

# Save results
results_df = DataFrame(results)

# Create symmetric distance matrix for later matrix algebra
# dist1: original word1 -> word2
dist1 = select(results_df, :word1, :word2, :phon_dist)

# dist2: mirror image word2 -> word1
dist2 = DataFrame(
    word1 = results_df.word2,
    word2 = results_df.word1,
    phon_dist = results_df.phon_dist
)

# dist3: self-distances (word -> word = 0)
all_words = unique(vcat(results_df.word1, results_df.word2))
dist3 = DataFrame(
    word1 = all_words,
    word2 = all_words,
    phon_dist = fill(0.0, length(all_words))
)

# combine all
dist_symmetric = vcat(dist1, dist2, dist3)

CSV.write("out/aligned_word_pairs_phonological_distance.tsv.gz", dist_symmetric; delim='\t', compress=true)
