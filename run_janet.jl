#!/usr/bin/env julia

# JANET: Joint Alignment and Nonparametric Estimation Toolkit
# Run the full pipeline from feature matrix to word distances

using Printf

# --- Colours for output --- #

const RED = "\033[0;31m"
const GREEN = "\033[0;32m"
const YELLOW = "\033[0;33m"
const NC = "\033[0m"

error_exit(msg) = (println(stderr, "$(RED)Error:$(NC) $msg"); exit(1))
warn(msg) = println(stderr, "$(YELLOW)Warning:$(NC) $msg")
info(msg) = println("$(GREEN)→$(NC) $msg")

# --- Usage --- #

function usage()
    println("""
    Usage: julia run_janet.jl <feature_matrix.tsv> <word_list.tsv>

    Arguments:
      feature_matrix.tsv  Tab-separated file with segment features
                          First column must be named 'segment'
                          See README.md for format details

      word_list.tsv       Tab-separated file with words (one per line)
                          If no header present, 'lemma' will be added
                          Words must use same alphabet as feature matrix

    Output:
      out/segment_similarity.tsv                            Pairwise segment similarities
      out/aligned_word_pairs_phonological_distance.tsv.gz   Pairwise word distances
    """)
    exit(1)
end

# --- Validation functions --- #

function validate_feature_matrix(path::String)
    info("Validating feature matrix...")
    
    isfile(path) || error_exit("Feature matrix not found: $path")
    
    lines = readlines(path)
    isempty(lines) && error_exit("Feature matrix is empty")
    
    # Check header
    header = split(lines[1], '\t')
    header[1] == "segment" || error_exit("Feature matrix first column must be named 'segment', found: '$(header[1])'")
    
    # Extract segments
    segments = String[]
    for (i, line) in enumerate(lines[2:end])
        fields = split(line, '\t')
        isempty(fields) && continue
        seg = fields[1]
        
        # Check single character
        length(seg) == 1 || error_exit("Segment '$seg' (line $(i+1)) is not a single character")
        
        # Check uniqueness
        seg in segments && error_exit("Duplicate segment '$seg' (line $(i+1))")
        
        push!(segments, seg)
    end
    
    info("Feature matrix OK: $(length(segments)) segments")
    return segments
end

function validate_word_list(path::String, segments::Vector{String})
    info("Validating word list...")
    
    isfile(path) || error_exit("Word list not found: $path")
    
    lines = readlines(path)
    isempty(lines) && error_exit("Word list is empty")
    
    # Check if header present
    has_header = strip(lines[1]) == "lemma"
    
    if !has_header
        warn("No 'lemma' header found, will add it")
    end
    
    # Build segment set for fast lookup
    segment_set = Set(segments)
    
    # Validate words
    start_line = has_header ? 2 : 1
    words = String[]
    
    for (i, line) in enumerate(lines[start_line:end])
        word = strip(line)
        isempty(word) && continue
        
        # Check each character
        for char in word
            char_str = string(char)
            if !(char_str in segment_set)
                error_exit("Unknown segment '$char_str' in word '$word' (line $(i + start_line - 1)). Not in feature matrix.")
            end
        end
        
        push!(words, word)
    end
    
    info("Word list OK: $(length(words)) words")
    return words, has_header
end

# --- Main --- #

function main()
    # Check arguments
    length(ARGS) == 2 || usage()
    
    feature_matrix_path = ARGS[1]
    word_list_path = ARGS[2]
    
    # Get script directory
    script_dir = dirname(@__FILE__)
    if isempty(script_dir)
        script_dir = "."
    end
    
    # Check required files exist
    gen_script = joinpath(script_dir, "script", "generate_segmental_distances.jl")
    align_script = joinpath(script_dir, "script", "phonological_distance_between_words.jl")
    
    isfile(gen_script) || error_exit("Missing: $gen_script")
    isfile(align_script) || error_exit("Missing: $align_script")
    
    # Create directories
    source_dir = joinpath(script_dir, "source")
    out_dir = joinpath(script_dir, "out")
    mkpath(source_dir)
    mkpath(out_dir)
    
    # Validate inputs
    segments = validate_feature_matrix(feature_matrix_path)
    words, has_header = validate_word_list(word_list_path, segments)
    
    # Copy feature matrix
    info("Copying input files to source/...")
    cp(feature_matrix_path, joinpath(source_dir, "features.tsv"); force=true)
    
    # Copy word list (adding header if needed)
    forms_path = joinpath(source_dir, "forms.tsv")
    open(forms_path, "w") do f
        println(f, "lemma")
        for word in words
            println(f, word)
        end
    end
    
    # Run pipeline
    info("Generating segment similarities...")
    cd(script_dir)
    run(`julia $gen_script`)
    
    info("Computing word distances...")
    run(`julia $align_script`)
    
    # Done
    println()
    println("$(GREEN)Done!$(NC)")
    println()
    println("Output files:")
    println("  out/segment_similarity.tsv")
    println("  out/aligned_word_pairs_phonological_distance.tsv.gz")
    println()
    println("To use with kernel ridge regression, see script/krr.R")
end

main()
