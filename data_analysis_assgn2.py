# DNA Sequence Analysis Program
# Analysis of chr1_GL383518v1_alt sequence

def read_fasta_sequence(filename):
    """
    Read a FASTA file and return the DNA sequence as a string.
    Assumes single sequence in file.
    """
    sequence = ""
    try:
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if not line.startswith('>'):  # Skip header lines
                    sequence += line.upper()  # Convert to uppercase for consistency
        return sequence
    except FileNotFoundError:
        print(f"Error: File {filename} not found.")
        return None
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

def create_reverse_complement(sequence):
    """
    Create the reverse complement of a DNA sequence.
    Watson-Crick-Franklin base pairs: A-T, G-C
    """
    complement_map = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'N': 'N'  # Handle ambiguous nucleotides
    }
    
    # Create complement
    complement = ""
    for base in sequence:
        complement += complement_map.get(base, 'N')
    
    # Reverse the complement
    reverse_complement = complement[::-1]
    return reverse_complement

def analyze_sequence_by_kilobase(sequence):
    """
    Create a nested dictionary with nucleotide counts for each kilobase.
    Returns dict where keys are kilobase positions (0, 1000, 2000, etc.)
    and values are dictionaries with nucleotide counts.
    """
    analysis_dict = {}
    sequence_length = len(sequence)
    
    # Process sequence in 1000 bp chunks
    for start_pos in range(0, sequence_length, 1000):
        end_pos = min(start_pos + 1000, sequence_length)
        chunk = sequence[start_pos:end_pos]
        
        # Count nucleotides in this chunk
        nucleotide_counts = {
            'A': chunk.count('A'),
            'C': chunk.count('C'),
            'G': chunk.count('G'),
            'T': chunk.count('T'),
            'N': chunk.count('N')  # Count ambiguous bases
        }
        
        analysis_dict[start_pos] = nucleotide_counts
    
    return analysis_dict

def print_kilobase_analysis(analysis_dict, sequence_length):
    """
    Print the kilobase analysis results in a readable format.
    """
    print(f"\n{'='*70}")
    print(f"NUCLEOTIDE ANALYSIS BY KILOBASE")
    print(f"Total sequence length: {sequence_length:,} bp")
    print(f"Number of kilobase segments: {len(analysis_dict)}")
    print(f"{'='*70}\n")
    
    # Print header
    print(f"{'Position':<12} {'A':<8} {'C':<8} {'G':<8} {'T':<8} {'N':<8} {'Total':<8}")
    print(f"{'-'*70}")
    
    # Sort keys to ensure proper order
    sorted_positions = sorted(analysis_dict.keys())
    
    # Print each kilobase segment
    for position in sorted_positions:
        counts = analysis_dict[position]
        total = counts['A'] + counts['C'] + counts['G'] + counts['T'] + counts['N']
        kb_number = position // 1000 + 1
        
        print(f"KB {kb_number:<4} ({position:<6}): "
              f"{counts['A']:<8} {counts['C']:<8} {counts['G']:<8} "
              f"{counts['T']:<8} {counts['N']:<8} {total:<8}")
    
    print(f"{'-'*70}\n")

def create_nucleotide_lists(analysis_dict):
    """
    Create lists of nucleotide counts for each kilobase.
    Returns a list of lists, where each inner list contains [A, C, G, T] counts.
    """
    nucleotide_lists = []
    
    # Sort keys to ensure proper order
    sorted_positions = sorted(analysis_dict.keys())
    
    for position in sorted_positions:
        counts = analysis_dict[position]
        # Create list in order [A, C, G, T]
        nucleotide_list = [
            counts['A'],
            counts['C'],
            counts['G'],
            counts['T']
        ]
        nucleotide_lists.append(nucleotide_list)
    
    return nucleotide_lists

# Main analysis program
def main():
    # Assume the FASTA file is named 'chr1_GL383518v1_alt.fa'
    filename = 'chr1_GL383518v1_alt.fa'
    
    print("=== Part 1: Reading DNA Sequence ===")
    sequence = read_fasta_sequence(filename)
    
    if sequence is None:
        print("Cannot proceed without sequence file.")
        return
    
    print(f"Sequence length: {len(sequence)} bp")
    
    # Part 1a: Print 10th letter
    if len(sequence) >= 10:
        print(f"10th letter: {sequence[9]}")  # 0-indexed, so 9th index is 10th letter
    
    # Part 1b: Print 758th letter
    if len(sequence) >= 758:
        print(f"758th letter: {sequence[757]}")  # 0-indexed, so 757th index is 758th letter
    
    print("\n=== Part 2: Reverse Complement Analysis ===")
    reverse_comp = create_reverse_complement(sequence)
    print(f"Reverse complement length: {len(reverse_comp)} bp")
    
    # Part 2a: Print 79th letter of reverse complement
    if len(reverse_comp) >= 79:
        print(f"79th letter of reverse complement: {reverse_comp[78]}")
    
    # Part 2b: Print 500th through 800th letters
    if len(reverse_comp) >= 800:
        substring = reverse_comp[499:800]  # 0-indexed: 499-799 gives positions 500-800
        print(f"Letters 500-800 of reverse complement: {substring}")
    
    print("\n=== Part 3: Kilobase Analysis ===")
    kilobase_dict = analyze_sequence_by_kilobase(sequence)
    print(f"Analysis completed for {len(kilobase_dict)} kilobase segments")
    
    # Print detailed kilobase analysis table
    print_kilobase_analysis(kilobase_dict, len(sequence))
    
    # Show example usage
    print("Example usage:")
    print(f"kilobase_dict[0] = {kilobase_dict[0]}")
    print(f"kilobase_dict[0]['A'] = {kilobase_dict[0]['A']}")
    if 5000 in kilobase_dict:
        print(f"kilobase_dict[5000] = {kilobase_dict[5000]}")
        print(f"kilobase_dict[5000]['A'] = {kilobase_dict[5000]['A']}")
    
    print("\n=== Part 4: Nucleotide List Analysis ===")
    
    # Part 4a: First 1000 base pairs
    first_kb_counts = kilobase_dict.get(0, {'A': 0, 'C': 0, 'G': 0, 'T': 0})
    first_kb_list = [first_kb_counts['A'], first_kb_counts['C'], 
                     first_kb_counts['G'], first_kb_counts['T']]
    print(f"First 1000 bp nucleotide counts [A, C, G, T]: {first_kb_list}")
    
    # Part 4b & 4c: All kilobases
    all_nucleotide_lists = create_nucleotide_lists(kilobase_dict)
    print(f"Created {len(all_nucleotide_lists)} nucleotide count lists")
    
    # Part 4d: Calculate sums
    list_sums = []
    for i, nuc_list in enumerate(all_nucleotide_lists):
        list_sum = sum(nuc_list)
        list_sums.append(list_sum)
        if i < 5:  # Show first 5 examples
            print(f"Kilobase {i+1} sum: {list_sum}")
    
    print(f"\nAll list sums: {list_sums}")
    
    print("\n=== Part 5: Analysis Questions ===")
    print("\n" + "="*70)
    print("QUESTION 1: Are there any lists whose sums are not equal to 1000?")
    print("="*70)
    
    expected_sum = 1000
    unexpected_sums = []
    
    # Check each list sum
    for i, sum_val in enumerate(list_sums):
        if sum_val != expected_sum:
            kb_num = i + 1
            unexpected_sums.append((kb_num, i*1000, sum_val))
    
    if unexpected_sums:
        print(f"\nANSWER: YES - Found {len(unexpected_sums)} segment(s) with sums ≠ 1000:\n")
        for kb_num, position, sum_val in unexpected_sums:
            print(f"  • Kilobase {kb_num} (position {position}): sum = {sum_val}")
            if kb_num == len(list_sums):
                print(f"    → This is the LAST segment (expected if sequence length not divisible by 1000)")
            else:
                print(f"    → UNEXPECTED: This is NOT the last segment!")
    else:
        print("\nANSWER: NO - All segments have the expected sum of 1000.")
    
    print("\n" + "="*70)
    print("QUESTION 2: General explanation for differences")
    print("="*70)
    
    print("\nANSWER:")
    print(f"\nSequence Information:")
    print(f"  • Total sequence length: {len(sequence):,} bp")
    print(f"  • Number of complete kilobases: {len(sequence) // 1000}")
    print(f"  • Remainder bases in last segment: {len(sequence) % 1000} bp")
    
    print(f"\nExpected vs Observed Results:")
    
    if len(sequence) % 1000 != 0:
        last_segment_size = len(sequence) % 1000
        print(f"  • Expected: All segments = 1000 bp, except last = {last_segment_size} bp")
        print(f"  • Observed: {len([s for s in list_sums if s == 1000])} segments with sum=1000")
        if unexpected_sums:
            print(f"             {len(unexpected_sums)} segment(s) with different sums")
    else:
        print(f"  • Expected: All segments = 1000 bp (sequence length divisible by 1000)")
        print(f"  • Observed: {len([s for s in list_sums if s == 1000])} segments with sum=1000")
    
    print(f"\nExplanation of Differences:")
    print(f"  1. INCOMPLETE FINAL SEGMENT:")
    print(f"     → The sequence has {len(sequence)} total bases")
    print(f"     → {len(sequence)} ÷ 1000 = {len(sequence) // 1000} complete kilobases + {len(sequence) % 1000} remaining bases")
    print(f"     → Therefore, the last segment contains only {len(sequence) % 1000} nucleotides")
    
    print(f"\n  2. NUCLEOTIDE COUNTING METHOD:")
    print(f"     → We count only A, C, G, T (NOT N) in the list sums")
    print(f"     → N bases are tracked separately but excluded from [A,C,G,T] lists")
    
    # Count total N's
    total_n = sum(counts['N'] for counts in kilobase_dict.values())
    if total_n > 0:
        print(f"     → Total N bases in sequence: {total_n}")
        print(f"     → This explains why some sums might be < 1000")
    
    print(f"\n  3. VERIFICATION:")
    segments_with_n = []
    for i, pos in enumerate(sorted(kilobase_dict.keys())):
        n_count = kilobase_dict[pos]['N']
        if n_count > 0:
            segments_with_n.append((i+1, pos, n_count, list_sums[i]))
    
    if segments_with_n:
        print(f"     → Segments containing N bases:")
        for kb_num, pos, n_count, sum_val in segments_with_n[:5]:  # Show first 5
            print(f"       KB {kb_num} (pos {pos}): {n_count} N's, A+C+G+T sum = {sum_val}")
        if len(segments_with_n) > 5:
            print(f"       ... and {len(segments_with_n) - 5} more segments with N bases")
    else:
        print(f"     → No N bases found in sequence")
    
    # Summary statistics
    print("\n" + "="*70)
    print("SUMMARY STATISTICS")
    print("="*70)
    
    total_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
    for counts in kilobase_dict.values():
        for nucleotide in total_counts:
            total_counts[nucleotide] += counts[nucleotide]
    
    total_bases = sum(total_counts.values())
    
    print(f"\nTotal nucleotide counts across entire sequence:")
    for nucleotide, count in total_counts.items():
        percentage = (count / total_bases) * 100 if total_bases > 0 else 0
        print(f"  {nucleotide}: {count:>8,} ({percentage:>5.2f}%)")
    print(f"  Total: {total_bases:>8,}")

if __name__ == "__main__":
    main()