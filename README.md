# DNA Sequence Analysis Program

## Overview
This Python program performs comprehensive analysis of DNA sequences from FASTA files. It analyzes nucleotide composition, creates reverse complements, and provides detailed kilobase-by-kilobase breakdowns of sequence characteristics.

1. **Name**: Venkatesh Pramod Joshi
2. **Programming Language**: Python
3. **Date**: 09/28/2025


## Features
The program provides five main analysis components:

1. **DNA Sequence Reading**: Reads and parses FASTA format files
2. **Reverse Complement Generation**: Creates Watson-Crick-Franklin complementary sequences
3. **Kilobase Analysis**: Breaks down nucleotide composition by 1000 bp segments
4. **Nucleotide List Generation**: Creates structured data arrays for statistical analysis
5. **Data Validation**: Compares expected vs observed results with detailed explanations

## Requirements

### System Requirements
- Python 3.6 or higher
- No external libraries required (uses only Python standard library)

### Input File Requirements
- FASTA format file named `chr1_GL383518v1_alt.fa`
- File must be in the same directory as the script
- Standard DNA nucleotides (A, C, G, T, N)

## Installation

1. Download the Python script to your working directory
2. Ensure your FASTA file is named `chr1_GL383518v1_alt.fa` and placed in the same directory
3. No additional installation steps required

## Usage

### Basic Usage
```bash
python data_analysis_assgn2.py
```

### Expected Output
The program generates five sections of output:

#### Part 1: Reading DNA Sequence
- Total sequence length
- 10th nucleotide position
- 758th nucleotide position

#### Part 2: Reverse Complement Analysis
- Reverse complement sequence length
- 79th letter of reverse complement
- Letters 500-800 of reverse complement

#### Part 3: Kilobase Analysis
- Detailed table showing nucleotide counts (A, C, G, T, N) for each kilobase segment
- Example usage of the nested dictionary structure
- Format: `kilobase_dict[position]['nucleotide']`

#### Part 4: Nucleotide List Analysis
- First 1000 bp nucleotide counts [A, C, G, T]
- Complete list of all kilobase segment sums
- First 5 kilobase sums displayed

#### Part 5: Analysis Questions
- **Question 1**: Lists whose sums ≠ 1000 with specific positions and values
- **Question 2**: Detailed explanation of observed vs expected results

### Sample Output Format
```
=== Part 1: Reading DNA Sequence ===
Sequence length: 182896 bp
10th letter: A
758th letter: G

=== Part 2: Reverse Complement Analysis ===
Reverse complement length: 182896 bp
79th letter of reverse complement: T
Letters 500-800 of reverse complement: ATCG...

=== Part 3: Kilobase Analysis ===
======================================================================
NUCLEOTIDE ANALYSIS BY KILOBASE
Total sequence length: 182,896 bp
Number of kilobase segments: 183
======================================================================

Position     A        C        G        T        N        Total   
----------------------------------------------------------------------
KB 1    (0     ): 245      255      250      250      0        1000    
KB 2    (1000  ): 260      240      245      255      0        1000    
...

=== Part 4: Nucleotide List Analysis ===
First 1000 bp nucleotide counts [A, C, G, T]: [245, 255, 250, 250]
Created 183 nucleotide count lists
Kilobase 1 sum: 1000
...

=== Part 5: Analysis Questions ===
======================================================================
QUESTION 1: Are there any lists whose sums are not equal to 1000?
======================================================================

ANSWER: YES - Found 1 segment(s) with sums ≠ 1000:

  • Kilobase 183 (position 182000): sum = 896
    → This is the LAST segment (expected if sequence length not divisible by 1000)

======================================================================
QUESTION 2: General explanation for differences
======================================================================

ANSWER:

Sequence Information:
  • Total sequence length: 182,896 bp
  • Number of complete kilobases: 182
  • Remainder bases in last segment: 896 bp

Expected vs Observed Results:
  • Expected: All segments = 1000 bp, except last = 896 bp
  • Observed: 182 segments with sum=1000
             1 segment(s) with different sums

...
```

## Program Structure

### Functions

#### `read_fasta_sequence(filename)`
- **Purpose**: Reads FASTA file and extracts DNA sequence
- **Input**: Filename (string)
- **Output**: DNA sequence (string) or None if error
- **Error Handling**: Catches FileNotFoundError and general exceptions

#### `create_reverse_complement(sequence)`
- **Purpose**: Generates reverse complement using Watson-Crick-Franklin base pairing
- **Input**: DNA sequence (string)
- **Output**: Reverse complement sequence (string)

#### `analyze_sequence_by_kilobase(sequence)`
- **Purpose**: Creates nested dictionary of nucleotide counts per kilobase
- **Input**: DNA sequence (string)
- **Output**: Dictionary with structure `{position: {'A': count, 'C': count, ...}}`
- **Key Feature**: Processes sequence in 1000 bp chunks

#### `print_kilobase_analysis(analysis_dict, sequence_length)`
- **Purpose**: Displays formatted table of kilobase analysis
- **Input**: Analysis dictionary and sequence length
- **Output**: Formatted table printed to console

#### `main()`
- **Purpose**: Orchestrates all analysis functions
- **Output**: Complete analysis report

## Data Structures

### Nested Dictionary (kilobase_dict)
```python
{
    0: {'A': 245, 'C': 255, 'G': 250, 'T': 250, 'N': 0},
    1000: {'A': 260, 'C': 240, 'G': 245, 'T': 255, 'N': 0},
    2000: {'A': 248, 'C': 252, 'G': 249, 'T': 251, 'N': 0},
    ...
}
```

### Nucleotide Lists
```python
[
    [245, 255, 250, 250],  # Kilobase 1: [A, C, G, T]
    [260, 240, 245, 255],  # Kilobase 2: [A, C, G, T]
    [248, 252, 249, 251],  # Kilobase 3: [A, C, G, T]
    ...
]
```

## Key Findings & Analysis

### Expected Results
- Each kilobase segment should contain exactly 1000 nucleotides
- Exception: The final segment may contain fewer than 1000 bp if the total sequence length is not divisible by 1000

### Common Observations
1. **Incomplete Last Segment**: If sequence length is not a multiple of 1000, the last kilobase will have a sum < 1000
2. **Ambiguous Nucleotides (N)**: N bases are counted separately and excluded from [A, C, G, T] list sums
3. **Verification**: Total of all segment sums should equal total sequence length (when N bases are excluded)

### Troubleshooting

**Issue**: "File not found" error
- **Solution**: Ensure `chr1_GL383518v1_alt.fa` is in the same directory as the script

**Issue**: Unexpected segment sums in the middle of the sequence
- **Solution**: Check for N bases in those segments; they are counted separately

**Issue**: Program runs but produces no output
- **Solution**: Verify FASTA file format (should have header line starting with '>')

## Modifications & Customization

### Changing the Input File
Modify line in `main()`:
```python
filename = 'your_sequence_file.fa'  # Change this line
```

### Adjusting Segment Size
Change the step value in `analyze_sequence_by_kilobase()`:
```python
for start_pos in range(0, sequence_length, 2000):  # Changed from 1000 to 2000
    end_pos = min(start_pos + 2000, sequence_length)  # Changed from 1000 to 2000
```

### Adding Additional Nucleotides
Update the complement map in `create_reverse_complement()`:
```python
complement_map = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G',
    'N': 'N',
    'R': 'Y',  # Add purine/pyrimidine mappings
    'Y': 'R',
    ...
}
```

## Limitations
- Assumes single sequence per FASTA file
- Requires entire sequence to fit in memory
- Does not handle multi-line FASTA headers with embedded sequences
- N bases are not included in nucleotide list sums


