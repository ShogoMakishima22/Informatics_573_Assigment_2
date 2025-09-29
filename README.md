# DNA Sequence Analysis Program

## Overview
This Python program performs comprehensive analysis of DNA sequences from FASTA files, specifically designed to analyze the chr1_GL383518v1_alt sequence. The program reads DNA sequences, generates reverse complements, and provides detailed nucleotide composition analysis broken down by kilobase segments. It includes automated validation to compare expected versus observed results and provides explanations for any discrepancies found.

## Author
- **Name:** Venkatesh Pramod Joshi
- **Programming Language:** Python
- **Date:** 09/28/2025

## Table of Contents
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Program Structure](#program-structure)
- [Data Structures](#data-structures)
- [Analysis Components](#analysis-components)
- [Output Description](#output-description)
- [Troubleshooting](#troubleshooting)
- [Customization](#customization)
- [Limitations](#limitations)
- [Future Enhancements](#future-enhancements)
- [References](#references)

## Features

### Core Functionality
1. **FASTA File Reading**: Robust parsing of FASTA format files with error handling
2. **Reverse Complement Generation**: Creates Watson-Crick-Franklin complementary DNA sequences (A↔T, G↔C)
3. **Kilobase Segmentation**: Divides sequences into 1000 base pair segments for detailed analysis
4. **Nucleotide Counting**: Tracks A, C, G, T, and N (ambiguous) nucleotides separately
5. **Data Validation**: Automatically identifies segments with unexpected nucleotide counts
6. **Statistical Summary**: Provides total counts and percentage composition across the entire sequence

### Analysis Outputs
- Specific nucleotide position queries (10th, 758th, 79th letters)
- Reverse complement substring extraction (positions 500-800)
- Formatted table of nucleotide counts per kilobase
- Nested dictionary structure for programmatic access
- List-based data structures for statistical analysis
- Comprehensive validation reporting with explanations

## Requirements

### System Requirements
- **Python Version**: Python 3.6 or higher
- **Operating System**: Cross-platform (Windows, macOS, Linux)
- **Memory**: Sufficient RAM to load entire sequence into memory (typically < 1 GB for most sequences)
- **Dependencies**: None (uses only Python standard library)

### Input File Requirements
- **Format**: FASTA format (.fa, .fasta)
- **Filename**: `chr1_GL383518v1_alt.fa` (default, can be modified)
- **Location**: Must be in the same directory as the Python script
- **Content**: Standard DNA nucleotides (A, C, G, T, N)
- **Structure**: Single sequence per file with header line starting with '>'

## Installation

### Step 1: Download the Script
Save the Python script as `data_analysis_assgn2.py` in your working directory.

### Step 2: Prepare Input File
1. Obtain your FASTA file (e.g., `chr1_GL383518v1_alt.fa`)
2. Place it in the same directory as `data_analysis_assgn2.py`
3. Verify the file is properly formatted (header line with '>' followed by sequence lines)

### Step 3: Verify Python Installation
```bash
python --version
# or
python3 --version
```
Ensure version is 3.6 or higher.

### No Additional Dependencies Required
This program uses only Python's standard library, so no pip installations are necessary.

## Usage

### Basic Execution
```bash
python data_analysis_assgn2.py
```

Or on some systems:
```bash
python3 data_analysis_assgn2.py
```

### Command Line Options
Currently, the program does not accept command-line arguments. To analyze different files, modify the `filename` variable in the `main()` function.

## Program Structure

### Function Documentation

#### `read_fasta_sequence(filename)`
**Purpose**: Reads and parses a FASTA format file.

**Parameters**:
- `filename` (str): Name of the FASTA file to read

**Returns**:
- `str`: DNA sequence with all letters in uppercase
- `None`: If file not found or error occurs

**Features**:
- Skips header lines (lines starting with '>')
- Converts all nucleotides to uppercase for consistency
- Comprehensive error handling

**Error Handling**:
- `FileNotFoundError`: Displays specific error message
- General exceptions: Catches and reports any other reading errors

```python
sequence = read_fasta_sequence('chr1_GL383518v1_alt.fa')
```

---

#### `create_reverse_complement(sequence)`
**Purpose**: Generates the reverse complement of a DNA sequence using Watson-Crick-Franklin base pairing rules.

**Parameters**:
- `sequence` (str): DNA sequence string

**Returns**:
- `str`: Reverse complement sequence

**Base Pairing Rules**:
- A (Adenine) ↔ T (Thymine)
- G (Guanine) ↔ C (Cytosine)
- N (Ambiguous) ↔ N (Ambiguous)

**Algorithm**:
1. Create complement by replacing each base with its pair
2. Reverse the complemented sequence
3. Return the reverse complement

```python
reverse_comp = create_reverse_complement(sequence)
```

---

#### `analyze_sequence_by_kilobase(sequence)`
**Purpose**: Creates a nested dictionary containing nucleotide counts for each 1000 bp segment.

**Parameters**:
- `sequence` (str): DNA sequence to analyze

**Returns**:
- `dict`: Nested dictionary structure
  - Keys: Starting positions (0, 1000, 2000, ...)
  - Values: Dictionaries with nucleotide counts

**Dictionary Structure**:
```python
{
    0: {'A': 245, 'C': 255, 'G': 250, 'T': 250, 'N': 0},
    1000: {'A': 260, 'C': 240, 'G': 245, 'T': 255, 'N': 0},
    ...
}
```

**Features**:
- Processes sequence in 1000 bp chunks
- Handles incomplete final segment automatically
- Counts all five nucleotide types (A, C, G, T, N)

**Usage Example**:
```python
kilobase_dict = analyze_sequence_by_kilobase(sequence)
a_count_in_kb6 = kilobase_dict[5000]['A']  # A's in positions 5000-5999
```

---

#### `print_kilobase_analysis(analysis_dict, sequence_length)`
**Purpose**: Displays formatted table of nucleotide composition per kilobase.

**Parameters**:
- `analysis_dict` (dict): Dictionary from `analyze_sequence_by_kilobase()`
- `sequence_length` (int): Total sequence length in base pairs

**Output Format**:
```
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
```

**Features**:
- Formatted columns with proper alignment
- Comma-separated numbers for readability
- Total nucleotide count per segment
- Sorted by position for sequential display

---

#### `create_nucleotide_lists(analysis_dict)`
**Purpose**: Converts nested dictionary to list of lists for statistical analysis.

**Parameters**:
- `analysis_dict` (dict): Dictionary from `analyze_sequence_by_kilobase()`

**Returns**:
- `list`: List of lists, each containing [A, C, G, T] counts (excludes N)

**Output Structure**:
```python
[
    [245, 255, 250, 250],  # Kilobase 1
    [260, 240, 245, 255],  # Kilobase 2
    [248, 252, 249, 251],  # Kilobase 3
    ...
]
```

**Note**: N bases are intentionally excluded from these lists.

---

#### `main()`
**Purpose**: Orchestrates all analysis functions and generates comprehensive report.

**Process Flow**:
1. Read FASTA sequence file
2. Extract and display specific nucleotide positions
3. Generate reverse complement and extract substrings
4. Perform kilobase analysis
5. Create nucleotide count lists
6. Calculate and validate sums
7. Answer analysis questions
8. Display summary statistics

**Returns**: None (prints all output to console)

## Data Structures

### Nested Dictionary (kilobase_dict)
**Structure**:
```python
{
    position (int): {
        'A': count (int),
        'C': count (int),
        'G': count (int),
        'T': count (int),
        'N': count (int)
    }
}
```

**Example**:
```python
{
    0: {'A': 245, 'C': 255, 'G': 250, 'T': 250, 'N': 0},
    1000: {'A': 260, 'C': 240, 'G': 245, 'T': 255, 'N': 0},
    2000: {'A': 248, 'C': 252, 'G': 249, 'T': 251, 'N': 0}
}
```

**Access Pattern**:
```python
# Get all counts for kilobase starting at position 5000
kb6_counts = kilobase_dict[5000]

# Get specific nucleotide count
a_count = kilobase_dict[5000]['A']
```

### Nucleotide Lists
**Structure**: List of lists `[[A, C, G, T], [A, C, G, T], ...]`

**Example**:
```python
[
    [245, 255, 250, 250],  # Kilobase 1: positions 0-999
    [260, 240, 245, 255],  # Kilobase 2: positions 1000-1999
    [248, 252, 249, 251]   # Kilobase 3: positions 2000-2999
]
```

**Purpose**: Simplified format for statistical calculations and list comprehensions.

## Analysis Components

### Part 1: Reading DNA Sequence
**Outputs**:
- Total sequence length in base pairs
- 10th nucleotide (position index 9)
- 758th nucleotide (position index 757)

**Example Output**:
```
=== Part 1: Reading DNA Sequence ===
Sequence length: 182896 bp
10th letter: A
758th letter: G
```

### Part 2: Reverse Complement Analysis
**Outputs**:
- Reverse complement sequence length
- 79th letter of reverse complement (position index 78)
- Substring from position 500 to 800 (301 nucleotides)

**Example Output**:
```
=== Part 2: Reverse Complement Analysis ===
Reverse complement length: 182896 bp
79th letter of reverse complement: T
Letters 500-800 of reverse complement: ATCGATCG...
```

### Part 3: Kilobase Analysis
**Outputs**:
- Number of kilobase segments
- Detailed table with nucleotide counts per segment
- Dictionary access examples

**Features**:
- Complete breakdown of all segments
- Visual table for quick reference
- Programmatic access examples

### Part 4: Nucleotide List Analysis
**Outputs**:
- First kilobase nucleotide counts [A, C, G, T]
- Total number of lists created
- Sum of each list (first 5 shown)
- Complete list of all sums

**Example Output**:
```
=== Part 4: Nucleotide List Analysis ===
First 1000 bp nucleotide counts [A, C, G, T]: [245, 255, 250, 250]
Created 183 nucleotide count lists
Kilobase 1 sum: 1000
Kilobase 2 sum: 1000
...
All list sums: [1000, 1000, 1000, ..., 896]
```

### Part 5: Analysis Questions

#### Question 1: Sums Not Equal to 1000
**Analysis Performed**:
- Identifies all segments with sums ≠ 1000
- Reports position and actual sum for each
- Distinguishes expected (last segment) from unexpected deviations

**Example Output**:
```
======================================================================
QUESTION 1: Are there any lists whose sums are not equal to 1000?
======================================================================

ANSWER: YES - Found 1 segment(s) with sums ≠ 1000:

  • Kilobase 183 (position 182000): sum = 896
    → This is the LAST segment (expected if sequence length not divisible by 1000)
```

#### Question 2: Explanation of Differences
**Analysis Provided**:
1. Sequence information (total length, complete kilobases, remainder)
2. Expected vs observed comparison
3. Detailed explanation of differences
4. N base analysis if present
5. Verification of segments containing N bases

**Example Output**:
```
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

Explanation of Differences:
  1. INCOMPLETE FINAL SEGMENT:
     → The sequence has 182896 total bases
     → 182896 ÷ 1000 = 182 complete kilobases + 896 remaining bases
     → Therefore, the last segment contains only 896 nucleotides

  2. NUCLEOTIDE COUNTING METHOD:
     → We count only A, C, G, T (NOT N) in the list sums
     → N bases are tracked separately but excluded from [A,C,G,T] lists

  3. VERIFICATION:
     → No N bases found in sequence
```

### Summary Statistics
**Outputs**:
- Total count for each nucleotide type (A, C, G, T, N)
- Percentage composition
- Total base pairs analyzed

**Example Output**:
```
======================================================================
SUMMARY STATISTICS
======================================================================

Total nucleotide counts across entire sequence:
  A:   45,724 (25.00%)
  C:   45,724 (25.00%)
  G:   45,724 (25.00%)
  T:   45,724 (25.00%)
  N:        0 ( 0.00%)
  Total: 182,896
```

## Output Description

### Complete Output Structure
The program generates output in five main sections, each clearly labeled with headers:

1. **Part 1**: Basic sequence information and specific position queries
2. **Part 2**: Reverse complement analysis with substring extraction
3. **Part 3**: Comprehensive kilobase-by-kilobase breakdown table
4. **Part 4**: Nucleotide list generation and sum calculations
5. **Part 5**: Validation questions with detailed answers
6. **Summary**: Overall statistics for the entire sequence

### Output Formatting
- Clear section headers with separator lines (=== or ===)
- Aligned columns in tables
- Comma-separated large numbers for readability
- Bullet points for lists
- Indentation for hierarchical information

## Troubleshooting

### Common Issues and Solutions

#### Issue 1: File Not Found Error
**Error Message**: `Error: File chr1_GL383518v1_alt.fa not found.`

**Solutions**:
1. Verify the file is in the same directory as the script
2. Check the filename spelling (case-sensitive on Unix/Linux)
3. Ensure the file extension is correct (.fa, .fasta)
4. Use absolute path if file is in different directory

**How to Fix**:
```python
# Modify this line in main() function
filename = '/full/path/to/chr1_GL383518v1_alt.fa'
```

---

#### Issue 2: Empty or No Output
**Possible Causes**:
- File is empty
- File contains only header with no sequence
- Sequence contains only unsupported characters

**Solutions**:
1. Open the FASTA file and verify it contains sequence data
2. Check for proper FASTA format (header line with '>', followed by sequence)
3. Ensure sequence contains A, C, G, T, or N characters

---

#### Issue 3: Unexpected Segment Sums in Middle of Sequence
**Symptom**: Segments other than the last one have sums < 1000

**Explanation**: This occurs when segments contain N (ambiguous) nucleotides

**Verification**:
Look for this in the output:
```
3. VERIFICATION:
   → Segments containing N bases:
     KB X (pos YYYY): Z N's, A+C+G+T sum = WWW
```

**Note**: This is not an error; N bases are correctly tracked separately

---

#### Issue 4: Memory Error with Large Files
**Error**: `MemoryError`

**Cause**: Sequence file is too large to load into memory

**Solutions**:
1. Process the file in chunks (requires code modification)
2. Use a machine with more RAM
3. Use streaming approach instead of loading entire sequence

---

#### Issue 5: Incorrect Reverse Complement
**Symptom**: Unexpected characters in reverse complement

**Cause**: Sequence contains non-standard nucleotides not in complement_map

**Solution**: Verify input sequence contains only A, C, G, T, N
- Unknown characters will be converted to 'N'

## Customization

### Modifying Input File Name
**Location**: Line in `main()` function

**Original**:
```python
filename = 'chr1_GL383518v1_alt.fa'
```

**Custom**:
```python
filename = 'my_sequence.fa'  # Your file name
filename = '/path/to/my_sequence.fa'  # Full path
```

---

### Changing Segment Size
To analyze by different segment sizes (e.g., 2000 bp instead of 1000 bp):

**Location**: `analyze_sequence_by_kilobase()` function

**Original**:
```python
for start_pos in range(0, sequence_length, 1000):
    end_pos = min(start_pos + 1000, sequence_length)
```

**Modified for 2kb segments**:
```python
for start_pos in range(0, sequence_length, 2000):
    end_pos = min(start_pos + 2000, sequence_length)
```

**Note**: Also update function name and variable names for clarity.

---

### Adding Additional Nucleotide Types
To support extended nucleotide codes (R, Y, K, M, S, W, B, D, H, V):

**Location**: `create_reverse_complement()` function

**Add to complement_map**:
```python
complement_map = {
    'A': 'T', 'T': 'A',
    'G': 'C', 'C': 'G',
    'N': 'N',
    'R': 'Y',  # Purine (A or G) ↔ Pyrimidine (C or T)
    'Y': 'R',  # Pyrimidine ↔ Purine
    'K': 'M',  # Keto (G or T) ↔ Amino (A or C)
    'M': 'K',  # Amino ↔ Keto
    'S': 'S',  # Strong (G or C) ↔ Strong
    'W': 'W',  # Weak (A or T) ↔ Weak
    'B': 'V',  # Not A ↔ Not T
    'D': 'H',  # Not C ↔ Not G
    'H': 'D',  # Not G ↔ Not C
    'V': 'B'   # Not T ↔ Not A
}
```

Also update `analyze_sequence_by_kilobase()` to count these nucleotides.

---

### Exporting Results to File
Add this at the end of `main()` function:

```python
# Export to text file
with open('analysis_results.txt', 'w') as f:
    f.write(f"Sequence Length: {len(sequence)}\n")
    f.write(f"Total A: {total_counts['A']}\n")
    # Add other data as needed
```

---

### Changing Output Format
To modify table formatting, edit `print_kilobase_analysis()`:

```python
# Change column widths
print(f"{'Position':<15} {'A':<10} {'C':<10}")  # Wider columns

# Change separator style
print(f"{'*'*80}")  # Use asterisks instead of dashes
```

## Limitations

### Current Limitations
1. **Single Sequence Only**: Assumes one sequence per FASTA file; multi-sequence files will concatenate all sequences
2. **Memory Constraints**: Entire sequence must fit in memory; not suitable for very large genomes (> several GB)
3. **No Parallel Processing**: Processes sequence sequentially; could be optimized for multi-core systems
4. **Fixed Output Format**: Results printed to console only; no built-in file export
5. **N Base Exclusion**: N bases excluded from list sums, which may be unexpected for some analyses
6. **No Quality Scores**: Does not process or report FASTQ quality scores (FASTA format only)
7. **No Annotation Support**: Does not process or display sequence annotations or features
8. **ASCII Text Only**: Does not support binary sequence formats

### Known Issues
- Multi-line FASTA headers with embedded newlines may cause parsing issues
- Extremely long individual lines in FASTA files may be slow to process
- No validation of nucleotide characters (invalid characters become 'N')

