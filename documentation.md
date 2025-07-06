# Project Documentation

This document provides a detailed overview of the classes and methods used in the Bio Seq Analyzer project.

## CRC Cards (Class, Responsibility, Collaborators)

---

**Class:** `Sequence`

**Responsibilities:** 
- Serve as abstract base class for biological sequences
- Enforce implementation of sequence and length properties
- Store raw sequence metadata (`__info`) from a DataFrame row

**Collaborators:**
- `MitochondrialDNA` (subclass)

---

**Class:** `MitochondrialDNA`

**Responsibilities:**
- Represent a single mitochondrial DNA sequence
- Implements abstract interface from `Sequence`
- Store sequence information (ID, name, description, sequence data)
- Calculate GC content
- Extract subsequences
- Identify irregular bases

**Collaborators:**
- `Sequence` (superclass)
- `Parser` (instantiated from data parsed by)
- `FastaManager` (manages collections of)

---

**Class:** `Tool`

**Responsibilities:**
- Abstract superclass for simple sequence manipulation tools
- Define abstract methods (run, report) to be present in all subclasses

**Collaborators:**
- `Parser` (subclass)
- `SequenceAligner` (subclass)
- `MotifFinder` (subclass)

----

**Class:** `Parser`

**Responsibilities:**
- Parse a SeqIO supported file, default FASTA
- Convert parsed data into pandas DataFrame
- Optionally return MitochondrialDNA objects
- Export data to CSV
- Generate sequence summary reports

**Collaborators:**

- `FastaManager` (delegates parsing to Parser)
- `MitochondrialDNA` (instantiated from parsed records)
- `Tool` (superclass)

----

**Class:** `SequenceAligner`

**Responsibilities:**
- Align two biological sequences using Needleman-Wunsch (global) or Smith-Waterman (local) algorithms 
- Calculate the alignment score
- Calculate the number of matches/mismatches
- Optionally print the alignment score matrix
- Provide graphical representation of the alignment

**Collaborators:**
- `MitochondrialDNA` (input type)
- `Tool` (superclass)
- `SequenceAlignWrapper` (delegates alignment tasks)

----

**Class:** `MotifFinder`

**Responsibilities:**
- Search for a given motif in a DNA sequence and return its positions
- Discover overrepresented k-mers (motifs) in a sequence above a frequency threshold
- Return and report the last search result (motif matches or discovered motifs)

**Collaborators:**
- `Tool` (Superclass)
- `MitochondrialDNA` (provides input sequences)
- `ConservedMotifAnalyzer` (uses MotifFinder for motif discovery across multiple sequences)


----

**Class:**

**Responsibilities:**

**Collaborators:**

---

## Method Documentation

### `MitochondrialDNA` Class

| Method | Input | Output | Example |
| :--- | :--- | :--- | :--- |
| `__init__(df)` | `df`: pandas DataFrame row | A `MitochondrialDNA` object | `dna = MitochondrialDNA(df_row)` |
| `sequence()` | None | The DNA sequence: `str` | `seq = dna.sequence()` |
| `length()` | None | Length of the sequence: `int` | `l = dna.length()` |
| `gc_content()` | None | Percentage of GC content: `float` | `gc = dna.gc_content()` |
| `get_subsequence(start, end)` | `start`: `int`, `end`: `int` | The subsequence: `str` | `sub = dna.get_subsequence(10, 50)` |
| `find_irregular_bases()` | None | A list of non-standard bases found: `list[str]` | `irregulars = dna.find_irregular_bases()` |

