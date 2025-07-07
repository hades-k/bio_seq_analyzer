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
- `MotifFinder`, `SequenceAligner`, `SequenceComparer` (used by)

---

**Class:** `Tool`

**Responsibilities:**
- Abstract superclass for low-level sequence manipulation tools
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

- `FastaManager` (calls Parser)
- `MitochondrialDNA` (instantiated from parsed records)
- `Tool` (superclass)

----

**Class:** `SequenceAligner`

**Responsibilities:**
- Align two biological sequences using Needleman-Wunsch (global) or Smith-Waterman (local) algorithms
- Calculate the alignment score
- Calculate the number of matches/mismatches/gaps
- Optionally print the alignment score matrix
- Provide graphical representation of the alignment

**Collaborators:**
- `MitochondrialDNA` (input type)
- `Tool` (superclass)
- `SequenceAlignWrapper` (used internally)
- `SequenceComparer`, `MultiAligner` (use `SequenceAligner`)

----

**Class:** `MotifFinder`

**Responsibilities:**
- Search for a given motif in a DNA sequence and return its positions
- Discover overrepresented k-mers (motifs) in a sequence above a frequency threshold
- Return and report the last search result (motif matches or discovered motifs)

**Collaborators:**
- `Tool` (Superclass)
- `MitochondrialDNA` (provides input sequences)

----

**Class:** `FastaManager`

**Responsibilities:**
- Load and manage multiple FASTA datasets, each containing `MitochondrialDNA` objects
- Set the currently active dataset
- Provide statistics (GC content, length range, names, count) for the current dataset

**Collaborators:**
- `Parser` (used to load sequences)
- `MitochondrialDNA` (objects managed)

----

**Class:** `SequenceAlignWrapper`

**Responsibilities:**
- Wrap `SequenceAligner` to simplify usage and expose a cleaner interface
- Perform alignment and return alignment data.

**Collaborators:**
- `SequenceAligner` (used internally)
- `AlignmentVisualizer`, `SequenceComparer` (use this wrapper)

----

**Class:** `AlignmentVisualizer`

**Responsibilities:**
- Display pairwise alignment between sequences in a formatted, readable layout
- Show alignment symbols, score, matches/mismatches/gaps

**Collaborators:**
- `SequenceAlignWrapper` (used to perform alignment)
- `MitochondrialDNA` (input data)

----

**Class:** `SequenceComparer`

**Responsibilities:**
- Compare all sequence pairs (or each to a reference) using alignments
- Store and return summary statistics for each comparison

**Collaborators:**
- `MitochondrialDNA` (input data)
- `SequenceAlignWrapper` (used to perform alignment)

---

**Class:** `MultiAligner`

**Responsibilities:**
- Perform pairwise alignments and return a matrix of scores

**Collaborators:**
- `MitochondrialDNA` (input data)
- `SequenceAligner` (performs alignments)

---



## Method Documentation

### `Sequence` Class

| Method/Property | Input | Output | Example |
| :--- | :--- | :--- | :--- |
| `__init__(df)` | `df`: pandas DataFrame row | Initializes the abstract base class. | `seq_obj = Sequence(df_row)` |
| `sequence` (abstract property) | None | The biological sequence: `str` | `seq = seq_obj.sequence` |
| `length` (abstract property) | None | Length of the sequence: `int` | `l = seq_obj.length` |


### `MitochondrialDNA` Class

| Method/Property | Input | Output | Example |
| :--- | :--- | :--- | :--- |
| `__init__(df)` | `df`: pandas DataFrame row | A `MitochondrialDNA` object | `dna = MitochondrialDNA(df_row)` |
| `sequence` (property) | None | The DNA sequence: `str` | `seq = dna.sequence` |
| `length` (property) | None | Length of the sequence: `int` | `l = dna.length` |
| `gc_content` (property) | None | Percentage of GC content: `float` | `gc = dna.gc_content` |
| `get_subsequence(start, end)` | `start`: `int`, `end`: `int` | The subsequence: `str` | `sub = dna.get_subsequence(10, 50)` |
| `find_irregular_bases()` | None | A list of non-standard bases found: `list[str]` | `irregulars = dna.find_irregular_bases()` |
| `name` (property) | None | The name of the sequence: `str` | `name = dna.name` |


