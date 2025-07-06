# Project Documentation

This document provides a detailed overview of the classes and methods used in the Bio Seq Analyzer project.

## CRC Cards (Class, Responsibility, Collaborators)

---

**Class:** `MitochondrialDNA`

**Responsibilities:**
- Represent a single mitochondrial DNA sequence.
- Store sequence information (ID, name, description, sequence data).
- Calculate sequence length.
- Calculate GC content.
- Extract subsequences.
- Identify irregular bases.

**Collaborators:**
- `Sequence` (inherits from)
- `Parser` (instantiated from data parsed by)
- `FastaManager` (manages collections of)

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

