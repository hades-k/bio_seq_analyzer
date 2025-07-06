# BioSeq Analyzer

## Project Overview
This project is a modular software solution for analyzing mitochondrial DNA (mtDNA) sequences from multiple species. It is designed using object-oriented programming (OOP) principles and provides both a command-line and a web-based interface (using Flask) for genomic data analysis.

The main goals are:
- Parse and organize FASTA files containing mtDNA sequences.
- Represent mitochondrial DNA and motifs as objects.
- Perform sequence analysis (GC content, motif search, alignments).
- Compare sequences and visualize results.
- Provide an easy-to-use web interface for all main functionalities.

---

## Project Structure

```
├── app.py                # Flask web application (Part 4)
├── tools.py              # Parsing, alignment, and motif-finding tools (Parts 1 & 2)
├── sequence.py           # MitochondrialDNA and sequence representation (Part 2)
├── comparer.py           # High-level analysis and integration (Part 3)
├── templates/            # HTML templates for the web interface
│   ├── index.html
│   ├── summary.html
│   ├── motif.html
│   ├── align.html
│   └── pairwise.html
├── synthetic_mtDNA_dataset.fasta  # Example input file
└── README.md             # This file
```

---

## How to Run the Web App

1. **Install dependencies:**
   - Make sure you have Python 3 and `pip` installed.
   - Install required packages:
     ```bash
     pip install flask pandas biopython numpy
     ```

2. **Start the Flask app:**
   ```bash
   python app.py
   ```
   The app will run on `http://127.0.0.1:5000/` by default.

3. **Open your browser:**
   - Go to `http://127.0.0.1:5000/`
   - Upload a FASTA file and explore the analysis features.

---

## Code Explanation

### Part 1: Parsing the Input (`tools.py`)
- **Parser class:** Reads FASTA files and stores the data in a Pandas DataFrame. It extracts sequence IDs, descriptions, and sequences, and calculates sequence lengths. The parser is reusable for other formats supported by Biopython.

### Part 2: Genomic Representation and Analysis
- **MitochondrialDNA class (`sequence.py`):** Represents a single mtDNA sequence. Provides methods to get the sequence, its length, GC content, extract subsequences, and find irregular bases.
- **MotifFinder class (`tools.py`):** Finds specific motifs or discovers frequent k-mers in a sequence. Can search for a motif or discover conserved motifs based on frequency.
- **SequenceAligner class (`tools.py`):** Performs global or local pairwise sequence alignment. Returns alignment score, matches, mismatches, and gaps.

### Part 3: Integration of Analytical Capabilities (`comparer.py`)
- **SequenceComparer:** Compares all pairs of sequences or compares each sequence to a reference. Uses the alignment tools to generate scores and statistics.
- **ConservedMotifAnalyzer:** Finds motifs that are conserved (present in all sequences) using the MotifFinder.
- **AlignmentVisualizer:** (for command-line use) Prints alignments in a readable format.

### Part 4: User Interface (`app.py` and `templates/`)
- **app.py:** Implements the Flask web app. Uses OOP controller classes to manage parsing, motif analysis, and alignments. Handles file uploads, user input, and routes for each feature.
- **HTML templates:**
  - `index.html`: Upload a FASTA file.
  - `summary.html`: View sequence statistics and GC content.
  - `motif.html`: Search for motifs or discover conserved motifs.
  - `align.html`: Perform and view pairwise alignments.
  - `pairwise.html`: See all pairwise alignment scores.

---

## Example Workflow
1. Upload a FASTA file on the main page.
2. View summary statistics and GC content for all sequences.
3. Search for a specific motif or discover conserved motifs.
4. Perform pairwise alignments between any two sequences.
5. View all pairwise alignment scores in a table.

---

## Notes
- The code is written to be clear and modular, following OOP principles like encapsulation, inheritance, and abstraction.
- You can extend the system by adding new analysis tools or web features easily.
- For any issues, check the terminal for error messages or use the command-line interface in the modules for debugging.

---

## Credits
- Developed for Advanced Programming 2024/2025
- Project Specification: Andrea Giovanni Nuzzolese
- Author: [Your Name Here] 