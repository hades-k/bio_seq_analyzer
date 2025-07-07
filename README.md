# BioSeq Analyzer

## Summary

**BioSeq Analyzer** is an interactive web application and Python toolkit for analyzing mitochondrial DNA (mtDNA) sequences. It provides intuitive ways to upload FASTA files, inspect DNA properties like GC content, identify conserved motifs, and perform global or local pairwise alignments. Built with modular object-oriented design, it is easily extendable for future genomic tools.

- Analyze mitochondrial genomes from any species  
- Visualize sequence statistics and motif frequencies  
- Align and compare DNA sequences using classic algorithms  
- Designed for extensibility, reusability, and teaching purposes


## Project Overview
BioSeq Analyzer is a modular software solution for analyzing mitochondrial DNA (mtDNA) sequences. Built with object-oriented programming (OOP) principles, it offers both command-line and web-based interfaces (using Flask) for genomic data analysis. Key functionalities include parsing FASTA files, representing DNA sequences and motifs as objects, performing sequence analysis (GC content, motif search, alignments), and comparing sequences.

---

## Project Structure

```
├── app.py                # Flask web application
├── tools.py              # Parsing, alignment, and motif-finding tools
├── sequence.py           # MitochondrialDNA and sequence representation
├── comparer.py           # High-level analysis and integration
├── templates/            # HTML templates for the web interface
│   ├── align.html
│   ├── base.html
│   ├── compare_reference.html
│   ├── index.html
│   ├── motif.html
│   └── summary.html
├── uploads/              # Directory for uploaded FASTA files
├── P04637.fasta          # Example FASTA file
├── synthetic_mtDNA_dataset.fasta  # Example input file
├── requirements.txt      # Python dependencies
├── documentation.md      # Project documentation
├── README.md             # This file
└── UML-diagram.png       # UML diagram visible in the documentation.md
```

---

## Getting Started (Installation & Usage)

To run the Bio Seq Analyzer web application locally:

### 1. Clone the Repository
```bash
git clone https://github.com/hades-k/bio-seq-analyzer.git
cd bio-seq-analyzer
```

### 2. Set Up a Virtual Environment (Optional but recommended)
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

### 3. Install Dependencies
```bash
pip install -r requirements.txt
```

If you don’t have a `requirements.txt` yet, create one with:
```
Flask
biopython
pandas
matplotlib
```

### 4. Run the Application
```bash
python app.py
```

Navigate to [http://127.0.0.1:5000](http://127.0.0.1:5000) in your browser.

### 5. Features via Web Interface
- Upload FASTA files
- View sequence statistics (length, GC content)
- Motif search and k-mer discovery
- Global & local sequence alignment
- View plots (GC content, motif distribution)

---

## Object-Oriented Design Principles

This project is structured around modern OOP principles:

### ✔ Encapsulation
Classes like `MitochondrialDNA` and `SequenceAligner` encapsulate internal data (e.g., sequences, alignment results), exposing functionality through clean public interfaces.

### ✔ Abstraction
`Sequence` and `Tool` are abstract base classes enforcing essential methods (`run`, `report`, etc.) in all subclasses, enabling modular and consistent logic.

### ✔ Inheritance
- `MitochondrialDNA` inherits from `Sequence`
- `MotifFinder`, `SequenceAligner`, and `Parser` all inherit from `Tool`

This promotes reuse and shared structure across tools.

### ✔ Polymorphism
Generic interfaces (`run()`, `report()`) allow tools like `Parser`, `MotifFinder`, and `SequenceAligner` to be used interchangeably in pipelines and the frontend.

## Example Workflow
1. Upload a FASTA file on the main page.
2. View summary statistics and GC content for all sequences.
3. Search for a specific motif or discover conserved motifs.
4. Perform pairwise alignments between any two sequences.
5. View all pairwise alignment scores in a table.

---

## Documentation

Full documentation can be found [here](https://github.com/hades-k/bio_seq_analyzer/blob/main/documentation.md)

---

## Notes
- The code is designed to be modular and follows OOP principles.
- The system can be extended with new analysis tools or web features.
- For debugging, check terminal error messages or use the command-line interface of individual modules.

---

## Credits
- Developed for Advanced Programming 2024/2025
- Professor: Andrea Giovanni Nuzzolese
- Authors: Salma Maria Mede, Aidana Kaukenova
