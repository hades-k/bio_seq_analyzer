# BioSeq Analyzer

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
│   ├── index.html
│   ├── motif.html
│   └── summary.html
├── uploads/              # Directory for uploaded FASTA files
├── P04637.fasta          # Example FASTA file
├── synthetic_mtDNA_dataset.fasta  # Example input file
├── requirements.txt      # Python dependencies
├── documentation.md      # Project documentation
└── README.md             # This file
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
