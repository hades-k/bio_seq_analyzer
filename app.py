from flask import Flask, render_template_string, request, redirect, url_for, flash, send_file
import os
import io
import base64
import matplotlib.pyplot as plt
import numpy as np
from tools import Parser, SequenceAligner, MotifFinder
from sequence import MitochondrialDNA
from comparer import SequenceComparer, ConservedMotifAnalyzer, MultiAligner

app = Flask(__name__)
app.secret_key = 'bioseq2024'

base_html_py = """<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <title>mtDNA Analyzer</title>
  <link href=\"https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css\" rel=\"stylesheet\">
</head>
<body>
  <nav class=\"navbar navbar-expand-lg navbar-dark bg-dark mb-4\">
    <div class=\"container-fluid\">
      <a class=\"navbar-brand\" href=\"/\">mtDNA Analyzer</a>
      <div class=\"collapse navbar-collapse\">
        <ul class=\"navbar-nav\">
          <li class=\"nav-item\"><a class=\"nav-link\" href=\"/summary\">Summary</a></li>
          <li class=\"nav-item\"><a class=\"nav-link\" href=\"/motif\">Motif Search</a></li>
          <li class=\"nav-item\"><a class=\"nav-link\" href=\"/align\">Alignment</a></li>
          <li class=\"nav-item\"><a class=\"nav-link\" href=\"/heatmap_view\">Heatmap</a></li>
        </ul>
      </div>
    </div>
  </nav>
  <div class=\"container\">
    {% with messages = get_flashed_messages() %}
      {% if messages %}
        <div class=\"alert alert-info\">{{ messages[0] }}</div>
      {% endif %}
    {% endwith %}
    {{ content|safe }}
  </div>
</body>
</html>"""

index_html_py = """<h2>Upload a FASTA File</h2>
<form method=\"POST\" enctype=\"multipart/form-data\">
  <div class=\"mb-3\">
    <input class=\"form-control\" type=\"file\" name=\"file\" required>
  </div>
  <button class=\"btn btn-primary\" type=\"submit\">Upload</button>
</form>"""

summary_html_py = """<h2>Summary Statistics</h2>
<ul class=\"list-group\">
  <li class=\"list-group-item\">Total Sequences: {{ stats.count }}</li>
  <li class=\"list-group-item\">Min Length: {{ stats.min_length }}</li>
  <li class=\"list-group-item\">Max Length: {{ stats.max_length }}</li>
  <li class=\"list-group-item\">Mean Length: {{ stats.mean_length | round(2) }}</li>
</ul>
<div class=\"row mt-4\">
    <div class=\"col-md-6\">
        <h3 class=\"mt-4\">GC Content Distribution</h3>
        <img src=\"/plot.png\" alt=\"GC Content Distribution\" class=\"img-fluid\">
    </div>
    <div class=\"col-md-6\">
        <h3 class=\"mt-4\">GC Content Frequency</h3>
        <img src=\"/gc_histogram.png\" alt=\"GC Content Frequency Histogram\" class=\"img-fluid\">
    </div>
</div>
<h3 class=\"mt-4\">GC Content (%)</h3>
<table class=\"table table-striped\">
  <thead><tr><th>Species</th><th>GC</th></tr></thead>
  <tbody>
    {% for name, gc in zip(names, gc_contents) %}
      <tr><td>{{ name }}</td><td>{{ gc | round(2) }}</td></tr>
    {% endfor %}
  </tbody>
</table>"""

motif_html_py = """<h2>Motif Search</h2>
<form method=\"POST\">
  <div class=\"mb-3\">
    <input class=\"form-control\" type=\"text\" name=\"motif\" placeholder=\"Enter motif (e.g., GATC)\">
  </div>
  <button class=\"btn btn-primary\" type=\"submit\">Search Motif</button>
</form>
{% if results %}
  <div class=\"row mt-4\">
    <div class=\"col-md-12\">
        <h3 class=\"mt-4\">Motif Distribution</h3>
        <img src=\"/motif_histogram.png\" alt=\"Motif Distribution Histogram\" class=\"img-fluid\">
    </div>
  </div>
  <h3 class=\"mt-4\">Results</h3>
  <ul class=\"list-group\">
    {% for r in results %}
      <li class=\"list-group-item\">
        <strong>{{ r.name }}</strong>: {{ r.count }} hits<br>
        Positions: {{ r.positions }}
      </li>
    {% endfor %}
  </ul>
{% endif %}"""

align_html_py = """<h2>Pairwise Sequence Alignment</h2>
<form method=\"POST\">
  <div class=\"row mb-3\">
    <div class=\"col\">
      <label>Sequence 1</label>
      <select class=\"form-control\" name=\"seq1\">
        {% for i in range(names|length) %}
          <option value=\"{{ i }}\">{{ names[i] }}</option>
        {% endfor %}
      </select>
    </div>
    <div class=\"col\">
      <label>Sequence 2</label>
      <select class=\"form-control\" name=\"seq2\">
        {% for i in range(names|length) %}
          <option value=\"{{ i }}\">{{ names[i] }}</option>
        {% endfor %}
      </select>
    </div>
  </div>
  <button class=\"btn btn-primary\" type=\"submit\">Align</button>
</form>
{% if result %}
  <h3 class=\"mt-4\">Alignment Result</h3>
  <p><strong>Score:</strong> {{ result.score }}</p>
    <pre>{{ result.aligned_seq1 | join('') }}</pre>
    <pre>{{ result.matches | join('') }}</pre>
    <pre>{{ result.aligned_seq2 | join('') }}</pre>
{% endif %}"""

class FastaManager:
    def __init__(self):
        self.mito_objs = []
        self.parser = Parser('fasta')
        self.motif_results = None

    def load_sequences(self, filepath):
        self.mito_objs = self.parser.run(filepath, return_objects=True)

    def get_stats(self):
        if not self.mito_objs:
            return {'count': 0, 'min_length': 0, 'max_length': 0, 'mean_length': 0}
        lengths = [m.length for m in self.mito_objs]
        return {
            'count': len(lengths),
            'min_length': min(lengths),
            'max_length': max(lengths),
            'mean_length': sum(lengths) / len(lengths)
        }

    def get_gc_contents(self):
        return [m.gc_content for m in self.mito_objs]

    def get_names(self):
        return [m._MitochondrialDNA__name for m in self.mito_objs]

    def get_sequences(self):
        return self.mito_objs

fasta_manager = FastaManager()

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        file = request.files['file']
        if file:
            filepath = 'uploaded.fasta'
            file.save(filepath)
            fasta_manager.load_sequences(filepath)
            flash("File uploaded and parsed successfully!")
            return redirect("/summary")
    return render_template_string(base_html_py, content=index_html_py)

@app.route('/summary')
def summary():
    stats = fasta_manager.get_stats()
    gc_contents = fasta_manager.get_gc_contents()
    names = fasta_manager.get_names()
    return render_template_string(
        base_html_py,
        content=render_template_string(summary_html_py, stats=stats, gc_contents=gc_contents, names=names, zip=zip)
    )

@app.route('/motif', methods=['GET', 'POST'])
def motif():
    results = None
    if request.method == 'POST':
        motif = request.form.get('motif')
        results = []
        finder = MotifFinder()
        for m in fasta_manager.get_sequences():
            res = finder.run(m.sequence, motif=motif)
            results.append({
                'name': m._MitochondrialDNA__name,
                'count': res.get('count', 0),
                'positions': res.get('positions', [])
            })
        fasta_manager.motif_results = results
    return render_template_string(base_html_py, content=render_template_string(motif_html_py, results=results))

@app.route('/align', methods=['GET', 'POST'])
def align():
    result = None
    names = fasta_manager.get_names()
    if request.method == 'POST':
        idx1 = int(request.form.get('seq1'))
        idx2 = int(request.form.get('seq2'))
        aligner = SequenceAligner()
        m1 = fasta_manager.get_sequences()[idx1]
        m2 = fasta_manager.get_sequences()[idx2]
        aligner.run(m1.sequence, m2.sequence)
        result = aligner.get_alignment_data()
    return render_template_string(base_html_py, content=render_template_string(align_html_py, names=names, result=result))

@app.route('/plot.png')
def plot_png():
    gc_contents = fasta_manager.get_gc_contents()
    names = fasta_manager.get_names()
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(names, gc_contents)
    ax.set_ylabel('GC Content (%)')
    ax.set_title('GC Content Distribution')
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.tight_layout()
    img = io.BytesIO()
    fig.savefig(img, format='png')
    img.seek(0)
    return send_file(img, mimetype='image/png')

@app.route('/gc_histogram.png')
def gc_histogram_png():
    gc_contents = fasta_manager.get_gc_contents()
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(gc_contents, bins=10, edgecolor='black')
    ax.set_xlabel('GC Content (%)')
    ax.set_ylabel('Frequency')
    ax.set_title('GC Content Frequency Histogram')
    plt.tight_layout()
    img = io.BytesIO()
    fig.savefig(img, format='png')
    img.seek(0)
    return send_file(img, mimetype='image/png')

@app.route('/motif_histogram.png')
def motif_histogram_png():
    if not fasta_manager.motif_results:
        return "", 404
    names = [r['name'] for r in fasta_manager.motif_results]
    counts = [r['count'] for r in fasta_manager.motif_results]
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(names, counts)
    ax.set_ylabel('Motif Count')
    ax.set_title('Motif Distribution')
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.tight_layout()
    img = io.BytesIO()
    fig.savefig(img, format='png')
    img.seek(0)
    return send_file(img, mimetype='image/png')

@app.route('/heatmap_view')
def heatmap_view():
    sequences = fasta_manager.get_sequences()
    names = fasta_manager.get_names()

    if not sequences or len(sequences) < 2:
        return render_template_string(base_html_py, content="<p>Not enough sequences to generate heatmap.</p>")

    try:
        aligner = MultiAligner(sequences)
        print(f"Starting pairwise alignment for {len(sequences)} sequences...")
        score_map = aligner.pairwise_scores()
        print("Pairwise alignment complete.")

        n = len(sequences)
        matrix = np.zeros((n, n))
        for (i, j), score in score_map.items():
            matrix[i][j] = score
            matrix[j][i] = score

        fig, ax = plt.subplots(figsize=(10, 8))
        im = ax.imshow(matrix, cmap='Blues')
        ax.set_xticks(np.arange(n))
        ax.set_yticks(np.arange(n))
        ax.set_xticklabels(names, rotation=90)
        ax.set_yticklabels(names)
        fig.colorbar(im, ax=ax, label='Alignment Score')
        plt.tight_layout()

        os.makedirs("static", exist_ok=True)
        fig.savefig("static/heatmap_output.png")
        print("Heatmap saved to static/heatmap_output.png")

        return render_template_string(base_html_py, content="""
            <h2>Alignment Heatmap</h2>
            <img src="/static/heatmap_output.png" class="img-fluid" alt="Alignment Heatmap">
        """)
    except Exception as e:
        print(f"Heatmap generation failed: {e}")
        return render_template_string(base_html_py, content=f"<p>Error: {e}</p>")
      
if __name__ == '__main__':
  app.run(debug=True)
