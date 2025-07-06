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
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>mtDNA Analyzer</title>
  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
</head>
<body>
  <nav class="navbar navbar-expand-lg navbar-dark bg-dark mb-4">
    <div class="container-fluid">
      <a class="navbar-brand" href="/">mtDNA Analyzer</a>
      <div class="collapse navbar-collapse">
        <ul class="navbar-nav">
          <li class="nav-item"><a class="nav-link" href="/summary">Summary</a></li>
          <li class="nav-item"><a class="nav-link" href="/motif">Motif Search</a></li>
          <li class="nav-item"><a class="nav-link" href="/align">Alignment</a></li>
        </ul>
      </div>
    </div>
  </nav>
  <div class="container">
    {% with messages = get_flashed_messages() %}
      {% if messages %}
        <div class="alert alert-info">{{ messages[0] }}</div>
      {% endif %}
    {% endwith %}
    {{ content|safe }}
  </div>
</body>
</html>"""

index_html_py = """<h2>Upload a FASTA File</h2>
<form method="POST" enctype="multipart/form-data">
  <div class="mb-3">
    <input class="form-control" type="file" name="file" required>
  </div>
  <button class="btn btn-primary" type="submit">Upload</button>
</form>"""

class FastaManager:
    def __init__(self):
        self.df = None
        self.records = None
        self.mito_objs = []

    def parse(self, filepath):
        parser = Parser('fasta')
        self.df = parser.run(filepath)
        self.records = self.df
        self.mito_objs = [MitochondrialDNA(self.df.loc[i]) for i in range(len(self.df))]
        return self.df

    def get_stats(self):
        if self.df is None:
            return {}
        lengths = self.df['length']
        stats = {
            'count': len(self.df),
            'min_length': lengths.min(),
            'max_length': lengths.max(),
            'mean_length': lengths.mean(),
        }
        return stats

    def get_gc_contents(self):
        if not self.mito_objs:
            return []
        return [obj.gc_content for obj in self.mito_objs]

    def get_names(self):
        if self.df is None:
            return []
        return [obj._MitochondrialDNA__name for obj in self.mito_objs]

    def get_sequences(self):
        return self.mito_objs

class HeatmapBuilder:
    def __init__(self, mito_objs):
        self.mito_objs = mito_objs
        self.names = [obj._MitochondrialDNA__name for obj in mito_objs]
        self.aligner = MultiAligner(mito_objs)

    def build_score_matrix(self):
        score_map = self.aligner.pairwise_scores()
        n = len(self.mito_objs)
        matrix = np.zeros((n, n))
        for (i, j), score in score_map.items():
            matrix[i][j] = score
            matrix[j][i] = score
        return matrix, self.names

fasta_manager = FastaManager()

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        file = request.files['file']
        if file:
            filepath = 'uploaded.fasta'
            file.save(filepath)
            fasta_manager.parse(filepath)
            flash("File uploaded and parsed successfully!")
            return redirect(url_for('summary'))
    return render_template_string(base_html_py, content=index_html_py)

@app.route('/heatmap.png')
def heatmap():
    builder = HeatmapBuilder(fasta_manager.get_sequences())
    matrix, names = builder.build_score_matrix()

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(matrix, cmap='Blues')

    ax.set_xticks(np.arange(len(names)))
    ax.set_yticks(np.arange(len(names)))
    ax.set_xticklabels(names, rotation=90)
    ax.set_yticklabels(names)

    fig.colorbar(im, ax=ax, label="Alignment Score")
    plt.tight_layout()

    img = io.BytesIO()
    fig.savefig(img, format='png')
    img.seek(0)
    return send_file(img, mimetype='image/png')

if __name__ == '__main__':
    app.run(debug=True)
