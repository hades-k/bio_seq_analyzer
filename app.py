from flask import Flask, render_template, request, redirect, url_for, flash, send_file
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import io
from tools import Parser, SequenceAligner, MotifFinder
from sequence import MitochondrialDNA
from comparer import SequenceComparer, ConservedMotifAnalyzer, MultiAligner  # ✅ comparer import

app = Flask(__name__)
app.secret_key = 'bioseq2024'

# --- Helper Classes ---
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
        return {
            'count': len(self.df),
            'min_length': lengths.min(),
            'max_length': lengths.max(),
            'mean_length': lengths.mean(),
        }

    def get_gc_contents(self):
        return [obj.gc_content for obj in self.mito_objs]

    def get_names(self):
        return [obj._MitochondrialDNA__name for obj in self.mito_objs]

    def get_sequences(self):
        return self.mito_objs

class HeatmapBuilder:  # ✅ new class
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
    return render_template('index.html')

@app.route('/summary')
def summary():
    stats = fasta_manager.get_stats()
    gc_contents = fasta_manager.get_gc_contents()
    names = fasta_manager.get_names()
    return render_template('summary.html', stats=stats, gc_contents=gc_contents, names=names)

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
    return render_template('motif.html', results=results)

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
    return render_template('align.html', names=names, result=result)

@app.route('/heatmap.png')  # ✅ uses HeatmapBuilder
def heatmap():
    builder = HeatmapBuilder(fasta_manager.get_sequences())
    matrix, names = builder.build_score_matrix()

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(matrix, cmap='Blues')
    ax.set_xticks(np.arange(len(names)))
    ax.set_yticks(np.arange(len(names)))
    ax.set_xticklabels(names, rotation=90)
    ax.set_yticklabels(names)
    fig.colorbar(im, ax=ax, label='Alignment Score')
    plt.tight_layout()

    img = io.BytesIO()
    fig.savefig(img, format='png')
    img.seek(0)
    return send_file(img, mimetype='image/png')

if __name__ == '__main__':
    app.run(debug=True)
