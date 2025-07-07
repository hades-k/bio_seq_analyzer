#app.py without sliders qaaa
from flask import Flask, render_template, request, redirect, url_for, flash, send_file
from werkzeug.utils import secure_filename
import os
import pandas as pd
from tools import Parser, SequenceAligner, MotifFinder
from sequence import MitochondrialDNA
from comparer import SequenceComparer, ConservedMotifAnalyzer

UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = {'fasta', 'fa', 'fna'}

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = 'bioseq2024'

if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)

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
        return list(self.df['name'])

    def get_sequences(self):
        return self.mito_objs

class MotifWebTool:
    def __init__(self, mito_objs):
        self.mito_objs = mito_objs
        self.finder = MotifFinder()

    def search_motif(self, motif):
        results = []
        for i, obj in enumerate(self.mito_objs):
            res = self.finder.run(obj.sequence, motif=motif)
            results.append({'index': i, 'name': obj._MitochondrialDNA__name, 'count': res['count'], 'positions': res['positions']})
        return results

    def discover_motifs(self, k=5, threshold=2):
        analyzer = ConservedMotifAnalyzer(self.mito_objs)
        return analyzer.find_conserved(k=k, threshold=threshold)

class AlignmentWebTool:
    def __init__(self, mito_objs):
        self.mito_objs = mito_objs
        self.comparer = SequenceComparer(mito_objs)

    def align_pair(self, idx1, idx2, method='global'):
        return self.comparer.compare_pair(idx1, idx2, method)

    def all_pairwise(self):
        return self.comparer.compare_all()

fasta_manager = FastaManager()

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)
            fasta_manager.parse(filepath)
            flash('File uploaded and parsed successfully!')
            return redirect(url_for('summary'))
        else:
            flash('Invalid file type!')
            return redirect(request.url)
    return render_template('index.html')

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/summary')
def summary():
    stats = fasta_manager.get_stats()
    gc_contents = fasta_manager.get_gc_contents()
    names = fasta_manager.get_names()
    return render_template('summary.html', stats=stats, names=names, gc_contents=gc_contents)

@app.route('/motif', methods=['GET', 'POST'])
def motif():
    results = None
    discovered = None
    if request.method == 'POST':
        motif = request.form.get('motif')
        k = int(request.form.get('k', 5))
        threshold = int(request.form.get('threshold', 2))
        tool = MotifWebTool(fasta_manager.get_sequences())
        if motif:
            results = tool.search_motif(motif)
        else:
            discovered = tool.discover_motifs(k=k, threshold=threshold)
    return render_template('motif.html', results=results, discovered=discovered)

@app.route('/align', methods=['GET', 'POST'])
def align():
    result = None
    names = fasta_manager.get_names()
    if request.method == 'POST':
        idx1 = int(request.form.get('seq1'))
        idx2 = int(request.form.get('seq2'))
        method = request.form.get('method', 'global')
        tool = AlignmentWebTool(fasta_manager.get_sequences())
        result = tool.align_pair(idx1, idx2, method)
    return render_template('align.html', result=result, names=names)

@app.route('/pairwise')
def pairwise():
    tool = AlignmentWebTool(fasta_manager.get_sequences())
    results = tool.all_pairwise()
    names = fasta_manager.get_names()
    return render_template('pairwise.html', results=results, names=names)

if __name__ == '__main__':
    app.run(debug=True)
