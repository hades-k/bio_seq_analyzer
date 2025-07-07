from flask import Flask, render_template, request, redirect, url_for, flash, send_file
from werkzeug.utils import secure_filename
import os
import io
import matplotlib.pyplot as plt
from tools import Parser, MotifFinder
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
        self.motif_results = None

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
    return render_template('summary.html', stats=stats, names=names, gc_contents=gc_contents, zip=zip)

@app.route('/motif', methods=['GET', 'POST'])
def motif():
    results = None
    discovered = None
    if request.method == 'POST':
        motif_str = request.form.get('motif')
        k = int(request.form.get('k', 5))
        threshold = int(request.form.get('threshold', 2))
        sequences = fasta_manager.get_sequences()
        if motif_str:
            finder = MotifFinder()
            results = []
            for i, obj in enumerate(sequences):
                res = finder.run(obj.sequence, motif=motif_str)
                results.append({'index': i, 'name': obj.name, 'count': res['count'], 'positions': res['positions']})
            fasta_manager.motif_results = results
        else:
            analyzer = ConservedMotifAnalyzer(sequences)
            discovered, _ = analyzer.find_conserved(k=k, threshold=threshold)
    return render_template('motif.html', results=results, discovered=discovered)

@app.route('/align', methods=['GET', 'POST'])
def align():
    result = None
    names = fasta_manager.get_names()
    if request.method == 'POST':
        idx1 = int(request.form.get('seq1'))
        idx2 = int(request.form.get('seq2'))
        method = request.form.get('method', 'global')
        comparer = SequenceComparer(fasta_manager.get_sequences())
        result = comparer.compare_pair(idx1, idx2, method)
    return render_template('align.html', result=result, names=names)

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

if __name__ == '__main__':
    app.run(debug=True)