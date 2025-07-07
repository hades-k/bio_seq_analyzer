from flask import Flask, render_template, request, redirect, url_for, flash, send_file
from werkzeug.utils import secure_filename
import os
import io
import matplotlib.pyplot as plt
from tools import Parser, MotifFinder
from sequence import MitochondrialDNA
from comparer import SequenceComparer

UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = {'fasta', 'fa', 'fna'}

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = 'bioseq2024'

if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)

class FastaManager:
    def __init__(self):
        self.datasets = {}  # {filename: list of MitochondrialDNA objects}
        self.current_file = None

    def parse(self, filepath):
        parser = Parser('fasta')
        df = parser.run(filepath)
        mito_objs = [MitochondrialDNA(df.loc[i]) for i in range(len(df))]
        filename = os.path.basename(filepath)
        self.datasets[filename] = {
            'df': df,
            'sequences': mito_objs
        }
        self.current_file = filename
        return df

    def set_current_file(self, filename):
        if filename in self.datasets:
            self.current_file = filename

    def get_stats(self):
        if not self.current_file:
            return {}
        df = self.datasets[self.current_file]['df']
        lengths = df['length']
        return {
            'count': len(df),
            'min_length': lengths.min(),
            'max_length': lengths.max(),
            'mean_length': lengths.mean()
        }

    def get_gc_contents(self):
        return [obj.gc_content for obj in self.datasets[self.current_file]['sequences']]

    def get_names(self):
        return list(self.datasets[self.current_file]['df']['name'])

    def get_sequences(self):
        return self.datasets[self.current_file]['sequences']

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
    return render_template('index.html', fasta_manager=fasta_manager)

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/summary')
def summary():
    stats = fasta_manager.get_stats()
    gc_contents = fasta_manager.get_gc_contents()
    names = fasta_manager.get_names()
    return render_template('summary.html', stats=stats, names=names, gc_contents=gc_contents, fasta_manager=fasta_manager, zip=zip)

@app.route('/set_file', methods=['POST'])
def set_file():
    filename = request.form.get('filename')
    fasta_manager.set_current_file(filename)
    return redirect(url_for('summary'))

@app.route('/motif', methods=['GET', 'POST'])
def motif():
    results = None
    discovered = None
    finder = MotifFinder()
    sequences = fasta_manager.get_sequences()

    if request.method == 'POST':
        motif_str = request.form.get('motif')
        k = int(request.form.get('k', 5))
        min_sequences_with_motif = int(request.form.get('threshold', 2))

        if motif_str:
            results = finder.run(sequences, motif=motif_str)
            fasta_manager.motif_results = results
        else:
            discovered = finder.run(sequences, k=k, threshold=min_sequences_with_motif)

    return render_template('motif.html', results=results,
                           discovered=discovered, sequence_names=fasta_manager.get_names())

@app.route('/align', methods=['GET', 'POST'])
def align():
    result = None
    names = fasta_manager.get_names()
    if request.method == 'POST':
        idx1 = int(request.form.get('seq1'))
        idx2 = int(request.form.get('seq2'))
        method = request.form.get('method', 'global')

        # Get custom scoring
        match = int(request.form.get('match', 2))
        mismatch = int(request.form.get('mismatch', -1))
        gap = int(request.form.get('gap', -2))

        from tools import SequenceAligner
        aligner = SequenceAligner(match=match, mismatch=mismatch, gap=gap)
        seq1 = fasta_manager.get_sequences()[idx1].sequence
        seq2 = fasta_manager.get_sequences()[idx2].sequence
        aligner.run(seq1, seq2, method)
        result = aligner.get_alignment_data()
        result['seq1_name'] = names[idx1]
        result['seq2_name'] = names[idx2]

        # Pre-join for template
        result['aligned_seq1'] = ''.join(result['aligned_seq1'])
        result['aligned_seq2'] = ''.join(result['aligned_seq2'])
        result['matches'] = ''.join(result['matches'])

    return render_template('align.html', result=result, names=names, indexed_names=list(enumerate(names)))

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
    plt.xticks()
    plt.tight_layout()
    img = io.BytesIO()
    fig.savefig(img, format='png')
    img.seek(0)
    return send_file(img, mimetype='image/png')

@app.route('/motif_histogram.png')
def motif_histogram_png():
    if not fasta_manager.motif_results:
        return "", 404
    names = [r['sequence_name'] for r in fasta_manager.motif_results]
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
