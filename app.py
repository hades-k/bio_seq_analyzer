from flask import Flask, render_template, request, redirect, url_for, flash, send_file, Response
from werkzeug.utils import secure_filename
import os
import pandas as pd
from tools import Parser, SequenceAligner, MotifFinder
from sequence import MitochondrialDNA
from comparer import SequenceComparer, ConservedMotifAnalyzer
import io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = {'fasta', 'fa', 'fna'}

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = 'bioseq2024'  # For flash messages

if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)

# --- OOP Web Controller Classes ---
class FastaManager:
    def __init__(self):
        self.files = []  # List of (filename, DataFrame, mito_objs)

    def parse(self, filepath, filename):
        parser = Parser('fasta')
        df = parser.run(filepath)
        mito_objs = [MitochondrialDNA(df.loc[i]) for i in range(len(df))]
        self.files.append((filename, df, mito_objs))
        return df

    def get_stats(self):
        if not self.files:
            return {}
        lengths = pd.concat([df['length'] for _, df, _ in self.files])
        stats = {
            'count': len(lengths),
            'min_length': lengths.min(),
            'max_length': lengths.max(),
            'mean_length': lengths.mean(),
        }
        return stats

    def get_gc_contents(self):
        if not self.files:
            return []
        return [obj.gc_content for _, _, mito_objs in self.files for obj in mito_objs]

    def get_names(self):
        if not self.files:
            return []
        return [obj._MitochondrialDNA__name for _, _, mito_objs in self.files for obj in mito_objs]

    def get_sequences(self):
        return [obj for _, _, mito_objs in self.files for obj in mito_objs]

    def get_file_labels(self):
        # Returns a list of (filename, [names])
        return [(filename, [obj._MitochondrialDNA__name for obj in mito_objs]) for filename, _, mito_objs in self.files]

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

# --- Flask Routes ---
fasta_manager = FastaManager()

# Store motif_counts for graphing
last_motif_counts = {}

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
            fasta_manager.parse(filepath, filename)
            flash(f'File {filename} uploaded and parsed successfully!')
            return redirect(url_for('summary'))
        else:
            flash('Invalid file type!')
            return redirect(request.url)
    return render_template('index.html', files=fasta_manager.get_file_labels())

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
    global last_motif_counts
    results = None
    discovered = None
    motif_counts = None
    k = 5
    threshold = 2
    if request.method == 'POST':
        motif = request.form.get('motif')
        k = int(request.form.get('k', 5))
        threshold = int(request.form.get('threshold', 2))
        tool = MotifWebTool(fasta_manager.get_sequences())
        if motif and motif.strip():
            results = tool.search_motif(motif)
            last_motif_counts = {}
        else:
            analyzer = ConservedMotifAnalyzer(fasta_manager.get_sequences())
            discovered, motif_counts = analyzer.find_conserved(k=k, threshold=threshold)
            last_motif_counts = motif_counts if motif_counts else {}
            if discovered is None:
                discovered = {}
    else:
        discovered = {}
        motif_counts = {}
        last_motif_counts = {}
    return render_template('motif.html', results=results, discovered=discovered, motif_counts=motif_counts, k=k, threshold=threshold)

@app.route('/motif_graph')
def motif_graph():
    global last_motif_counts
    motif_counts = last_motif_counts
    if not motif_counts:
        return Response(status=204)
    # Sort motifs by frequency
    motifs, counts = zip(*sorted(motif_counts.items(), key=lambda x: x[1], reverse=True))
    plt.figure(figsize=(10, 4))
    plt.bar(motifs, counts)
    plt.xlabel('Motif')
    plt.ylabel('Frequency')
    plt.title('Motif Frequency')
    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    plt.close()
    buf.seek(0)
    return Response(buf.getvalue(), mimetype='image/png')

@app.route('/gc_graph')
def gc_graph():
    gc_contents = fasta_manager.get_gc_contents()
    if not gc_contents:
        return Response(status=204)
    plt.figure(figsize=(8, 4))
    plt.hist([gc for gc in gc_contents if gc is not None], bins=20, color='skyblue', edgecolor='black')
    plt.xlabel('GC Content (%)')
    plt.ylabel('Number of Sequences')
    plt.title('GC Content Distribution')
    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    plt.close()
    buf.seek(0)
    return Response(buf.getvalue(), mimetype='image/png')

@app.route('/align', methods=['GET', 'POST'])
def align():
    result = None
    names = fasta_manager.get_names()
    seq1_name = seq2_name = method = None
    if request.method == 'POST':
        idx1 = int(request.form.get('seq1'))
        idx2 = int(request.form.get('seq2'))
        method = request.form.get('method', 'global')
        tool = AlignmentWebTool(fasta_manager.get_sequences())
        result = tool.align_pair(idx1, idx2, method)
        result['method'] = method.capitalize()
        result['seq1_name'] = names[idx1]
        result['seq2_name'] = names[idx2]
        result['aligned_seq1'] = ''.join(result['aligned_seq1'])
        result['aligned_seq2'] = ''.join(result['aligned_seq2'])
        result['matches'] = ''.join(result['matches'])
        seq1_name = names[idx1]
        seq2_name = names[idx2]
    return render_template('align.html', result=result, names=names, seq1_name=seq1_name, seq2_name=seq2_name, method=method)

if __name__ == '__main__':
    app.run(debug=True) 