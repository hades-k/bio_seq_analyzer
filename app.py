from flask import Flask, render_template, request, redirect, url_for, flash
from werkzeug.utils import secure_filename
import os
from tools import Parser
from sequence import MitochondrialDNA
from comparer import SequenceComparer
import matplotlib.pyplot as plt
import io
import base64
import pandas as pd
from alignment_window import AlignmentWindow

UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = {'fasta', 'fa', 'fna'}

print(">>> THIS FILE IS RUNNING <<<")

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.secret_key = 'bioseq2025'

if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)

class SequenceManager:
    def __init__(self):
        self.df = None
        self.sequences = []

    def load_fasta(self, file_path):
        parser = Parser()
        df = parser.run(file_path)
        if self.df is None:
            self.df = df
        else:
            self.df = pd.concat([self.df, df], ignore_index=True)
        new_seqs = [MitochondrialDNA(df.loc[i]) for i in range(len(df))]
        self.sequences.extend(new_seqs)

    def get_gc_plot(self):
        if not self.sequences:
            return None
        names = [seq._MitochondrialDNA__name for seq in self.sequences]
        values = [seq.gc_content for seq in self.sequences]
        plt.clf()
        plt.figure(figsize=(10, 4))
        plt.bar(names, values)
        plt.xticks(rotation=90)
        plt.title("GC Content per Sequence")
        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        return base64.b64encode(buf.getvalue()).decode()

    def get_motif_plot(self, motif):
        if not self.sequences:
            return None
        counts = []
        names = []
        for seq in self.sequences:
            positions = seq.sequence.count(motif)
            counts.append(positions)
            names.append(seq._MitochondrialDNA__name)
        plt.clf()
        plt.figure(figsize=(10, 4))
        plt.bar(names, counts)
        plt.xticks(rotation=90)
        plt.title(f"Occurrences of motif '{motif}'")
        plt.tight_layout()
        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        return base64.b64encode(buf.getvalue()).decode()

    def get_sequence_names(self):
        return [seq._MitochondrialDNA__name for seq in self.sequences]

    def get_sequences(self):
        return self.sequences

seq_manager = SequenceManager()

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        files = request.files.getlist('file')
        if not files or files[0].filename == '':
            flash('No selected file')
            return redirect(request.url)
        for file in files:
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                file.save(file_path)
                seq_manager.load_fasta(file_path)
        flash('FASTA file(s) uploaded and parsed successfully!')
        return redirect(url_for('summary'))
    return render_template('index.html')

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/summary')
def summary():
    gc_plot = seq_manager.get_gc_plot()
    return render_template('summary.html', gc_plot=gc_plot)

@app.route('/motif', methods=['GET', 'POST'])
def motif():
    motif_plot = None
    if request.method == 'POST':
        motif = request.form.get('motif')
        if motif:
            motif_plot = seq_manager.get_motif_plot(motif)
    return render_template('motif.html', motif_plot=motif_plot)



from alignment_window import AlignmentWindow  # ⬅ make sure this import is at the top

@app.route('/alignment', methods=['GET', 'POST'])
def alignment():
    names = seq_manager.get_sequence_names()
    info = {}
    result_block = ""
    slider = 0
    max_slider = 0

    if request.method == 'POST':
        idx1 = int(request.form['seq1'])
        idx2 = int(request.form['seq2'])
        method = request.form['method']
        slider = int(request.form.get('slider', 0))

        comparer = SequenceComparer(seq_manager.get_sequences())
        raw = comparer.compare_pair(idx1, idx2, method)

        win = AlignmentWindow(
            seq1_name=names[idx1],
            seq2_name=names[idx2],
            aligned_seq1=raw['aligned_seq1'],
            aligned_seq2=raw['aligned_seq2'],
            matches=raw['matches'],
            score=raw['score'],
            match_count=raw['match_count'],
            mismatch_count=raw['mismatch_count'],
            gap_count=raw['gap_count']
        )

        result = win.get_window(slider)
        result_block = result.window_text
        info = {
            'score': result.score,
            'matches': result.matches,
            'mismatches': result.mismatches,
            'gaps': result.gaps,
            'name1': result.name1,
            'name2': result.name2
        }
        max_slider = result.max_slider

    return render_template("alignment.html", names=names, result_block=result_block,
                           info=info, slider=slider, max_slider=max_slider)


if __name__ == '__main__':
    app.run(debug=True)
