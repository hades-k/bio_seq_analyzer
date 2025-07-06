from flask import Flask, render_template_string, request, redirect, url_for, flash

app = Flask(__name__)
app.secret_key = 'bioseq2024'

# --- Template Strings ---
base_html_py="""<!doctype html>
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

summary_html_py = """<h2>Summary Statistics</h2>
<ul class="list-group">
  <li class="list-group-item">Total Sequences: {{ stats.count }}</li>
  <li class="list-group-item">Min Length: {{ stats.min_length }}</li>
  <li class="list-group-item">Max Length: {{ stats.max_length }}</li>
  <li class="list-group-item">Mean Length: {{ stats.mean_length | round(2) }}</li>
</ul>
<h3 class="mt-4">GC Content (%)</h3>
<table class="table table-striped">
  <thead><tr><th>Species</th><th>GC</th></tr></thead>
  <tbody>
    {% for name, gc in zip(names, gc_contents) %}
      <tr><td>{{ name }}</td><td>{{ gc | round(2) }}</td></tr>
    {% endfor %}
  </tbody>
</table>"""

# --- Dummy Backend ---
class FastaManager:
    def get_stats(self):
        return {'count': 5, 'min_length': 800, 'max_length': 1400, 'mean_length': 1023.6}
    def get_gc_contents(self):
        return [42.1, 43.2, 41.8, 44.5, 45.1]
    def get_names(self):
        return ['Species A', 'Species B', 'Species C', 'Species D', 'Species E']

fasta_manager=FastaManager()

# --- Routes ---
@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        flash("File uploaded successfully!")  # Stub action
        return redirect(url_for('summary'))
    return render_template_string(base_html_py, content=index_html_py)

@app.route('/summary')
def summary():
    stats=fasta_manager.get_stats()
    gc_contents=fasta_manager.get_gc_contents()
    names=fasta_manager.get_names()
    return render_template_string(
        base_html_py,
        content=render_template_string(summary_html_py, stats=stats, gc_contents=gc_contents, names=names)
    )

if __name__=='__main__':
    app.run(debug=True)
