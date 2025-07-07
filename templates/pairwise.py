{# templates/pairwise.html #}
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>All Pairwise Alignments</title>
</head>
<body>
    <h1>All Pairwise Sequence Alignments</h1>
    <p>Total comparisons: {{ results|length }}</p>

    {% for res in results %}
        <h2>{{ names[res.idx1] }} vs {{ names[res.idx2] }} ({{ res.type|capitalize }})</h2>
        <p><strong>Score:</strong> {{ res.score }}</p>
        <p><strong>Matches:</strong> {{ res.match_count }} | <strong>Mismatches:</strong> {{ res.mismatch_count }} | <strong>Gaps:</strong> {{ res.gap_count }}</p>
        <pre style="font-family: monospace; white-space: pre-wrap;">
{% for i in range(0, res.aligned_seq1|length, 80) %}
{{ res.aligned_seq1[i:i+80] }}
{{ res.matches[i:i+80] }}
{{ res.aligned_seq2[i:i+80] }}

{% endfor %}
        </pre>
        <hr>
    {% endfor %}
</body>
</html>
