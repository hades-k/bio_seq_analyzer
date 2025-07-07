{# templates/align.html #}
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Pairwise Alignment</title>
</head>
<body>
    <h1>Pairwise Sequence Alignment</h1>
    <form method="POST">
        <label for="seq1">Select Sequence 1:</label>
        <select name="seq1" required>
            {% for idx, name in enumerate(names) %}
                <option value="{{ idx }}">{{ name }}</option>
            {% endfor %}
        </select>

        <label for="seq2">Select Sequence 2:</label>
        <select name="seq2" required>
            {% for idx, name in enumerate(names) %}
                <option value="{{ idx }}">{{ name }}</option>
            {% endfor %}
        </select>

        <label for="method">Alignment Method:</label>
        <select name="method">
            <option value="global">Global</option>
            <option value="local">Local</option>
        </select>

        <input type="submit" value="Align">
    </form>

    {% if result %}
    <h2>Alignment Result ({{ result.type|capitalize }}):</h2>
    <p><strong>Score:</strong> {{ result.score }}</p>
    <p><strong>Matches:</strong> {{ result.match_count }} | <strong>Mismatches:</strong> {{ result.mismatch_count }} | <strong>Gaps:</strong> {{ result.gap_count }}</p>
    <hr>
    <pre style="font-family: monospace; white-space: pre-wrap;">
{% for i in range(0, result.aligned_seq1|length, 80) %}
{{ result.aligned_seq1[i:i+80] }}
{{ result.matches[i:i+80] }}
{{ result.aligned_seq2[i:i+80] }}

{% endfor %}
    </pre>
    {% endif %}
</body>
</html>