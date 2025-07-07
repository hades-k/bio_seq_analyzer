<!DOCTYPE html>
<html>
<head>
    <title>Alignment Viewer</title>
    <style>
        body {
            font-family: monospace;
            padding: 20px;
        }
        pre {
            background: #f8f8f8;
            padding: 10px;
            border: 1px solid #ddd;
            white-space: pre;
        }
    </style>
</head>
<body>
    <h1>Pairwise Sequence Alignment</h1>

    <form method="POST">
        <label>Sequence 1:</label>
        <select name="seq1">
            {% for i in range(names|length) %}
                <option value="{{ i }}">{{ names[i] }}</option>
            {% endfor %}
        </select>

        <label>Sequence 2:</label>
        <select name="seq2">
            {% for i in range(names|length) %}
                <option value="{{ i }}">{{ names[i] }}</option>
            {% endfor %}
        </select>

        <label>Method:</label>
        <select name="method">
            <option value="global">Global</option>
            <option value="local">Local</option>
        </select>

        <label>Alignment Window:</label>
        <input type="range" name="slider" min="0" max="{{ max_slider }}" value="{{ slider }}" oninput="this.form.submit()">
    </form>

    {% if result_block %}
    <h2>Result</h2>
    <p><strong>Score:</strong> {{ info.score }} |
       <strong>Matches:</strong> {{ info.matches }} |
       <strong>Mismatches:</strong> {{ info.mismatches }} |
       <strong>Gaps:</strong> {{ info.gaps }}</p>
    <pre>{{ result_block }}</pre>
    {% endif %}
</body>
</html>
