{# templates/motif.html #}
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Motif Analysis</title>
</head>
<body>
    <h1>Motif Search</h1>
    <form method="POST">
        <label>Enter motif (or leave blank for discovery):</label>
        <input type="text" name="motif">
        <label>k (for k-mer):</label>
        <input type="number" name="k" value="5">
        <label>Threshold:</label>
        <input type="number" name="threshold" value="2">
        <input type="submit" value="Search">
    </form>

    {% if results %}
        <h2>Motif Occurrences</h2>
        <ul>
            {% for res in results %}
                <li>{{ res.name }} — Count: {{ res.count }}, Positions: {{ res.positions }}</li>
            {% endfor %}
        </ul>
    {% endif %}

    {% if discovered %}
        <h2>Discovered Motifs</h2>
        <ul>
            {% for motif, count in discovered.items() %}
                <li>{{ motif }} — {{ count }}</li>
            {% endfor %}
        </ul>
    {% endif %}
</body>
</html>
