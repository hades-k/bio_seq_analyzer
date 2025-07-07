{# templates/summary.html #}
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Summary</title>
</head>
<body>
    <h1>FASTA Summary</h1>
    <p><strong>Total sequences:</strong> {{ stats.count }}</p>
    <p><strong>Min length:</strong> {{ stats.min_length }}</p>
    <p><strong>Max length:</strong> {{ stats.max_length }}</p>
    <p><strong>Average length:</strong> {{ stats.mean_length }}</p>

    <h2>GC Content</h2>
    <ul>
        {% for name, gc in zip(names, gc_contents) %}
            <li>{{ name }}: {{ gc|round(2) }}%</li>
        {% endfor %}
    </ul>
</body>
</html>
