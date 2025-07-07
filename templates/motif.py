{# templates/motif.html #}
<!DOCTYPE html>
<html>
<head>
    <title>Motif Analysis</title>
</head>
<body>
    <h1>Motif Search</h1>
    <form method="POST">
        <label for="motif">Enter motif:</label>
        <input type="text" name="motif">
        <input type="submit" value="Search">
    </form>
    {% if motif_plot %}
        <h2>Motif Occurrence Chart</h2>
        <img src="data:image/png;base64,{{ motif_plot }}"/>
    {% endif %}
</body>
</html>