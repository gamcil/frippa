def get_maximum_row_lengths(rows):
    """Finds the longest lengths of fields per column in a collection of rows."""
    lengths, total = [], len(rows[0])
    for index in range(total):
        largest = max(len(row[index]) for row in rows)
        lengths.append(largest)
    return lengths


def add_field_whitespace(rows, lengths):
    """Fills table fields with whitespace to specified lengths."""
    result = []
    for row in rows:
        fmt = [f"{row[index]:{length}}" for index, length in enumerate(lengths)]
        result.append(fmt)
    return result


def humanise(rows):
    """Formats a collection of fields as human-readable."""
    lengths = get_maximum_row_lengths(rows)
    table = add_field_whitespace(rows, lengths)
    return table


def check_index(index, cluster):
    return [
        "X" if index in array else "-"
        for array in (cluster.duf, cluster.precursors, cluster.signalp)
    ]


def format_cluster(cluster, delimiter=None, show_headers=True):
    rows = []

    if show_headers:
        headers = [
            "Protein",
            "Start",
            "End",
            "Strand",
            "DUF3328",
            "Precursor",
            "Signal Peptide",
        ]
        rows.append(headers)

    for index, protein in enumerate(cluster.proteins):
        row = [
            protein.name,
            str(protein.start),
            str(protein.end),
            "+" if protein.strand == 1 else "-",
            *check_index(index, cluster)
        ]
        rows.append(row)

    rows = humanise(rows)

    repeats = "\n\n".join(
        cluster.proteins[idx].name
        + "\n"
        + cluster.proteins[idx].summarise_repeats()
        for idx in cluster.precursors
    ) or "No repeats found..."

    if not delimiter:
        delimiter = "  "

    rows = [delimiter.join(row) for row in rows]

    if show_headers:
        separator = "-" * len(rows[0])
        rows[0] = "{}\n{}\n{}".format(separator, rows[0], separator)

    return "{}\n\n{}\n\n{}\n{}".format(
        cluster,
        repeats,
        "\n".join(rows),
        separator,
    )
