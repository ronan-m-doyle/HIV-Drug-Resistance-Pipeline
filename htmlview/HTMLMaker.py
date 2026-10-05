import json
import re
import csv
from collections import defaultdict
from pathlib import Path
import argparse

# ---------------------------------------------------------------------------
# Drug code -> full generic name lookup.
# The Sierra JSON only returns the short code (e.g. "ATV") and the clinical
# abbreviation (e.g. "ATV/r"); it doesn't return the full generic name, so we
# maintain this small fixed table ourselves. Falls back to the short code if
# an unrecognized drug ever appears (e.g. a newly added drug not yet listed
# here).
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Gene code -> full region name lookup, for the scored genes that appear in
# drugResistance blocks (PR/RT/IN, and CA for capsid inhibitors on newer
# HIVDB algorithm versions). Falls back to the raw code if an unrecognized
# gene ever appears.
# ---------------------------------------------------------------------------
GENE_FULL_NAMES = {
    'PR': 'Protease',
    'RT': 'Reverse Transcriptase',
    'IN': 'Integrase',
    'CA': 'Capsid',
}

DRUG_FULL_NAMES = {
    # PIs
    'ATV': 'Atazanavir', 'DRV': 'Darunavir', 'FPV': 'Fosamprenavir',
    'IDV': 'Indinavir', 'LPV': 'Lopinavir', 'NFV': 'Nelfinavir',
    'SQV': 'Saquinavir', 'TPV': 'Tipranavir',
    # NRTIs
    'ABC': 'Abacavir', 'AZT': 'Zidovudine', 'D4T': 'Stavudine',
    'DDI': 'Didanosine', 'FTC': 'Emtricitabine', 'ISL': 'Islatravir',
    'LMV': 'Lamivudine', 'TDF': 'Tenofovir',
    # NNRTIs
    'DOR': 'Doravirine', 'DPV': 'Dapivirine', 'EFV': 'Efavirenz',
    'ETR': 'Etravirine', 'NVP': 'Nevirapine', 'RPV': 'Rilpivirine',
    # INSTIs
    'BIC': 'Bictegravir', 'CAB': 'Cabotegravir', 'DTG': 'Dolutegravir',
    'EVG': 'Elvitegravir', 'RAL': 'Raltegravir',
}

# Standard codon -> amino acid translation table (single-letter code, '*' = stop)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

MUTATION_RE = re.compile(r'^[A-Za-z\*]*?(\d+)([A-Za-z\*]+)$')

# ---------------------------------------------------------------------------
# Top-level "analysis options" fields present in seqreads report JSON.
# All four rate/prevalence fields (minPrevalence, actualMinPrevalence,
# maxMixtureRate, mixtureRate) are fractions (0-1 scale) and are multiplied
# by 100 for display -- confirmed: 0.1 means 10%, 0.02 means 2%, etc.
# ---------------------------------------------------------------------------
ANALYSIS_PARAM_FIELDS = [
    ('minPositionReads', 'Minimum read depth per position', None, 1),
    ('minCodonReads', 'Minimum codon reads', None, 1),
    ('minPrevalence', 'Minimum prevalence threshold (requested)', '%', 100),
    ('actualMinPrevalence', 'Minimum prevalence threshold (actual)', '%', 100),
    ('maxMixtureRate', 'Maximum mixture rate (threshold)', '%', 100),
    ('mixtureRate', 'Observed mixture rate', '%', 100),
]


def format_analysis_params(entry):
    """
    Pull out the seqreads analysis-options fields present on the entry
    (if any) as a list of (label, formatted_value) pairs.
    """
    rows = []
    for key, label, suffix, scale in ANALYSIS_PARAM_FIELDS:
        if key not in entry:
            continue
        value = entry[key]
        if isinstance(value, float):
            value_str = "{:.3f}".format(value * scale)
        else:
            value_str = str(value)
        if suffix:
            value_str += suffix
        rows.append((label, value_str))
    return rows


def build_gene_position_ranges(entry):
    """
    Map gene code -> (firstAA, lastAA) using the seqreads entry's
    allGeneSequenceReads block, so coverage plots/tables can span the full
    sequenced range of each gene (including gaps, which get 0 depth).
    """
    ranges = {}
    for g in entry.get('allGeneSequenceReads', []):
        ranges[g['gene']['name']] = (g['firstAA'], g['lastAA'])
    return ranges


def get_position_depths(gene_code, first_aa, last_aa, codfreq_table):
    """Return (positions, depths) lists for a gene's full sequenced range."""
    positions = list(range(first_aa, last_aa + 1))
    depths = [codfreq_table.get((gene_code, p), {}).get('total_reads', 0) for p in positions]
    return positions, depths


def compute_coverage_stats(depths, min_reads=None):
    """Summary stats for a gene's per-position depth list."""
    n = len(depths)
    if n == 0:
        return None
    sorted_depths = sorted(depths)
    mean_depth = sum(depths) / n
    median_depth = sorted_depths[n // 2] if n % 2 else \
        (sorted_depths[n // 2 - 1] + sorted_depths[n // 2]) / 2
    stats = {
        'n_positions': n,
        'min': min(depths),
        'max': max(depths),
        'mean': mean_depth,
        'median': median_depth,
    }
    if min_reads:
        below = sum(1 for d in depths if d < min_reads)
        stats['pct_below_threshold'] = 100.0 * below / n
    return stats


def render_coverage_table(gene_stats_rows, min_reads=None):
    """
    gene_stats_rows: list of (gene_label, stats_dict) tuples.
    Renders one small summary table across all genes.
    """
    threshold_col = "<th>% Positions &lt; Min Depth</th>" if min_reads else ""
    rows_html = []
    for gene_label, stats in gene_stats_rows:
        if stats is None:
            continue
        threshold_cell = ""
        if min_reads:
            threshold_cell = "<td>{:.1f}%</td>".format(stats.get('pct_below_threshold', 0.0))
        rows_html.append(
            "<tr><td>{}</td><td>{}</td><td>{:,}</td><td>{:,.0f}</td>"
            "<td>{:,.0f}</td><td>{:,}</td>{}</tr>".format(
                gene_label, stats['n_positions'], stats['min'],
                stats['mean'], stats['median'], stats['max'], threshold_cell
            )
        )
    if not rows_html:
        return ''
    return (
        "<table class='table table-striped table-condensed' style='max-width:640px;'>"
        "<thead><tr><th>Gene</th><th>Positions</th><th>Min</th><th>Mean</th>"
        "<th>Median</th><th>Max</th>{}</tr></thead>"
        "<tbody>{}</tbody></table>".format(threshold_col, "".join(rows_html))
    )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('json', type=argparse.FileType('r'), help='Input JSON file.')
    parser.add_argument('-o', '--output', default='results.html',
                         help='Name of the output HTML file (default: results.html).')
    parser.add_argument('--codfreq', default=None,
                         help='Path to the original .codfreq file, used to add per-position '
                              'read-support (relative abundance + coverage) to the report. '
                              'If omitted, the script tries to auto-locate it using the '
                              '"name" field embedded in seqreads JSON output.')
    args = parser.parse_args()
    return args


def translate_codon(codon):
    """
    Translate a single codon token to an amino acid. Handles the standard
    3-nt codons; anything else (indel markers containing '-', insertions
    longer than 3nt, ambiguous/partial codons) is bucketed as a labelled
    non-substitution category rather than silently dropped or shown as '?'.
    """
    codon = codon.upper().replace('U', 'T')
    if '-' in codon:
        return 'del'          # partial/complete deletion within this codon
    if len(codon) > 3:
        return 'ins'          # insertion (extra nucleotides read at this codon)
    if len(codon) < 3:
        return '?'            # incomplete/partial read, can't translate
    return CODON_TABLE.get(codon, '?')


def load_codfreq_table(path):
    """
    Parse a .codfreq file into:
        {(gene, position): {'total_reads': int, 'aa_counts': {AA: count}}}
    The file is comma-delimited with a header row:
        gene,position,total,codon,count,total_quality_score
    and may start with a UTF-8 byte-order-mark. Handles files with or
    without the header row.
    """
    table = defaultdict(lambda: {'total_reads': 0, 'aa_counts': defaultdict(int)})
    with open(path, encoding='utf-8-sig', newline='') as fh:
        reader = csv.reader(fh)
        for fields in reader:
            if len(fields) < 5:
                continue
            gene, pos, total_reads, codon, count = [f.strip() for f in fields[:5]]
            if not pos.isdigit():
                continue  # header row, skip
            pos = int(pos)
            entry = table[(gene, pos)]
            entry['total_reads'] = int(total_reads)
            entry['aa_counts'][translate_codon(codon)] += int(count)
    return table


def format_read_support(gene, mutation_text, codfreq_table):
    """
    Given a gene (e.g. 'RT') and a mutation string like 'V106VI', look up the
    read support at that position for just the AA(s) named in the mutation
    call (here, V and I) and return a summary string like:
        "V:57.6% (28,777) I:40.6% (20,300) [coverage 49,980]"
    Low-frequency background noise at the position (sequencing errors, rare
    minority variants not part of the call) is intentionally excluded.
    Returns '' if no matching position data is found (e.g. codfreq
    unavailable, or the mutation is an insertion/deletion that doesn't map
    to a single tracked position).
    """
    if codfreq_table is None:
        return ''
    m = MUTATION_RE.match(mutation_text)
    if not m:
        return ''
    pos = int(m.group(1))
    called_aas = set(m.group(2).upper())  # e.g. {'V', 'I'} from "VI" in "V106VI"
    entry = codfreq_table.get((gene, pos))
    if not entry or entry['total_reads'] == 0:
        return ''
    total = entry['total_reads']
    parts = []
    for aa, count in sorted(entry['aa_counts'].items(), key=lambda kv: -kv[1]):
        if aa not in called_aas:
            continue
        pct = 100.0 * count / total
        parts.append("{}:{:.1f}% ({:,})".format(aa, pct, count))
    if not parts:
        return ''
    return " ".join(parts) + " [coverage {:,}]".format(total)


def main():
    args = parse_args()

    html_template = [
        "<html>",
        "<head>",
        '<link rel="stylesheet" href="css/bootstrap.css" crossorigin="anonymous">',
        "</head>",
        "<body>",
    ]

    js = json.load(args.json)
    if isinstance(js, dict):
        js = [js]  # seqreads output is a single object, not a list

    html_template += ["<div class='col-md-8'>"]
    for idx, entry in enumerate(js):
        # --- handle both fasta-style and seqreads-style shapes ---
        is_seqreads = 'inputSequence' not in entry

        if not is_seqreads:
            header = entry['inputSequence']['header']
            subtype = entry.get('subtypeText', '')
        else:
            raw_name = entry.get('name', 'Unknown sample')
            # keep just the filename, and drop a trailing .codfreq extension
            header = Path(raw_name).name
            if header.endswith('.codfreq'):
                header = header[:-len('.codfreq')]
            subtype = entry.get('bestMatchingSubtype', {}).get('display', '')

        validationResults = entry.get('validationResults', [])
        drugResistance = entry['drugResistance']  # list of per-gene blocks, same in both

        # --- locate / load the codfreq file for read-support data (seqreads only) ---
        codfreq_table = None
        if is_seqreads:
            codfreq_path = args.codfreq or entry.get('name')
            if codfreq_path and Path(codfreq_path).is_file():
                try:
                    codfreq_table = load_codfreq_table(codfreq_path)
                except (OSError, ValueError):
                    codfreq_table = None

        html_template += ["<div class='page-header' tabindex={}><h2>{}</h2>".format(idx+1, header)]

        if is_seqreads:
            analysis_rows = format_analysis_params(entry)
            if analysis_rows:
                html_template += ["<div class='analysis-params'><h4>Analysis Parameters</h4><ul>"]
                for label, value_str in analysis_rows:
                    html_template += ["<li><strong>{}:</strong> {}</li>".format(label, value_str)]
                html_template += ["</ul></div>"]

        html_template += ["<p><strong>Subtype:</strong> {}</p>".format(subtype)]

        if is_seqreads and codfreq_table is None:
            html_template += [
                "<div class='alert alert-info' role='alert'>"
                "Read-support data unavailable (original .codfreq file not found "
                "at the recorded path). Pass --codfreq to specify its location."
                "</div>"
            ]

        # --- validation / QC warnings block ---
        if validationResults:
            html_template += ["<div class='validation-results'>"]
            for v in validationResults:
                level = v.get('level', 'INFO')
                message = v.get('message', '')
                alert_class = {
                    'ERROR': 'alert-danger',
                    'WARNING': 'alert-warning',
                    'SEVERE_WARNING': 'alert-warning',
                }.get(level.upper(), 'alert-info')
                html_template += [
                    "<div class='alert {}' role='alert'><strong>{}:</strong> {}</div>".format(
                        alert_class, level, message
                    )
                ]
            html_template += ["</div>"]

        # --- coverage table (before the drug-resistance tables) ---
        if is_seqreads and codfreq_table is not None:
            gene_ranges = build_gene_position_ranges(entry)
            min_reads_threshold = entry.get('minPositionReads')
            table_rows = []
            for geneBlock in drugResistance:
                gcode = geneBlock['gene']['name']
                if gcode not in gene_ranges:
                    continue
                first_aa, last_aa = gene_ranges[gcode]
                gene_label = "{} ({})".format(GENE_FULL_NAMES.get(gcode, gcode), gcode)
                positions, depths = get_position_depths(gcode, first_aa, last_aa, codfreq_table)
                stats = compute_coverage_stats(depths, min_reads_threshold)
                table_rows.append((gene_label, stats))

            if table_rows:
                html_template += ["<h3>Coverage</h3>"]
                html_template += [render_coverage_table(table_rows, min_reads_threshold)]

        # one table per gene block (PR, RT, IN) instead of only drugResistance[0]
        for geneBlock in drugResistance:
            geneName = geneBlock['gene']['name']
            geneLabel = "{} ({})".format(GENE_FULL_NAMES.get(geneName, geneName), geneName)
            drugScores = geneBlock['drugScores']

            html_template += ["<h3>{}</h3>".format(geneLabel)]
            html_template += ["<table class='table table-striped'>"]
            header_cols = "<th>Drug</th><th>Score</th><th>Breakdown</th>"
            if codfreq_table is not None:
                header_cols += "<th>Read Support</th>"
            html_template += ['<thead><tr>{}</tr></thead>'.format(header_cols)]
            html_template += ['<tbody>']
            for drugScore in drugScores:
                drug_code = drugScore['drug']['name']
                drug_abbr = drugScore['drug']['displayAbbr']
                full_name = DRUG_FULL_NAMES.get(drug_code, drug_code)
                drug_label = "{} ({})".format(full_name, drug_abbr)

                html_template += ["<tr><td>{}</td>".format(drug_label)]
                html_template += ["<td>{} ({})</td>".format(str(drugScore['score']), drugScore['text'])]

                pscoretext = ''
                read_support_text = ''
                for partialScore in drugScore['partialScores']:
                    pscoretext += str(partialScore['score']) + ' '
                    mut_texts = [str(m['text']) for m in partialScore['mutations']]
                    pscoretext += '+'.join(mut_texts) + ' '

                    if codfreq_table is not None:
                        for mut_text in mut_texts:
                            support = format_read_support(geneName, mut_text, codfreq_table)
                            if support:
                                read_support_text += "{}: {}<br>".format(mut_text, support)

                html_template += ['<td>{}</td>'.format(pscoretext)]
                if codfreq_table is not None:
                    html_template += ['<td>{}</td>'.format(read_support_text)]
                html_template += ["</tr>"]
            html_template += ['</tbody></table>']

        html_template += ['</div>']

    html_template.append("</div>")
    html_template.append("<br>")
    html_template += ["</body>", "</html>"]

    with open(args.output, "w") as htmlout:
        htmlout.write("\n".join(html_template))

if __name__ == "__main__":
    main()