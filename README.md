# Causal Sentence Extraction from PubMed Papers

A Python command-line tool for retrieving **PubMed abstracts** and extracting sentences that mention multiple user-defined biomedical concepts alongside **causal or associative language**.

The extractor combines domain-specific synonym matching, cue-phrase detection, and a token-distance **window-hit filter** to identify potentially relevant statements for literature review and evidence screening. It also records article metadata, estimates study design, and extracts reported sample sizes when possible.

## Features

- **PubMed search:** Search via a query string or a `.txt` query file, or supply a list of PubMed IDs (PMIDs).
- **Abstract-only extraction:** Retrieve titles, abstracts, publication types, journal information, DOIs, and article links using NCBI E-utilities.
- **Custom domain groups:** Define biomedical concepts and their synonyms directly in `pubmed_extractor.py`.
- **Sentence-level filtering:** Require mentions from a specified number of distinct domain groups and at least one causal or associative cue.
- **Proximity-aware filtering:** Require a cue within or near a short token window spanning tokens from two different domain groups.
- **Language signals:** Capture strong causal expressions, associative expressions, hedging terms, and negation terms. Negated sentences are **retained**, not discarded.
- **Study design labeling:** Use PubMed publication types first, with keyword-based fallbacks.
- **Sample size extraction:** Identify common sample-size expressions in titles or abstracts.
- **CSV output:** Export matched sentences with citation links and supporting metadata.
- **Optional spaCy integration:** Use a spaCy English model for sentence splitting; otherwise use a regular-expression fallback.

## How it works

```text
PubMed query / query file / PMID list
                |
                v
       Retrieve PubMed records
                |
                v
       Split abstracts into sentences
                |
                v
       Match domain-group synonyms
                |
                v
       Detect causal / associative cues
                |
                v
       Apply optional window-hit filter
                |
                v
  Add study design + sample-size metadata
                |
                v
             Export CSV
```

By default, a sentence is included only when **all** of these conditions are met:

1. It matches synonyms from at least **2 distinct domain groups**.
2. It contains at least **1 strong causal or associative cue**.
3. It passes the **window-hit check**: tokens associated with two different domain groups appear no more than **6 tokens apart**, with a causal/associative cue between or near them (up to **2 tokens of slack** on either side).

These thresholds are configurable in the script. The window-hit check uses a token vocabulary built from group synonyms; it is a heuristic proximity filter, not a syntactic or semantic relationship parser.

## Requirements

- Python 3.8+ (recommended)
- Internet access to NCBI PubMed E-utilities
- [`requests`](https://pypi.org/project/requests/)
- **Optional:** [`spaCy`](https://spacy.io/) with the `en_core_web_sm` English model

The tool also uses Python standard-library modules; no other third-party packages are required by the script.

## Installation

1. Download or clone this repository and enter its directory.
2. Install the required dependency:

   ```bash
   python -m pip install requests
   ```

3. (Optional) Install spaCy and its English model for improved sentence splitting:

   ```bash
   python -m pip install spacy
   python -m spacy download en_core_web_sm
   ```

If spaCy or the English model is unavailable, the script automatically falls back to regular-expression sentence splitting.

## Configure domain groups

**Before running the extractor, review `DOMAIN_GROUPS` in `pubmed_extractor.py`.** Only uncommented entries are active. In the supplied version, the active groups are `AGE` and `DEATH`; other biomedical groups are provided as commented-out examples.

Synonym matching is case-insensitive and supports single-word terms, multiword phrases with spaces/hyphens, and a trailing `*` for prefix matching (for example, `outcome*`).

Update the filter settings in the same file as needed:

```python
MIN_DISTINCT_GROUPS = 2
WINDOW_SIZE = 6
REQUIRE_WINDOW_HIT = True
CUE_SLACK_TOKENS = 2
```

| Setting | Default | Meaning |
| --- | --- | --- |
| `MIN_DISTINCT_GROUPS` | `2` | Minimum number of distinct domain groups matched in a sentence |
| `WINDOW_SIZE` | `6` | Maximum token-index distance between tokens from two groups |
| `REQUIRE_WINDOW_HIT` | `True` | Whether the proximity-and-cue condition is mandatory |
| `CUE_SLACK_TOKENS` | `2` | Extra tokens allowed around the matched group-token pair when locating a cue |

At least as many domain groups must be enabled as `MIN_DISTINCT_GROUPS`, or the script exits with an error.

## Usage

### 1. Read a PubMed query from a text file

Create a `query.txt` file containing a PubMed search expression. Example queries are available in the folder named `sample_queries`. **NOTE: Clinical Context:** Additional oxygenation-related and ICU-specific keywords were included in the PubMed search query to restrict the literature to the clinical context of oxygen therapy and critical care. Although each time the primary relationship of interest is between the variables of concer, such as **age and mortality** in the current code, these additional terms helped retrieve studies relevant to oxygenation and mortality outcomes among critically ill patients.

Run:

```bash
python pubmed_extractor.py --query-file query.txt --max-articles 10000 --out results.csv
```

### 2. Use an NCBI API key (optional)

```bash
python pubmed_extractor.py \
  --query-file query.txt \
  --api-key YOUR_NCBI_API_KEY \
  --out results.csv
```

Treat your API key as a secret: do not commit it to your repository. Observe [NCBI E-utilities usage guidelines](https://www.ncbi.nlm.nih.gov/books/NBK25497/) when running larger searches.

### Command-line options

| Argument | Description | Default |
| --- | --- | --- |
| `--query` | PubMed query string, or path to an existing `.txt` query file | Empty |
| `--query-file` | Path to a UTF-8 text file containing a PubMed query | Empty |
| `--pmids` | Comma-separated PMIDs; bypasses query search | Empty |
| `--max-articles` | Maximum PMIDs requested from PubMed search | `200` |
| `--api-key` | Optional NCBI API key | Empty |
| `--out` | Path to the output CSV file | `results.csv` |

At least one of `--query`, `--query-file`, or `--pmids` must be provided.

## Output format

The script writes a UTF-8 CSV file containing **one row per retained sentence**, with these columns:

| Column | Description |
| --- | --- |
| `edge` | Semicolon-separated names of the domain groups detected in the sentence; **not** a verified or directed causal edge |
| `pmid` | PubMed ID of the source article |
| `study_design` | Heuristically assigned study category, such as `RCT`, `Cohort`, `Case-Control`, `Systematic Review`, `Meta-analysis`, `Case Report`, or `Other` |
| `sample_size` | Sample size inferred from recognized text patterns, or blank if none is found |
| `strong_cues` | Matched terms from the strong causal-cue dictionary |
| `assoc_cues` | Matched terms from the associative-cue dictionary |
| `hedge_or_neg` | Matched hedging and negation expressions |
| `sentence` | Extracted abstract sentence |
| `url` | DOI link when available, otherwise the PubMed article link |

A CSV containing only the header row is valid when no sentences pass the filters. Progress and row counts are printed to the terminal.

### Study design classification

The script first maps publication types returned by PubMed into a predefined set of categories. If no mapped publication type is present, it looks for study-design keywords in the title and abstract. This is a convenient **screening label**, not a validated methodological assessment. Some publication types, including clinical-trial phase labels, are deliberately mapped to `Other` in the current implementation.

### Sample size extraction

The script recognizes patterns such as `n=250`, `patients (n=250)`, and `250 patients`. When several candidate numbers occur, it uses the **largest match**. This can misidentify the actual analytic sample size, so verify it against the full paper before using it in evidence weighting.

## Important Points

- **Abstracts only:** The extractor does not download or examine full-text articles.
- **No causal verification:** Strong causal wording and associative wording are recorded separately; neither guarantees a true causal effect.
- **Heuristic proximity:** A window hit may link unrelated terms, while valid relationships expressed across longer spans may be missed.
- **Negation and hedging:** These are flagged in `hedge_or_neg`, but do not prevent a sentence from being included.
- **Study design and sample size:** Both are heuristic, and their accuracy depends on the available metadata and text.
- **PubMed XML parsing:** The implementation uses regular expressions to extract fields from returned XML and may not capture every complex or nested formatting case.
- **Rate limits and connectivity:** Requests require network access and can fail when PubMed is unavailable or limits are exceeded. Fetch requests are made in batches of up to 200 PMIDs, with a short pause between batches.

## Project files

```text
pubmed_extractor.py       # PubMed search, sentence filtering, and CSV export
query.txt        # User-created PubMed query file
README.md        # Project documentation
```

## Intended use

This tool is designed to support **biomedical literature screening**, **causal-knowledge discovery research**, and the collection of **candidate evidence sentences** for subsequent human review or downstream analysis.

The extracted sentences and labels should be treated as leads for further investigation rather than clinical conclusions, validated evidence grades, or ground-truth causal relationships.


