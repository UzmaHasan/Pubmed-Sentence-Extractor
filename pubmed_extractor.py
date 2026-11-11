#!/usr/bin/env python3
"""
Abstract-only causal sentence extractor (NO EDGE LIST) with MANUAL DOMAIN GROUPS + WINDOW HIT.

UPDATED:
- Supports --query-file query.txt (reads PubMed query from a file).
- If --query looks like a .txt file path and exists, it will auto-read it too.

What it does:
- Uses a PubMed query to retrieve PMIDs (or accepts explicit --pmids).
- Fetches title + abstract + publication types + venue + DOI/URL from PubMed (abstract-only).
- Splits abstracts into sentences.
- Keeps sentences that:
    (1) contain synonyms from at least MIN_DISTINCT_GROUPS distinct domain groups, AND
    (2) contain at least one causal/associative cue, AND
    (3) have WINDOW HIT (Option A): at least one pair of domain-group tokens within WINDOW_SIZE tokens
        AND a causal/associative cue occurs between (or very near) those tokens.
- Negation is NOT dropped; it is only flagged.

Usage:
  python Pubmed5.py --query-file query.txt --max-articles 10000 --out results.csv
  python Pubmed5.py --query "(age OR elderly) AND (death OR mortality)" --max-articles 250 --out results.csv
"""

import argparse
import csv
import html
import os
import re
import sys
import time
from typing import List, Dict, Tuple, Optional
from urllib.parse import urlencode

import requests

# ---------- Optional spaCy for better sentence splitting ----------
_SPACY_OK = False
try:
    import spacy  # type: ignore
    _nlp = spacy.load("en_core_web_sm", disable=["ner", "tagger"])
    _SPACY_OK = True
except Exception:
    _SPACY_OK = False

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# ================================================================
# 1) MANUAL DOMAIN GROUPS (EDIT THIS EVERY TIME YOU WANT)
# ================================================================
DOMAIN_GROUPS = {
    # "COPD": [
    #     "COPD",
    #     "chronic obstructive pulmonary disease",
    #     "chronic obstructive lung disease",
    #     "emphysema",
    #     "chronic bronchitis",
    # ],
    # "PACO2": [
    #     "PaCO2",
    #     "arterial carbon dioxide tension",
    #     "arterial partial pressure of carbon dioxide",
    #     "partial pressure of carbon dioxide",
    #     "carbon dioxide tension",
    #     "hypercapnia",
    #     "hypercapnoea",
    #     "CO2 retention",
    #     "carbon dioxide retention",
    #     "elevated carbon dioxide",
    #     "increased carbon dioxide",
    #     "high carbon dioxide",
    # ],
    # "ARDS": [
    #     "ARDS",
    #     "acute respiratory distress syndrome",
    #     "acute lung injury",
    #     "ALI",
    # ],
    # "PH": [
    #     "pH",
    #     "blood pH",
    #     "acidemia",
    #     "acidaemia",
    #     "alkalemia",
    #     "alkalaemia",
    #     "acid-base",
    #     "acid base",
    # ],
    # "LACTATE": [
    #     "lactate",
    #     "lactic acid",
    #     "hyperlactatemia",
    #     "hyperlactataemia",
    # ],
    # "BMI": [
    #     "BMI",
    #     "body mass index",
    # ],
    # "HEMOGLOBIN": [
    #     "hemoglobin",
    #     "haemoglobin",
    #     "Hb",
    #     "Hgb",
    # ],
    # "AGE": [
    #     "age",
    #     "older age",
    #     "advanced age",
    #     "elderly",
    #     "aging",
    #     "aged",
    # ],
    "GENDER": [
        "gender",
        "sex",
        "male",
        "female",
        "men",
        "women",
    ],
    # "PEEP": [
    #     "PEEP",
    #     "positive end-expiratory pressure",
    #     "positive end expiratory pressure",
    # ],
    # "SAO2": [
    #     "SaO2",
    #     "arterial oxygen saturation",
    #     "oxygen saturation",
    #     "blood oxygen saturation",
    # ],
    # "FIO2": [
    #     "FiO2",
    #     "fraction of inspired oxygen",
    #     "inspired oxygen fraction",
    # ],
    # "DEATH": [
    #     "death",
    #     "mortality",
    #     "in-hospital mortality",
    #     "hospital mortality",
    #     "survival",
    #     "30-day mortality",
    #     "ICU mortality",
    # ],
    # "SPO2": [
    #     "SpO2",
    #     "oxygen saturation",
    #     "pulse oximetry",
    #     "oximetry",
    #     "peripheral oxygen saturation",
    # ],
    # "PAO2": [
    #     "PaO2",
    #     "arterial oxygen tension",
    #     "arterial partial pressure of oxygen",
    #     "partial pressure of oxygen",
    # ],
    # "IHD": [
    #     "ischemic heart disease",
    #     "ischaemic heart disease",
    #     "IHD",
    #     "coronary artery disease",
    #     "CAD",
    #     "myocardial ischemia",
    #     "myocardial ischaemia",
    # ],
    # "COPD": [
    #     "COPD",
    #     "chronic obstructive pulmonary disease",
    #     "chronic obstructive lung disease",
    #     "emphysema",
    #     "chronic bronchitis",
    # ],
    # "MINUTE_VENTILATION": [
    #     "minute ventilation",
    #     "minute volume",
    #     "minute ventilator",
    #     "minute ventilators",
    #     "minute vent volume",
    # ],
    # "SMOKING": [
    #     "smoker",
    #     "smokers",
    #     "smoking",
    #     "tobacco use",
    #     "cigarette",
    #     "cigarettes",
    #     "nicotine",
    #     "tobacco smoking",
    # ],
    # "Trauma": [
    #     "trauma",
    #     "traumatic",
    #     "injury",
    #     "polytrauma",
    #     "multiple trauma",
    #     "trauma patient*",
    # ],
    # "SURGERY": [
    #     "surgery",
    #     "surgical",
    #     "postoperative",
    #     "perioperative",
    #     "post-operative",
    # ],
    # "APSIII": [
    #     "APSIII",
    #     "Acute Physiology Score",
    #     "Acute Physiology Score III",
    #     "APACHE III",
    #     "APACHE",
    # ],
    # "OXYGENATION": [
    #     "oxygenation",
    #     "hypoxemia",
    #     "hypoxaemia",
    #     "hyperoxia",
    #     "oxygen therapy",
    #     "supplemental oxygen",
    #     "oxygen supplementation",
    #     "conservative oxygen therapy",
    #     "liberal oxygen therapy",
    #     "conservative oxygen strategy",
    #     "liberal oxygen strategy",
    #     "conservative oxygenation strategy",
    #     "liberal oxygenation strategy",
    #     "oxygenation strategy",
    # ],
    # "SOFA": [
    #     "SOFA",
    #     "sequential organ failure assessment",
    # ],
    # "PAP": [
    #     "peak airway pressure",
    #     "peak inspiratory pressure",
    #     "PIP",
    #     "peak pressure",
    # ],
    "TIDAL_VOLUME": [
        "tidal volume",
        "VT",
        "lung-protective ventilation",
        "low tidal volume",
    ],
}

MIN_DISTINCT_GROUPS = 2

# ================================================================
# 2) WINDOW-HIT SETTINGS (OPTION A)
# ================================================================
WINDOW_SIZE = 6                 # token distance between group tokens
REQUIRE_WINDOW_HIT = True       # enforce window hit for keeping a sentence

# Option A: require cue between (or near) the two group tokens that are within the window
CUE_SLACK_TOKENS = 2            # allow the cue to be within +/- this many tokens around the pair span

# ---------- Causal / hedge / negation cues ----------
CAUSAL_STRONG = {
    "cause", "causes", "caused", "causal", "causally",
    "lead to", "leads to", "led to",
    "result in", "results in", "resulted in",
    "increase the risk", "increases the risk", "increased the risk",
    "reduce the risk", "reduces the risk", "reduced the risk",
    "effect on", "effects on", "affect", "affects", "affected", "effect of",
    "induce", "induces", "induced", "trigger", "triggers", "triggered",
    "determinant of", "determinants of", "responsible for",
}

CAUSAL_ASSOCIATIVE = {
    "associated with", "association with", "association between",
    "linked to", "link between", "link with",
    "correlated with", "correlation with", "relationship with", "relationship between",
    "predictor of", "risk factor for", "contributor to",
}

NEGATION = {
    "no", "not", "neither", "nor", "lack of", "failed to", "did not", "was not", "were not",
    "does not", "do not", "without", "non-", "absence of", "lack", "unlikely", "didn’t",
    "cannot", "unable to",
}

HEDGES = {
    "may", "might", "could", "possibly", "potentially", "suggest", "suggests", "suggested",
    "appears to", "appeared to", "likely", "unlikely", "tend to", "tended to",
    "evidence of", "evidence suggests", "consistent with",
}

# ---------- Text utils ----------
def normalize(s: str) -> str:
    return re.sub(r"\s+", " ", (s or "").strip())

def sentence_split(text: str) -> List[str]:
    text = normalize(text)
    if not text:
        return []
    if _SPACY_OK:
        return [normalize(s.text) for s in _nlp(text).sents if normalize(s.text)]
    return [normalize(p) for p in re.split(r"(?<=[\.\?\!])\s+", text) if normalize(p)]

def slug_tokens(s: str) -> List[str]:
    """
    Tokenize similarly to your edge-based script:
    - lowercase
    - remove punctuation (keep letters/numbers/hyphen/space)
    - split on whitespace
    """
    s = (s or "").lower()
    s = re.sub(r"[^a-z0-9\s\-]+", " ", s)
    return [t for t in s.split() if t]

def cue_hits(sentence: str):
    s = sentence.lower()
    strong = [c for c in CAUSAL_STRONG if c in s]
    assoc  = [c for c in CAUSAL_ASSOCIATIVE if c in s]
    hedges = [h for h in HEDGES if h in s]
    #negs   = [n for n in NEGATION if n in s]
    negs = []
    for n in NEGATION:
        if " " in n:  # phrase like "did not"
            if n in s:
                negs.append(n)
        else:         # single word like "no"
            if re.search(rf"\b{re.escape(n)}\b", s):
                negs.append(n)

    return strong, assoc, hedges, negs

# ================================================================
# OPTION A HELPERS (cue must be between the two group tokens)
# ================================================================
def find_cue_spans(sentence: str) -> List[Tuple[int, int, str]]:
    """
    Returns list of (start_char, end_char, cue_string) for all strong/assoc cues found.
    Uses literal substring search (like cue_hits).
    """
    s = (sentence or "").lower()
    spans: List[Tuple[int, int, str]] = []
    cues = sorted((CAUSAL_STRONG | CAUSAL_ASSOCIATIVE), key=len, reverse=True)
    for cue in cues:
        start = 0
        while True:
            i = s.find(cue, start)
            if i == -1:
                break
            spans.append((i, i + len(cue), cue))
            start = i + 1
    return spans

def token_char_spans(sentence: str) -> List[Tuple[str, int, int]]:
    """
    Tokenize similar to slug_tokens, but also keep char spans in the original lowercased sentence.
    Pattern keeps hyphenated tokens as one token.
    """
    s = (sentence or "").lower()
    spans: List[Tuple[str, int, int]] = []
    for m in re.finditer(r"[a-z0-9]+(?:-[a-z0-9]+)?", s):
        spans.append((m.group(0), m.start(), m.end()))
    return spans

# ---------- Domain groups + synonyms matching ----------
def build_synonym_regex(term: str) -> re.Pattern:
    """
    - wildcard prefixes: outcome* => \\boutcome\\w*
    - multiword phrases: allow spaces/hyphens between words
    - acronyms/single tokens: word-boundary match
    """
    t = normalize(term)
    if not t:
        return re.compile(r"$^")

    if t.endswith("*") and len(t) > 1:
        base = re.escape(t[:-1])
        return re.compile(rf"\b{base}\w*", flags=re.I)

    if " " in t or "-" in t:
        parts = re.split(r"[\s\-]+", t)
        parts = [re.escape(p) for p in parts if p]
        if not parts:
            return re.compile(r"$^")
        pat = r"\b" + r"(?:[\s\-]+)".join(parts) + r"\b"
        return re.compile(pat, flags=re.I)

    return re.compile(rf"\b{re.escape(t)}\b", flags=re.I)

def compile_domain_regex(domain_groups: Dict[str, List[str]]) -> Dict[str, List[Tuple[str, re.Pattern]]]:
    compiled: Dict[str, List[Tuple[str, re.Pattern]]] = {}
    for group, syns in domain_groups.items():
        if not isinstance(group, str) or not group.strip():
            continue
        group_name = group.strip()

        cleaned: List[str] = []
        seen = set()
        for s in syns:
            if not isinstance(s, str):
                continue
            t = normalize(s)
            if not t:
                continue
            key = t.lower()
            if key in seen:
                continue
            seen.add(key)
            cleaned.append(t)

        compiled[group_name] = [(syn, build_synonym_regex(syn)) for syn in cleaned]
    return compiled

def compile_group_token_sets(domain_groups: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """
    For window-hit: build a token vocabulary for each group from its synonyms.
    Example: "in-hospital mortality" -> ["in-hospital","mortality"] -> tokens after slugging: ["in-hospital","mortality"]
    Wildcards like outcome* contribute the base token "outcome".
    """
    out: Dict[str, set] = {}
    for group, syns in domain_groups.items():
        g = group.strip()
        out[g] = set()
        for syn in syns:
            t = normalize(syn)
            if not t:
                continue
            if t.endswith("*") and len(t) > 1:
                t = t[:-1]
            for tok in slug_tokens(t):
                out[g].add(tok)
    return {g: sorted(list(v)) for g, v in out.items()}

def domain_group_hits(sentence: str, domain_regex: Dict[str, List[Tuple[str, re.Pattern]]]) -> Dict[str, List[str]]:
    hits: Dict[str, List[str]] = {}
    for group, syn_rxs in domain_regex.items():
        matched_syns: List[str] = []
        for syn, rx in syn_rxs:
            if rx.search(sentence):
                matched_syns.append(syn)
        if matched_syns:
            seen = set()
            uniq = []
            for s in matched_syns:
                k = s.lower()
                if k in seen:
                    continue
                seen.add(k)
                uniq.append(s)
            hits[group] = uniq
    return hits

# ---------- Window hit between groups (OPTION A) ----------
def window_hit_for_groups(sentence: str,
                          present_groups: List[str],
                          group_tokens: Dict[str, List[str]],
                          win: int,
                          cue_slack_tokens: int = CUE_SLACK_TOKENS) -> Tuple[bool, Optional[int], str]:
    """
    Option A (refined):
    window_hit=True if there exists a pair of group tokens (from two different groups)
    within <= win tokens AND at least one causal/associative cue occurs between (or near)
    the two token positions (with +/- cue_slack_tokens slack).

    Returns: (window_hit, global_min_token_distance, hit_pairs_string)
    """
    if len(present_groups) < 2:
        return False, None, ""

    toks_spans = token_char_spans(sentence)
    if not toks_spans:
        return False, None, ""

    toks = [t for (t, _, _) in toks_spans]
    cue_spans = find_cue_spans(sentence)
    if not cue_spans:
        return False, None, ""

    # Precompute token indices where each group appears
    group_occ: Dict[str, List[int]] = {}
    for g in present_groups:
        vocab = set(group_tokens.get(g, []))
        idxs = [i for i, tok in enumerate(toks) if tok in vocab]
        group_occ[g] = idxs

    def cue_in_token_range(lo_idx: int, hi_idx: int) -> bool:
        lo = max(0, lo_idx - cue_slack_tokens)
        hi = min(len(toks_spans) - 1, hi_idx + cue_slack_tokens)
        left_char = toks_spans[lo][1]
        right_char = toks_spans[hi][2]
        for cs, ce, _ in cue_spans:
            if cs >= left_char and ce <= right_char:
                return True
        return False

    best_dist: Optional[int] = None
    hit_pairs: List[str] = []

    for i in range(len(present_groups)):
        for j in range(i + 1, len(present_groups)):
            g1 = present_groups[i]
            g2 = present_groups[j]
            occ1 = group_occ.get(g1, [])
            occ2 = group_occ.get(g2, [])
            if not occ1 or not occ2:
                continue

            for idx1 in occ1:
                for idx2 in occ2:
                    d = abs(idx1 - idx2)
                    if best_dist is None or d < best_dist:
                        best_dist = d
                    if d <= win:
                        lo_idx, hi_idx = (idx1, idx2) if idx1 <= idx2 else (idx2, idx1)
                        if cue_in_token_range(lo_idx, hi_idx):
                            hit_pairs.append(f"{g1}<->{g2}(d={d})")

    return (len(hit_pairs) > 0), best_dist, "; ".join(sorted(set(hit_pairs)))

def sentence_is_causal(sent: str,
                       domain_regex: Dict[str, List[Tuple[str, re.Pattern]]],
                       group_tokens: Dict[str, List[str]],
                       min_groups: int,
                       win: int,
                       require_window: bool) -> Tuple[bool, Dict]:
    """
    Keep sentence if:
      - it contains synonyms from >= min_groups DISTINCT domain groups, AND
      - it contains >=1 causal/associative cue, AND
      - (optionally) window_hit is True for at least one pair of present groups within win tokens
        AND cue-between condition (Option A).
    Negation is NOT dropped; only flagged.
    """
    strong, assoc, hedges, negs = cue_hits(sent)
    g_hits = domain_group_hits(sent, domain_regex)
    group_count = len(g_hits)
    present_groups = sorted(g_hits.keys(), key=lambda x: x.lower())

    has_cue = bool(strong) or bool(assoc)

    win_hit, min_dist, hit_pairs = window_hit_for_groups(sent, present_groups, group_tokens, win)

    keep = (group_count >= min_groups) and has_cue
    if require_window:
        keep = keep and win_hit

    group_syn_str = "; ".join(
        f"{g}:[{', '.join(syns)}]" for g, syns in sorted(g_hits.items(), key=lambda kv: kv[0].lower())
    )

    meta = {
        "domain_groups": "; ".join(present_groups),
        "domain_group_count": group_count,
        "domain_group_synonyms": group_syn_str,
        "strong_cues": "; ".join(strong),
        "assoc_cues": "; ".join(assoc),
        "hedge_or_neg": "; ".join(hedges + negs),
        "negated": bool(negs),

        # window-hit outputs
        "window_size": win,
        "window_hit": win_hit,
        "min_token_distance": "" if min_dist is None else str(min_dist),
        "window_pairs": hit_pairs,
    }
    return keep, meta

# ---------- Study design classification ----------
PUBTYPE_TO_BUCKET = {
    "Randomized Controlled Trial": "RCT",
    "Clinical Trial, Phase I": "Other",
    "Clinical Trial, Phase II": "Other",
    "Clinical Trial, Phase III": "Other",
    "Clinical Trial, Phase IV": "Other",
    "Controlled Clinical Trial": "Other",
    "Pragmatic Clinical Trial": "Other",
    "Clinical Trial": "Other",
    "Cohort Studies": "Cohort",
    "Case-Control Studies": "Case-Control",
    "Cross-Sectional Studies": "Other",
    "Observational Study": "Cohort",
    "Prospective Studies": "Cohort",
    "Retrospective Studies": "Cohort",
    "Case Reports": "Case Report",
    "Case Series": "Case Series",
    "Systematic Review": "Systematic Review",
    "Meta-Analysis": "Meta-analysis",
    "Review": "Other",
    "Umbrella Review": "Systematic Review",
    "Comparative Study": "Other",
    "Evaluation Study": "Other",
    "Multicenter Study": "Other",
    "Letter": "Other",
    "Editorial": "Other",
    "Critically Appraised Topic": "CAT",
    "Critically Appraised Individual Article": "CAIA",
}

FALLBACK_REGEX = [
    (re.compile(r"\brandomi[sz]ed\b|\brandomly assigned\b|\bdouble-?blind\b|\bsingle-?blind\b"), "RCT"),
    (re.compile(r"\b(cohort|prospective|retrospective)\b"), "Cohort"),
    (re.compile(r"\bcase[- ]control\b"), "Case-Control"),
    (re.compile(r"\bcross[- ]sectional\b"), "Other"),
    (re.compile(r"\bcase (series|report)\b"), "Case Report"),
    (re.compile(r"\bsystematic review\b"), "Systematic Review"),
    (re.compile(r"\bmeta[- ]analysis\b"), "Meta-analysis"),
    (re.compile(r"\bobservational\b"), "Cohort"),
    (re.compile(r"\bpragmatic trial\b"), "Other"),
    (re.compile(r"\bpilot study\b"), "Other"),
    (re.compile(r"\bmendelian randomization\b"), "Cohort"),
    (re.compile(r"\bcritically[-\s]+appraised\s+topic(s)?\b", re.I), "CAT"),
    (re.compile(r"\bcritically[-\s]+appraised\s+individual\s+article(s)?\b", re.I), "CAIA"),
    (re.compile(r"\bcritically[-\s]+appraised\s+article(s)?\b", re.I), "CAIA"),
]

def classify_study_design(pub_types: List[str], title: str, abstract: str):
    signals = []
    buckets = []
    for pt in pub_types:
        pt_norm = pt.strip()
        if pt_norm in PUBTYPE_TO_BUCKET:
            buckets.append(PUBTYPE_TO_BUCKET[pt_norm])
            signals.append(f"PubType:{pt_norm}")
    if buckets:
        priority = [
            "Meta-analysis", "Systematic Review", "CAT", "CAIA", "RCT", 
            "Cohort", "Case-Control", "Case Series", "Case Report", "Other"
        ]
        for p in priority:
            if p in buckets:
                return p, 0.98 if p in {"RCT", "Meta-analysis", "Systematic Review"} else 0.9, signals

    text = f"{title} {abstract}".lower()
    for rx, label in FALLBACK_REGEX:
        if rx.search(text):
            signals.append(f"Regex:{label}")
            conf = 0.95 if label in {"RCT", "Meta-analysis", "Systematic Review"} else 0.75
            return label, conf, signals

    return "Other", 0.4, signals

# ---------- PubMed helpers ----------
def pubmed_search(query: str, retmax: int = 200, api_key: str = "") -> List[str]:
    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "sort": "most+recent",
        "retmax": retmax,
    }
    if api_key:
        params["api_key"] = api_key
    url = f"{NCBI_BASE}/esearch.fcgi?" + urlencode(params)
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.json().get("esearchresult", {}).get("idlist", [])

def pubmed_fetch(pmids: List[str], api_key: str = "") -> Dict[str, Dict]:
    out: Dict[str, Dict] = {}
    if not pmids:
        return out

    CHUNK = 200
    for i in range(0, len(pmids), CHUNK):
        chunk = pmids[i:i + CHUNK]
        params = {"db": "pubmed", "id": ",".join(chunk), "retmode": "xml"}
        if api_key:
            params["api_key"] = api_key

        url = f"{NCBI_BASE}/efetch.fcgi"
        r = requests.post(url, data=params, timeout=60)
        r.raise_for_status()
        xml = r.text

        articles = re.split(r"<PubmedArticle>", xml)[1:]
        for art in articles:
            pmid_m = re.search(r"<PMID.*?>(\d+)</PMID>", art)
            if not pmid_m:
                continue
            pmid = pmid_m.group(1)

            title = html.unescape(" ".join(re.findall(r"<ArticleTitle>(.*?)</ArticleTitle>", art, flags=re.S)))
            abs_chunks = re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>", art, flags=re.S)
            abstract = html.unescape(" ".join(abs_chunks))

            pub_types = re.findall(r"<PublicationType[^>]*>(.*?)</PublicationType>", art)
            pub_types = [html.unescape(normalize(pt)) for pt in pub_types]

            journal_title = ""
            m = re.search(r"<Journal>.*?<Title>(.*?)</Title>.*?</Journal>", art, flags=re.S)
            if m:
                journal_title = html.unescape(normalize(m.group(1)))
            journal_iso = ""
            m = re.search(r"<Journal>.*?<ISOAbbreviation>(.*?)</ISOAbbreviation>.*?</Journal>", art, flags=re.S)
            if m:
                journal_iso = html.unescape(normalize(m.group(1)))

            doi = ""
            m = re.search(r"<ELocationID[^>]*EIdType=\"doi\"[^>]*>(.*?)</ELocationID>", art, flags=re.S | re.I)
            if m:
                doi = normalize(html.unescape(m.group(1)))
            if not doi:
                m = re.search(r"<ArticleId IdType=\"doi\">(.*?)</ArticleId>", art, flags=re.S | re.I)
                if m:
                    doi = normalize(html.unescape(m.group(1)))

            paper_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            if doi:
                paper_url = f"https://doi.org/{doi}"

            out[pmid] = {
                "title": normalize(title),
                "abstract": normalize(abstract),
                "pub_types": pub_types,
                "journal_title": journal_title,
                "journal_iso": journal_iso,
                "doi": doi,
                "paper_url": paper_url,
            }

        time.sleep(0.34)

    return out

# ---------- Sample size extraction ----------
_SAMPLE_PATTERNS = [
    re.compile(r"\b[Nn]\s*=\s*(\d{2,5})\b"),
    re.compile(r"\b(?:participants|patients|subjects|adults|children)\s*\(\s*[Nn]\s*=\s*(\d{2,5})\s*\)"),
    re.compile(r"\b(?:enrolled|included|analyzed|randomi[sz]ed|recruited)\s+(\d{2,5})\s+(?:participants|patients|subjects)\b"),
    re.compile(r"\b(\d{2,5})\s+(?:participants|patients|subjects)\b"),
]

def extract_sample_size(title: str, abstract: str) -> Optional[int]:
    txt = f"{title} {abstract}"
    hits: List[int] = []
    for rx in _SAMPLE_PATTERNS:
        for m in rx.finditer(txt):
            try:
                hits.append(int(m.group(1)))
            except Exception:
                pass
    return max(hits) if hits else None

def read_query_file(path: str) -> str:
    with open(path, "r", encoding="utf-8") as f:
        return f.read().strip()

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results.csv")
    ap.add_argument("--query", required=False, default="", help="PubMed query string OR a path to a .txt file")
    ap.add_argument("--query-file", default="", help="Path to a .txt file containing the PubMed query")
    ap.add_argument("--pmids", default="", help="Comma-separated PMIDs (optional; overrides query search)")
    ap.add_argument("--max-articles", type=int, default=200)
    ap.add_argument("--api-key", default="")

    args = ap.parse_args()

    # ---- Resolve query ----
    query = args.query.strip()

    if args.query_file:
        if not os.path.exists(args.query_file):
            print(f"ERROR: --query-file not found: {args.query_file}", file=sys.stderr)
            sys.exit(1)
        query = read_query_file(args.query_file)

    # Backward-friendly: if user passed --query query.txt and file exists, read it
    if query.lower().endswith(".txt") and os.path.exists(query):
        query = read_query_file(query)

    if not query and not args.pmids.strip():
        print("ERROR: Provide --query \"...\" or --query-file query.txt (or use --pmids).", file=sys.stderr)
        sys.exit(1)

    if MIN_DISTINCT_GROUPS < 1:
        print("ERROR: MIN_DISTINCT_GROUPS must be >= 1.", file=sys.stderr)
        sys.exit(1)
    if len(DOMAIN_GROUPS) < MIN_DISTINCT_GROUPS:
        print(
            f"ERROR: You defined {len(DOMAIN_GROUPS)} DOMAIN_GROUPS but MIN_DISTINCT_GROUPS={MIN_DISTINCT_GROUPS}.",
            file=sys.stderr
        )
        sys.exit(1)

    domain_regex = compile_domain_regex(DOMAIN_GROUPS)
    group_tokens = compile_group_token_sets(DOMAIN_GROUPS)

    # PMIDs via query search or explicit list
    if args.pmids.strip():
        pmids = [p.strip() for p in args.pmids.split(",") if p.strip().isdigit()]
        print(f"Using {len(pmids)} provided PMIDs.", file=sys.stderr)
    else:
        pmids = pubmed_search(query, retmax=args.max_articles, api_key=args.api_key)
        print(f"Found {len(pmids)} PubMed IDs.", file=sys.stderr)

    meta = pubmed_fetch(pmids, api_key=args.api_key)
    print(f"Fetched metadata for {len(meta)} PMIDs.", file=sys.stderr)

    rows: List[Dict] = []
    for pmid in pmids:
        m = meta.get(pmid)
        if not m:
            continue

        title = m.get("title", "")
        abstract = m.get("abstract", "")
        pub_types = m.get("pub_types", [])
        journal_title = m.get("journal_title", "")
        journal_iso = m.get("journal_iso", "")
        paper_url = m.get("paper_url", "")

        design_label, design_conf, design_signals = classify_study_design(pub_types, title, abstract)
        sample_n = extract_sample_size(title, abstract)

        for sent in sentence_split(abstract):
            keep, meta_bits = sentence_is_causal(
                sent=sent,
                domain_regex=domain_regex,
                group_tokens=group_tokens,
                min_groups=MIN_DISTINCT_GROUPS,
                win=WINDOW_SIZE,
                require_window=REQUIRE_WINDOW_HIT,
            )
            if not keep:
                continue

            rows.append({
                "edge": meta_bits["domain_groups"],
                "pmid": pmid,
                #"venue": journal_title or journal_iso,
                "study_design": design_label,
                #"design_signals": "; ".join(design_signals),
                "sample_size": sample_n if sample_n is not None else "",
                #"domain_group_synonyms": meta_bits["domain_group_synonyms"],
                "strong_cues": meta_bits["strong_cues"],
                "assoc_cues": meta_bits["assoc_cues"],
                "hedge_or_neg": meta_bits["hedge_or_neg"],
                #"negated": meta_bits["negated"],
                "sentence": sent,
                "url": paper_url,
            })

    rows.sort(key=lambda r: (r["study_design"], r["pmid"]))

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            # "edge", "pmid", "venue",
            # "study_design", "design_signals",
            # "sample_size", "url",
            # "domain_group_synonyms",
            # "strong_cues", "assoc_cues", "hedge_or_neg", "negated",
            # "sentence",
            "edge", "pmid",
            "study_design",
            "sample_size", 
            "strong_cues", "assoc_cues", "hedge_or_neg",
            "sentence", "url",
        ])
        for r in rows:
            w.writerow([
                # r["edge"], r["pmid"], r["venue"],
                # r["study_design"], r["design_signals"],
                # r["sample_size"], r["url"],
                # r["domain_group_synonyms"],
                # r["strong_cues"], r["assoc_cues"], r["hedge_or_neg"], r["negated"],
                # r["sentence"],
                r["edge"], r["pmid"],
                r["study_design"],
                r["sample_size"], 
                r["strong_cues"], r["assoc_cues"], r["hedge_or_neg"],
                r["sentence"], r["url"],
            ])

    print(f"Wrote {len(rows)} rows to {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
