r"""Two independent, deterministic claim extractors and their agreement with the AI labels.

    python -m track2.engine_extractors [--out out]

(A) `rule_extract(sentence, drug)`: documented lexicon, negation and scope handling. Returns direction in
    {benefit, harm, no_effect, unclear}, the system, and the study-design cue.
(B) `LearnedExtractor`: TF-IDF (1-2 grams) + L2-regularised logistic regression (scikit-learn if installed, otherwise a pure
    NumPy fallback), evaluated out-of-fold with folds GROUPED BY PMID (a paper is never in training and test), seed 20261015.
(C) `consensus(sentence, drug, ai_label=None)`: majority of rule, learned and (optional) AI label; 'unresolved' if no majority.

Mapping of the AI labels (declared before the results): may_help -> benefit; may_harm -> harm; no_effect -> no_effect;
unclear -> unclear; selectively_kills_aneuploid_cells -> harm (it kills cells; the tumour-context benefit is a use
decision, not a direction of effect on the cell). The AI `system` field is free text in mixed Spanish/English and is not
scored against the system lexicon; only the distribution of the rule-based system is reported.

The lexicon below was written once, before looking at any discordant sentence; the agreement reported is that of this
first version. Nothing was tuned afterwards (no post hoc version was made).
Writes `out/engine/extractor_agreement.json` (sorted keys). No network, no patient data.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import numpy as np

SEED = 20261015
LABELS = ["benefit", "harm", "no_effect", "unclear"]
AI_MAP = {"may_help": "benefit", "may_harm": "harm", "no_effect": "no_effect", "unclear": "unclear", "selectively_kills_aneuploid_cells": "harm"}

# ----------------------------------------------------------------------------------------------------------------
# (A) Lexicon of the rule-based extractor
# ----------------------------------------------------------------------------------------------------------------
# Verbs that are good for the patient whatever the object.
BENEFIT_CUES = r"improv\w*|rescu\w*|ameliorat\w*|restor\w*|alleviat\w*|attenuat\w*|prevent\w*|protect\w*|beneficial|benefit\w*|rejuvenat\w*|revers(?:e|ed|es|al)|normali[sz]\w*|effective|extends?|extended|lifespan extension|overwhelmingly positive|positive effect"
# Words that are bad for the patient whatever the object.
HARM_CUES = (r"toxic\w*|cytotoxic\w*|adverse|side[- ]effects?|deleterious|detrimental|harm\w*|impair\w*|worsen\w*|exacerbat\w*|atrophy|retardation|lethal|"
             r"immunosuppress\w*|hyperlipid\w*|hyperglyc\w*|kill(?:s|ed|ing)?|cell death|negative(?:ly)?|sensitiv\w*|hypersensitiv\w*|favou?red placebo|resistan(?:t|ce)|"
             r"tumou?rigen\w*|fatal|death|mortality increase")
# Verbs of decrease / increase: the polarity depends on the object (a disease measure or a health measure).
DOWN_VERBS = r"reduc\w*|decreas\w*|lower\w*|diminish\w*|inhibit\w*|suppress\w*|declin\w*|abolish\w*|ablat\w*|loss of"
UP_VERBS = r"increas\w*|elevat\w*|rais\w*|induc\w*|enhanc\w*|promot\w*|accumulat\w*|higher|greater"
DISEASE_MEASURES = r"tumou?rs?|cancers?|senescen\w*|apoptos\w*|aggregat\w*|inflammat\w*|oxidative|damage|mortality|lesions?|symptoms?|progression|proliferation|aneuploid\w*|toxicity|stress|ros\b|death|fibrosis|pathology|defects?|instability|neuropath\w*|seizures?|hyperactivation|hyperactivity|microcephaly|infection"
HEALTH_MEASURES = r"growth|height|weight|survival|lifespan|function|muscle|testosterone|wound healing|bone|density|viability|differentiation|fertility|mass|size|autophagy|mitophagy|clearance|strength|performance|bmi"
# Absence of an effect.
NO_EFFECT_CUES = (r"no (?:significant |statistically significant |measurable |clear |obvious )?(?:difference|differences|effect|effects|change|changes|impact|improvement|association|benefit)|"
                  r"did not|does not|do not|is not|was not|were not|not (?:significantly |substantially )?(?:affect\w*|alter\w*|chang\w*|differ\w*|requir\w*|improv\w*|reduc\w*|increas\w*)|"
                  r"failed to|fail to|unaffected|unchanged|without (?:significant |any )?(?:effect|change|benefit)|not necessary|unlikely to|neither|nor\b|no evidence|not associated")
NEGATORS = r"no|not|n't|without|never|neither|nor|failed to|fail to|unable to|lack of|absence of|cannot|did not|does not|do not"
# Hedges and statements of uncertainty.
UNCLEAR_STRONG = r"remains? (?:lacking|unclear|unknown|uncertain|to be)|unclear|unknown|uncertain|inconclusive|insufficient|limited evidence|conflicting|mixed|whether|further (?:studies|research|work)|not (?:yet )?(?:known|clear|established)"
HEDGES = r"may|might|could|suggest\w*|potential\w*|possibl\w*|appear\w*|seem\w*|likely"

# System lexicon (checked in this order; first match wins). review is checked first as a statement-level cue.
SYSTEM_LEXICON: list[tuple[str, str]] = [
    ("review or unspecified", r"\b(?:review|meta-analys\w*|narrative|overview|literature)\b"),
    ("invertebrate or yeast", r"\b(?:drosophila|fruit fl\w+|c\. elegans|caenorhabditis|nematode|worms?|yeast|saccharomyces|larva\w*|insect\w*)\b"),
    ("other vertebrate", r"\b(?:zebrafish|danio|xenopus|chick\w*|pigs?|porcine|dogs?|canine|rabbits?|monkeys?|primates?|sheep|bovine|cattle|fish)\b"),
    ("mouse or rat", r"\b(?:mice|mouse|murine|rats?|mefs?|rodents?|knockout|c57bl/6)\b"),
    ("in vitro biochemical", r"\b(?:purified|recombinant|cell-free|kinase assay|enzym\w+ assay|biochemical|in vitro kinase|binding assay)\b"),
    ("human patients", r"\b(?:patients?|children|child|boys?|girls?|infants?|individuals|participants|volunteers|adults|men|women|adolescents?|cohort|trial|case|subjects|toddlers?|neonates?|humans)\b"),
    ("human cells", r"\b(?:hela|rpe-?1|hek293\w*|fibroblasts?|cell lines?|cells|a375|a549|myoblasts?|keratinocytes?|lymphocytes|organoids?|cultured)\b"),
]
DESIGN_LEXICON: list[tuple[str, str]] = [
    ("randomised trial", r"\b(?:randomi[sz]ed|placebo|rct|double-blind|phase (?:i{1,3}|[1-3]) (?:trial|study))\b"),
    ("case report", r"\b(?:case report|case study|-year-old|-month-old|a (?:boy|girl|patient|child|woman|man) (?:with|who|presented))\b"),
    ("cohort", r"\b(?:cohort|retrospective|prospective|registry|observational|case-control)\b"),
    ("review", r"\b(?:review|meta-analys\w*|narrative|overview)\b"),
    ("in vitro", r"\b(?:in vitro|cell lines?|cultured|cultures?|fibroblasts?|clones?|hela|rpe-?1)\b"),
]
DRUG_ALIASES = {"SIROLIMUS": ["rapamycin", "mtor inhibitor", "rapalog"], "PHENYLBUTANOIC ACID": ["phenylbutyrate", "4-pba", "pba"], "NIACINAMIDE": ["nicotinamide", "nam"],
                "SELINEXOR": ["kpt-330"], "BOSUTINIB": ["ski606", "ski-606"], "TRAMETINIB": ["mek inhibitor"], "EVEROLIMUS": ["rad001"]}

_CLAUSE_SPLIT = re.compile(r";|,? (?:but|whereas|while|however|although|though)\b,?")


def _clauses(sentence: str, drug: str) -> list[str]:
    """Scope: split into clauses and keep those that mention the drug (or an alias); all clauses if none does."""
    parts = [p.strip() for p in _CLAUSE_SPLIT.split(sentence.lower()) if p and p.strip()]
    names = [drug.lower()] + DRUG_ALIASES.get(drug.upper(), [])
    names = [n for n in names if n]
    keep = [p for p in parts if any(re.search(r"\b" + re.escape(n) + r"\b", p) for n in names)]
    return keep or parts or [sentence.lower()]


def _negated(clause: str, start: int) -> bool:
    """A negator in the 4 words before position `start` (but not 'not only')."""
    before = clause[:start].split()[-4:]
    text = " ".join(before)
    if re.search(r"not only", text):
        return False
    return bool(re.search(r"(?:^|\s)(?:" + NEGATORS + r")(?:\s|$)", text))


def _object_class(clause: str, end: int) -> str | None:
    window = " ".join(clause[end:].split()[:5])
    d, h = re.search(DISEASE_MEASURES, window), re.search(HEALTH_MEASURES, window)
    if d and (not h or d.start() <= h.start()):
        return "disease"
    if h:
        return "health"
    return None


def rule_direction(sentence: str, drug: str) -> str:
    b = h = n = u = 0.0
    for clause in _clauses(sentence, drug):
        hedge = 0.7 if re.search(r"\b(?:" + HEDGES + r")\b", clause) else 1.0
        for m in re.finditer(r"\b(?:" + NO_EFFECT_CUES + r")", clause):
            n += 1.0
        for m in re.finditer(r"\b(?:" + UNCLEAR_STRONG + r")\b", clause):
            u += 1.5
        for pol, pattern in (("b", BENEFIT_CUES), ("h", HARM_CUES)):
            for m in re.finditer(r"\b(?:" + pattern + r")\b", clause):
                if _negated(clause, m.start()):
                    n += 1.0
                elif pol == "b":
                    b += hedge
                else:
                    h += hedge
        for verbs, down in ((DOWN_VERBS, True), (UP_VERBS, False)):
            for m in re.finditer(r"\b(?:" + verbs + r")\b", clause):
                obj = _object_class(clause, m.end())
                if obj is None:
                    continue
                good = (obj == "disease") == down  # decreasing a disease measure or increasing a health measure
                if _negated(clause, m.start()):
                    n += 1.0
                elif good:
                    b += hedge
                else:
                    h += hedge
    top = max(b, h, n, u)
    if top < 0.5:
        return "unclear"
    if n >= max(b, h, u) and n > 0:
        return "no_effect"
    if u > max(b, h, n):
        return "unclear"
    if b == h:
        return "unclear"
    return "benefit" if b > h else "harm"


def rule_system(sentence: str, context: str = "") -> str:
    text = (sentence + " " + context).lower()
    for name, pattern in SYSTEM_LEXICON:
        if re.search(pattern, text):
            return name
    return "review or unspecified"


def rule_design(sentence: str) -> str | None:
    text = sentence.lower()
    for name, pattern in DESIGN_LEXICON:
        if re.search(pattern, text):
            return name
    return None


def rule_extract(sentence: str, drug: str) -> dict[str, Any]:
    return {"direction": rule_direction(sentence, drug), "system": rule_system(sentence), "design": rule_design(sentence)}


# ----------------------------------------------------------------------------------------------------------------
# (B) Learned extractor
# ----------------------------------------------------------------------------------------------------------------
try:  # pragma: no cover - depends on the environment
    from sklearn.feature_extraction.text import TfidfVectorizer
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedGroupKFold
    HAVE_SKLEARN = True
except Exception:  # pragma: no cover
    HAVE_SKLEARN = False


class _NumpyTfidfLR:
    """Fallback: unigram+bigram TF-IDF (min_df=2) and multinomial L2 logistic regression by full-batch gradient descent."""

    def __init__(self, C: float = 1.0, iters: int = 400, lr: float = 0.5):
        self.C, self.iters, self.lr = C, iters, lr

    @staticmethod
    def _grams(text: str) -> list[str]:
        w = re.findall(r"[a-z0-9\-]+", text.lower())
        return w + [a + " " + b for a, b in zip(w, w[1:])]

    def _tfidf(self, texts: list[str]) -> np.ndarray:
        X = np.zeros((len(texts), len(self.vocab)))
        for i, t in enumerate(texts):
            for g in self._grams(t):
                j = self.vocab.get(g)
                if j is not None:
                    X[i, j] += 1
        X = (1 + np.log(np.maximum(X, 1))) * (X > 0) * self.idf
        return X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)

    def fit(self, texts: list[str], y: list[str]) -> "_NumpyTfidfLR":
        df: dict[str, int] = {}
        for t in texts:
            for g in set(self._grams(t)):
                df[g] = df.get(g, 0) + 1
        self.vocab = {g: i for i, g in enumerate(sorted(g for g, c in df.items() if c >= 2))}
        self.idf = np.array([np.log((1 + len(texts)) / (1 + df[g])) + 1 for g in sorted(self.vocab, key=self.vocab.get)])
        self.classes = sorted(set(y))
        X = self._tfidf(texts)
        Y = np.array([[c == v for c in self.classes] for v in y], dtype=float)
        W, b, n = np.zeros((X.shape[1], len(self.classes))), np.zeros(len(self.classes)), len(texts)
        for _ in range(self.iters):
            Z = X @ W + b
            P = np.exp(Z - Z.max(1, keepdims=True))
            P /= P.sum(1, keepdims=True)
            G = P - Y
            W -= self.lr * (X.T @ G + W / self.C / 1.0) / n * 1.0
            b -= self.lr * G.mean(0)
        self.W, self.b = W, b
        return self

    def predict(self, texts: list[str]) -> list[str]:
        Z = self._tfidf(texts) @ self.W + self.b
        return [self.classes[i] for i in Z.argmax(1)]


class LearnedExtractor:
    """TF-IDF + regularised logistic regression on the sentence text."""

    def __init__(self, C: float = 10.0):
        self.C = C

    def fit(self, texts: list[str], y: list[str]) -> "LearnedExtractor":
        if HAVE_SKLEARN:
            self.vec = TfidfVectorizer(ngram_range=(1, 2), min_df=2, sublinear_tf=True, lowercase=True)
            X = self.vec.fit_transform(texts)
            self.clf = LogisticRegression(C=self.C, max_iter=2000, random_state=SEED)
            self.clf.fit(X, y)
        else:  # pragma: no cover
            self.clf = _NumpyTfidfLR(C=self.C).fit(texts, y)
        return self

    def predict(self, texts: list[str]) -> list[str]:
        if HAVE_SKLEARN:
            return [str(v) for v in self.clf.predict(self.vec.transform(texts))]
        return self.clf.predict(texts)  # pragma: no cover


def grouped_folds(y: list[str], groups: list[str], k: int = 5, seed: int = SEED) -> list[tuple[np.ndarray, np.ndarray]]:
    """Folds in which no group (PMID) is on both sides; stratified by label when scikit-learn is available."""
    idx = np.arange(len(y))
    if HAVE_SKLEARN:
        return [(a, b) for a, b in StratifiedGroupKFold(n_splits=k, shuffle=True, random_state=seed).split(idx, y, groups)]
    rng = np.random.default_rng(seed)  # pragma: no cover
    uniq = sorted(set(groups))
    rng.shuffle(uniq)
    fold_of = {g: i % k for i, g in enumerate(uniq)}
    f = np.array([fold_of[g] for g in groups])
    return [(idx[f != i], idx[f == i]) for i in range(k)]


def out_of_fold(texts: list[str], y: list[str], groups: list[str], k: int = 5) -> list[str]:
    pred = [""] * len(y)
    for tr, te in grouped_folds(y, groups, k):
        assert not ({groups[i] for i in tr} & {groups[i] for i in te}), "same PMID in training and test"
        m = LearnedExtractor().fit([texts[i] for i in tr], [y[i] for i in tr])
        for i, p in zip(te, m.predict([texts[i] for i in te])):
            pred[i] = p
    return pred


# ----------------------------------------------------------------------------------------------------------------
# (C) Consensus
# ----------------------------------------------------------------------------------------------------------------
_DEFAULT_MODEL: LearnedExtractor | None = None


def _default_model() -> LearnedExtractor | None:
    """Trained once on all labelled sentences (only for consensus() calls without an explicit model)."""
    global _DEFAULT_MODEL
    if _DEFAULT_MODEL is None:
        rows = load_entries(Path("out"))
        if len(rows) < 10:
            return None
        _DEFAULT_MODEL = LearnedExtractor().fit([r["sentence"] for r in rows], [r["ai"] for r in rows])
    return _DEFAULT_MODEL


def consensus(sentence: str, drug: str, ai_label: str | None = None, model: LearnedExtractor | None = None) -> dict[str, Any]:
    """Majority of rule, learned and optional AI label (AI label may be in the AI vocabulary or in the 4-class vocabulary).
    `label` is 'unresolved' when no label has a majority (all three differ; with two sources only, when they differ)."""
    votes = {"rule": rule_direction(sentence, drug)}
    model = model or _default_model()
    if model is not None:
        votes["learned"] = model.predict([sentence])[0]
    if ai_label is not None:
        votes["ai"] = AI_MAP.get(ai_label, ai_label)
    counts: dict[str, int] = {}
    for v in votes.values():
        counts[v] = counts.get(v, 0) + 1
    top = max(counts.values())
    winners = sorted(k for k, c in counts.items() if c == top)
    label = winners[0] if len(winners) == 1 and top >= 2 else "unresolved"
    return {"label": label, "votes": votes}


# ----------------------------------------------------------------------------------------------------------------
# Agreement statistics
# ----------------------------------------------------------------------------------------------------------------
def kappa(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) == 0:
        return float("nan")
    po = float(np.mean(a == b))
    pe = sum(float(np.mean(a == c)) * float(np.mean(b == c)) for c in set(a.tolist()) | set(b.tolist()))
    return 1.0 if pe >= 1 else (po - pe) / (1 - pe)


def _metrics(ref: np.ndarray, hyp: np.ndarray) -> dict[str, float]:
    out = {"accuracy": float(np.mean(ref == hyp)), "kappa": kappa(ref, hyp)}
    for c in LABELS:
        r, h = ref == c, hyp == c
        out[f"recall_{c}"] = float(np.mean(h[r])) if r.any() else float("nan")
        out[f"kappa_{c}"] = kappa(r.astype(int), h.astype(int)) if (r.any() or h.any()) else float("nan")
    return out


def agreement(ref: list[str], hyp: list[str], groups: list[str], n_boot: int = 2000, seed: int = SEED) -> dict[str, Any]:
    """Point estimates and bootstrap 95 % intervals resampling PMIDs (papers) with replacement."""
    ref_a, hyp_a, g = np.array(ref), np.array(hyp), np.array(groups)
    point = _metrics(ref_a, hyp_a)
    papers = sorted(set(groups))
    members = {p: np.where(g == p)[0] for p in papers}
    rng = np.random.default_rng(seed)
    boots: dict[str, list[float]] = {k: [] for k in point}
    for _ in range(n_boot):
        pick = rng.integers(0, len(papers), len(papers))
        idx = np.concatenate([members[papers[i]] for i in pick])
        m = _metrics(ref_a[idx], hyp_a[idx])
        for k, v in m.items():
            boots[k].append(v)
    res: dict[str, Any] = {}
    for k, v in point.items():
        arr = np.array(boots[k], dtype=float)
        arr = arr[~np.isnan(arr)]
        lo, hi = (float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))) if len(arr) else (float("nan"), float("nan"))
        res[k] = {"estimate": round(v, 4) if v == v else None, "ci95": [round(lo, 4) if lo == lo else None, round(hi, 4) if hi == hi else None]}
    return res


def load_entries(out: Path) -> list[dict[str, Any]]:
    """The labelled sentences, read exactly as `track2.literature_summary.collect` reads them (verified entries only)."""
    try:
        from track2.literature_summary import collect
        by_drug = collect(out / "literature_hypotheses")
    except Exception:
        return []
    rows = []
    for drug in sorted(by_drug):
        for r in by_drug[drug]["entries"]:
            if r.get("sentence") and r.get("direction") in AI_MAP:
                rows.append({"drug": drug, "pmid": str(r["pmid"]) if r.get("pmid") else f"nopmid-{len(rows)}", "sentence": r["sentence"], "ai_raw": r["direction"], "ai": AI_MAP[r["direction"]],
                             "ai_system": r.get("system"), "ai_level": r.get("level")})
    return rows


def run(out: Path) -> dict[str, Any]:
    rows = load_entries(out)
    texts, groups = [r["sentence"] for r in rows], [r["pmid"] for r in rows]
    ai = [r["ai"] for r in rows]
    rule = [rule_direction(r["sentence"], r["drug"]) for r in rows]
    learned = out_of_fold(texts, ai, groups)
    pairs = {"rule_vs_ai": (ai, rule), "learned_oof_vs_ai": (ai, learned), "rule_vs_learned_oof": (rule, learned)}
    result: dict[str, Any] = {
        "n_sentences": len(rows), "n_papers": len(set(groups)), "seed": SEED, "bootstrap_resamples": 2000, "folds": "StratifiedGroupKFold(5) grouped by PMID" if HAVE_SKLEARN else "numpy grouped 5-fold by PMID",
        "label_mapping": AI_MAP, "lexicon_version": "v1, written before looking at any discordant sentence; not tuned afterwards",
        "label_counts": {"ai": {c: ai.count(c) for c in LABELS}, "rule": {c: rule.count(c) for c in LABELS}, "learned_oof": {c: learned.count(c) for c in LABELS}},
        "agreement": {name: agreement(r, h, groups) for name, (r, h) in pairs.items()},
        "rule_system_distribution": {},
        "rule_design_distribution": {},
    }
    systems = [rule_system(r["sentence"]) for r in rows]
    designs = [rule_design(r["sentence"]) or "none" for r in rows]
    result["rule_system_distribution"] = {s: systems.count(s) for s in sorted(set(systems))}
    result["rule_design_distribution"] = {s: designs.count(s) for s in sorted(set(designs))}
    scored = []
    for r, a, ru, le in zip(rows, ai, rule, learned):
        distinct = len({a, ru, le})
        ai_odd = int(distinct == 2 and ru == le)
        scored.append((-(distinct * 10 + ai_odd), r["pmid"], r["drug"], r["sentence"], {"ai": r["ai_raw"], "ai_4class": a, "rule": ru, "learned_oof": le}, distinct))
    scored.sort(key=lambda t: (t[0], t[1], t[2], t[3]))
    result["most_discordant_20"] = [{"drug": d, "pmid": p, "sentence_120": s[:120], "labels": lab, "distinct_labels": nd} for _, p, d, s, lab, nd in scored[:20]]
    result["n_all_three_differ"] = sum(1 for t in scored if t[5] == 3)
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out")
    args = ap.parse_args()
    out = Path(args.out)
    result = run(out)
    dest = out / "engine"
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "extractor_agreement.json").write_text(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    for name, a in result["agreement"].items():
        print(name, {k: (v["estimate"], v["ci95"]) for k, v in a.items() if k in ("accuracy", "kappa")})


if __name__ == "__main__":
    main()
