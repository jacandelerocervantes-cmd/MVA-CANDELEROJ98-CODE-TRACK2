"""Scripted PubMed search for the literature rationale of drugs supported by one usable line (rule: `protocol/addendum_four_lines.md`, amendment 3).

    python -m track2.literature_search --out out

For every drug of `out/four_lines.json` that is supported by one usable line and passes the gates, and for each curated human target gene of the drug, runs

    GENE[Title/Abstract] AND ("spindle assembly checkpoint" OR BUBR1 OR BUB1B OR "mosaic variegated aneuploidy" OR "chromosome missegregation" OR aneuploidy)

with NCBI E-utilities and stores the count, the first PMIDs and their titles. The count is not the rationale: a drug has one only if the author reads the abstracts and records the
supported claim in `config/literature_rationale.json`. Writes `out/literature_search.json`.
"""

from __future__ import annotations

import argparse
import json
import time
import urllib.parse
from pathlib import Path
from typing import Any

from track2.http_cache import get_json

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
TOPIC = '("spindle assembly checkpoint"[Title/Abstract] OR BUBR1[Title/Abstract] OR BUB1B[Title/Abstract] OR "mosaic variegated aneuploidy"[Title/Abstract] OR "chromosome missegregation"[Title/Abstract] OR aneuploidy[Title/Abstract])'
RETMAX = 5


def gene_symbols(accessions: list[str]) -> dict[str, str]:
    """UniProt accession -> primary gene symbol (reviewed human entries), in batches."""
    out: dict[str, str] = {}
    for i in range(0, len(accessions), 40):
        batch = accessions[i:i + 40]
        query = urllib.parse.quote("(" + " OR ".join(f"accession:{a}" for a in batch) + ")")
        for r in get_json(f"{UNIPROT}?query={query}&fields=accession,gene_primary&format=json&size=100")["results"]:
            genes = r.get("genes") or []
            if genes and genes[0].get("geneName"):
                out[r["primaryAccession"]] = genes[0]["geneName"]["value"]
    return out


def search_gene(gene: str) -> dict[str, Any]:
    term = urllib.parse.quote(f"{gene}[Title/Abstract] AND {TOPIC}")
    res = get_json(f"{ESEARCH}?db=pubmed&term={term}&retmax={RETMAX}&retmode=json&sort=relevance")["esearchresult"]
    ids = res.get("idlist", [])
    titles = {}
    if ids:
        summary = get_json(f"{ESUMMARY}?db=pubmed&id={','.join(ids)}&retmode=json")["result"]
        titles = {i: summary[i]["title"] for i in ids if i in summary}
    time.sleep(0.4)  # NCBI allows three requests per second without a key
    return {"count": int(res.get("count", 0)), "pmids": ids, "titles": titles}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    args = parser.parse_args()
    four = json.loads((args.out / "four_lines.json").read_text(encoding="utf-8"))
    universe = {d["drug_id"]: d for d in json.loads(args.universe.read_text(encoding="utf-8"))["drugs"]}
    drugs = [r for r in four["tier_c_one_line"]]
    accessions = sorted({a for r in drugs for a in universe[r["drug_id"]]["human_targets"]})
    symbols = gene_symbols(accessions)
    genes = sorted(set(symbols.values()))
    print(len(drugs), "drugs with one usable line;", len(accessions), "targets;", len(genes), "genes", flush=True)
    by_gene = {g: search_gene(g) for g in genes}
    result = {"query_topic": TOPIC, "genes": by_gene, "drugs": {}}
    for r in drugs:
        targets = sorted({symbols[a] for a in universe[r["drug_id"]]["human_targets"] if a in symbols})
        result["drugs"][r["drug_id"]] = {"name": r["name"], "supported_by": r["supported_by"], "target_genes": targets,
                                         "counts": {g: by_gene[g]["count"] for g in targets}, "n_papers_any_target": sum(by_gene[g]["count"] for g in targets)}
    (args.out / "literature_search.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    ranked = sorted(result["drugs"].values(), key=lambda d: -d["n_papers_any_target"])
    print([(d["name"], d["n_papers_any_target"]) for d in ranked[:25]])


if __name__ == "__main__":
    main()
