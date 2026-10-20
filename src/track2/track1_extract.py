"""Gene-level extract of the Track 1 output, so that the modules that read it can be rerun without the Track 1 repository.

    python -m track2.track1_extract          # needs the Track 1 repository next to this one; writes config/track1_gene_level_extract.json

The extract holds only what `track2.differential` and `track2.variant_level_consequence` already use, at gene level and in categorical or population terms that the Track 1 report publishes
(candidate type, number of high-impact variants, log10 joint rarity, phenotype score and rank, finding types, consequence classes). It holds no position, allele, genotype call or read metric:
the loaders of those two modules read only the whitelisted columns, and `tests/test_track1_extract.py` asserts that the extract carries no such string and that rerunning both modules from
the extract reproduces their outputs exactly.
"""

from __future__ import annotations

import json
from pathlib import Path

from track2 import differential as diff
from track2 import variant_level_consequence as vlc

ROOT = Path(__file__).resolve().parents[2]
EXTRACT = ROOT / "config" / "track1_gene_level_extract.json"


def build() -> dict:
    t1 = diff.TRACK1
    ranking = diff.load_ranking(t1 / "candeleroj98_genomewide_ranking.csv")
    structural = {g: sorted(v) for g, v in diff.load_structural(t1 / "candeleroj98_gatk-manta-fusion.csv").items()}
    pheno = diff.load_phenotype(t1 / "candeleroj98_phenotype_summary.json")
    snv = vlc.read_ranking_classes(t1 / "candeleroj98_genomewide_ranking.csv")
    sv, decay = vlc.read_structural(t1 / "candeleroj98_gatk-manta-fusion.csv")
    pheno_decay = vlc.read_phenotype_decay(t1 / "candeleroj98_phenotype_summary.json")
    return {
        "description": "Gene-level extract of the Track 1 output used by track2.differential and track2.variant_level_consequence; no positions, alleles, genotype calls or read metrics.",
        "differential": {"ranking": ranking, "structural": structural, "phenotype": pheno},
        "variant_level": {"snv_classes": snv, "sv_classes": sv, "decay_states": decay, "phenotype_has_decay": pheno_decay},
    }


def load() -> dict:
    return json.loads(EXTRACT.read_text(encoding="utf8"))


def differential_inputs() -> tuple[dict, dict, dict]:
    d = load()["differential"]
    return d["ranking"], {g: set(v) for g, v in d["structural"].items()}, d["phenotype"]


def variant_level_inputs() -> tuple[dict, dict, dict, bool]:
    v = load()["variant_level"]
    return v["snv_classes"], v["sv_classes"], v["decay_states"], v["phenotype_has_decay"]


def main() -> None:
    EXTRACT.write_text(json.dumps(build(), indent=2, sort_keys=True) + "\n", encoding="utf8")
    print("wrote", EXTRACT)


if __name__ == "__main__":
    main()
