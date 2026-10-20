"""State, for each candidate, what it is proposed for (rules in `config/repositioning_rationale.json`).

The text comes from fixed rules over the drug's own ChEMBL mechanism text and its measured links, tied to the five consequences of the variant. Nothing is written drug by
drug. A hypothesis for follow-up, not medical advice.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

CONFIG = Path("config/repositioning_rationale.json")


def load(path: Path = CONFIG) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _antigen(moa: list[str], aliases: dict[str, str] | None = None) -> str:
    """The non-tubulin binding-agent target of an antibody-drug conjugate, taken from its mechanism text (readable alias when one is listed in the config)."""
    names = [m.replace(" binding agent", "") for m in moa if "binding agent" in m.lower()]
    return ", ".join((aliases or {}).get(n, n) for n in names) or "its target antigen"


def _matches(rule: dict[str, Any], text: str, drug_list: str | None = None) -> bool:
    if rule.get("only_for_list") and rule["only_for_list"] != drug_list:
        return False
    if "match_all_groups" in rule:
        return all(any(k in text for k in group) for group in rule["match_all_groups"])
    return any(k in text for k in rule["match_any"])


def proposed_use(moa: list[str], proximal_genes: list[str], novelty: str, config: dict[str, Any], drug_list: str | None = None) -> dict[str, Any]:
    text = " ".join(moa).lower()
    for rule in config["rules"]:
        if _matches(rule, text, drug_list):
            return {"rule": rule["id"], "consequences": rule["consequences"], "context": rule["context"],
                    "new_use": rule["new_use"].format(antigen=_antigen(moa, config.get("antigen_aliases"))), "caveat": rule["caveat"], "validation": rule.get("validation", "")}
    if novelty == "new_via_expanded_links" and proximal_genes:
        d = config["default_measured_links"]
        return {"rule": "measured_links", "consequences": d["consequences"], "context": d["context"],
                "new_use": d["new_use"].format(genes=", ".join(proximal_genes[:3])), "caveat": d["caveat"], "validation": d.get("validation", "")}
    d = config["default_none"]
    return {"rule": "none", "consequences": [], "context": d["context"], "new_use": d["new_use"], "caveat": d["caveat"], "validation": d.get("validation", "")}
