"""Match LINCS perturbagen names to approved-drug names (ChEMBL parent names) after removing salts."""

from __future__ import annotations

import re
from typing import Any

SALTS = {
    "HYDROCHLORIDE", "HYDROBROMIDE", "SODIUM", "POTASSIUM", "CALCIUM", "MAGNESIUM", "MESYLATE", "MESILATE", "CITRATE", "TARTRATE",
    "MALEATE", "FUMARATE", "SUCCINATE", "SULFATE", "SULPHATE", "PHOSPHATE", "ACETATE", "BESYLATE", "BESILATE", "TOSYLATE",
    "BROMIDE", "CHLORIDE", "NITRATE", "OXALATE", "LACTATE", "GLUCONATE", "HCL", "HBR", "DIHYDROCHLORIDE", "MONOHYDRATE",
    "DIHYDRATE", "HYDRATE", "ANHYDROUS", "HEMIHYDRATE", "TRIHYDRATE", "MEGLUMINE", "PAMOATE", "DIBASIC", "MONOBASIC", "TRIBASIC",
    "EDISYLATE", "ESYLATE", "CAMSYLATE", "NAPSYLATE", 
    "HEMISULFATE", "DISODIUM", "TROMETHAMINE", "TRIHYDRATE", "SESQUIHYDRATE", "ETHANOLATE", "FORMATE", "XINAFOATE",
}


def clean(name: str) -> str:
    """Upper-case, drop parenthesised text and salt/hydrate words, keep letters and digits only."""
    text = re.sub(r"\([^)]*\)", " ", str(name).upper())
    tokens = [t for t in re.split(r"[\s,/]+", text) if t and t not in SALTS]
    return re.sub(r"[^A-Z0-9]", "", "".join(tokens))


def drug_index(universe: list[dict[str, Any]]) -> dict[str, list[str]]:
    """Cleaned name -> drug_ids (more than one id means the match is ambiguous)."""
    index: dict[str, list[str]] = {}
    for drug in universe:
        key = clean(drug["name"])
        if key and drug["drug_id"] not in index.setdefault(key, []):
            index[key].append(drug["drug_id"])
    return index


def match(pert_name: str, index: dict[str, list[str]]) -> list[str]:
    return index.get(clean(pert_name), [])
