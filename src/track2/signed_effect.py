"""Line D: direction of effect on checkpoint strength from a signed, directed network (rule: `protocol/addendum_four_lines.md`, notes for line D).

    python -m track2.signed_effect --out out

Reads the OmniPath interactions (`data/sources/omnipath/omnipath_interactions.tsv`), checks the six textbook perturbations, and, for every approved drug with a curated
mechanism of known action type, reports whether its action may strengthen, weaken or leave unclear the checkpoint. Writes `out/signed_effect.json`.
A signed-path heuristic on a literature-curated network; not a pharmacodynamic model and never a recommendation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse

from track2.direction import ACTIVATING, INHIBITORY

MCC = ("BUB1B", "BUB3", "MAD2L1")  # +1
APCC = ("CDC20", "ANAPC1", "ANAPC2", "ANAPC4", "ANAPC5", "ANAPC10", "ANAPC11", "CDC16", "CDC23", "CDC26", "CDC27")  # -1
MAX_STEPS = 3
STRENGTHEN, WEAKEN, UNCLEAR = "strengthens_checkpoint", "weakens_checkpoint", "unclear"
MIN_RECOVERED = 5
# (gene, perturbation sign, expected call): the six textbook perturbations fixed in the addendum
CONTROLS = (("TTK", -1, WEAKEN), ("BUB1", -1, WEAKEN), ("AURKB", -1, WEAKEN), ("CDC20", -1, STRENGTHEN), ("MAD2L1", -1, WEAKEN), ("TRIP13", -1, STRENGTHEN))


class SignedNetwork:
    def __init__(self, path: Path) -> None:
        edges: dict[tuple[str, str], int] = {}
        self.symbol: dict[str, str] = {}
        lines = path.read_text(encoding="utf-8").splitlines()
        col = {name: i for i, name in enumerate(lines[0].split("\t"))}
        for line in lines[1:]:
            c = line.split("\t")
            if c[col["is_directed"]] != "True" or c[col["source"]].startswith("COMPLEX") or c[col["target"]].startswith("COMPLEX"):
                continue
            stim, inh = c[col["is_stimulation"]] == "True", c[col["is_inhibition"]] == "True"
            if stim == inh:  # both or neither: the sign is not defined
                continue
            s, t = c[col["source"]], c[col["target"]]
            self.symbol.setdefault(s, c[col["source_genesymbol"]])
            self.symbol.setdefault(t, c[col["target_genesymbol"]])
            edges[(s, t)] = edges.get((s, t), 0) + (1 if stim else -1)
        self.nodes = sorted(self.symbol)
        self.index = {n: i for i, n in enumerate(self.nodes)}
        self.by_symbol = {sym: acc for acc, sym in self.symbol.items()}
        rows = [self.index[s] for (s, _), v in edges.items() if v != 0]
        cols = [self.index[t] for (_, t), v in edges.items() if v != 0]
        vals = [np.sign(v) for v in edges.values() if v != 0]
        n = len(self.nodes)
        self.adj = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=np.float64)
        self.readout = {self.by_symbol[g]: w for g, w in [(g, 1.0) for g in MCC] + [(g, -1.0) for g in APCC] if g in self.by_symbol}

    def effect(self, accession: str, sign: int) -> float | None:
        """Net effect of changing the node by `sign` (-1 inhibition/loss, +1 activation) on checkpoint strength; None if the node is not in the network."""
        if accession not in self.index:
            return None
        walks = []
        vec = np.zeros(len(self.nodes))
        vec[self.index[accession]] = 1.0
        for _ in range(MAX_STEPS):
            vec = self.adj.T.dot(vec)
            walks.append(vec)
        total = 0.0
        for acc, weight in self.readout.items():
            if acc == accession:
                total += weight * sign
                continue
            for layer in walks:  # shortest length at which a walk reaches the readout node
                if layer[self.index[acc]] != 0:
                    total += weight * sign * float(np.sign(layer[self.index[acc]]))
                    break
        return total


def call_of(net_effect: float | None) -> str:
    return UNCLEAR if not net_effect else STRENGTHEN if net_effect > 0 else WEAKEN


def validate(net: SignedNetwork) -> dict[str, Any]:
    rows = []
    for gene, sign, expected in CONTROLS:
        acc = net.by_symbol.get(gene)
        got = call_of(net.effect(acc, sign)) if acc else "not_in_network"
        rows.append({"gene": gene, "perturbation": "inhibition or loss", "expected": expected, "got": got, "recovered": got == expected,
                     "direct_readout_node": gene in MCC or gene in APCC})
    n = sum(r["recovered"] for r in rows)
    return {"controls": rows, "n_recovered": n, "used": n >= MIN_RECOVERED}


def drug_calls(net: SignedNetwork, universe: dict[str, Any]) -> dict[str, dict[str, Any]]:
    uniprot = {tid: t.get("uniprot", []) for tid, t in universe["targets"].items()}
    out: dict[str, dict[str, Any]] = {}
    for drug in universe["drugs"]:
        parts = []
        for m in drug["mechanisms"]:
            action = (m.get("action") or "").upper()
            sign = -1 if action in INHIBITORY else 1 if action in ACTIVATING else 0
            if not sign:
                continue
            for acc in uniprot.get(m["target_id"], []):
                eff = net.effect(acc, sign)
                if eff:
                    parts.append({"gene": net.symbol.get(acc, acc), "action": action, "net": eff})
        total = sum(p["net"] for p in parts)
        out[drug["drug_id"]] = {"name": drug["name"], "call": call_of(total), "net": total, "targets": parts}
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("out"))
    parser.add_argument("--omnipath", type=Path, default=Path("data/sources/omnipath/omnipath_interactions.tsv"))
    parser.add_argument("--universe", type=Path, default=Path("data/universe.json"))
    args = parser.parse_args()
    net = SignedNetwork(args.omnipath)
    validation = validate(net)
    universe = json.loads(args.universe.read_text(encoding="utf-8"))
    calls = drug_calls(net, universe)
    summary = {c: sum(1 for v in calls.values() if v["call"] == c) for c in (STRENGTHEN, WEAKEN, UNCLEAR)}
    result = {"network": {"nodes": len(net.nodes), "signed_edges": int(net.adj.nnz), "readout_nodes_found": sorted(net.symbol[a] for a in net.readout)},
              "validation": validation, "summary": summary, "drugs": calls}
    (args.out / "signed_effect.json").write_text(json.dumps(result, indent=1, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    for r in validation["controls"]:
        print(r["gene"], "expected", r["expected"], "got", r["got"], "OK" if r["recovered"] else "NOT recovered")
    print("recovered", validation["n_recovered"], "of", len(CONTROLS), "| line D used:", validation["used"], "| drugs:", summary)


if __name__ == "__main__":
    main()
