r"""General (disease-agnostic) mechanism classes of a gene from UniProt keywords and Gene Ontology term names.

Pre-registration: protocol/addendum_engine_v2.md (section 4). Complements, and does not replace, the narrow rules of
`disease_landscape.classify_terms` (which this module never imports or edits).

Input: the list of annotation terms of the reviewed human UniProt entry (UniProt keywords + GO term names), as returned by
`disease_landscape.uniprot_annotation(symbol)["terms"]`.

Rule table (fixed, deterministic, no gene or disease names)
  * RULES maps a class to a list of (regular expression, weight). Weight 2 = strong term, weight 1 = weak term.
  * Every distinct term is tested against every pattern of a class (case-insensitive search); a term contributes the
    HIGHEST weight among the patterns of that class that match it (a term is counted once per class).
  * The class score is the sum of its term contributions, capped at MAX_CLASS_SCORE (6). The class with the highest score wins; ties are broken by the
    order of ORDER (specific classes before generic ones). A best score below MIN_SCORE (2) gives 'other'.
"""
from __future__ import annotations

import re
from typing import Dict, Iterable, List, Sequence, Tuple

MIN_SCORE = 2
MAX_CLASS_SCORE = 6  # cap: pleiotropic GO annotation must not let a generic class outweigh a specific one by sheer term count
OTHER = "other"

# order = tie-break order (specific first)
ORDER: Tuple[str, ...] = (
    "mitochondrial",
    "dna_repair_replication",
    "cell_cycle_mitotic_checkpoint",
    "rna_processing_translation",
    "lysosomal_vesicular_trafficking",
    "ion_channel_transporter",
    "transcription_chromatin_regulation",
    "structural_cytoskeletal",
    "extracellular_secreted",
    "metabolic_enzyme",
    "receptor_signalling",
)
CLASSES: Tuple[str, ...] = ORDER + (OTHER,)

_ENZYME_ACTIVITY = (r"^(?!.*(protein|kinase|ubiquitin|DNA|RNA|peptid|nucle|GTP|ATP)).*"
                    r"(transferase|hydrolase|oxidoreductase|dehydrogenase|synthase|synthetase|lyase|isomerase|ligase|"
                    r"carboxylase|reductase|oxidase|phosphatase|mutase|epimerase|racemase|aminotransferase|esterase|"
                    r"lipase|glycosidase|sulfatase|catalytic) activity")

RULES: Dict[str, List[Tuple[str, int]]] = {
    "mitochondrial": [
        (r"respiratory chain|oxidative phosphorylation|electron transport|respirasome|ATP synthesis coupled", 2),
        (r"mitochondri\w* (translation|respiratory|membrane organization|DNA|fission|fusion|ribosome|"
         r"protein import|gene expression|tRNA|electron|complex|calcium ion|transport|organization|"
         r"localization|ATP)", 2),
        (r"protein (import into|localization to) mitochondri", 2),
        (r"^NADH dehydrogenase|cytochrome-c oxidase|cytochrome c oxidase|respiratory chain complex", 2),
        (r"^mitochondri(on|al)\b", 1),
        (r"^Mitochondrion (inner|outer) membrane$|^Transit peptide$", 1),
    ],
    "dna_repair_replication": [
        (r"DNA repair|DNA replication|double-strand break|homologous recombination|interstrand cross-link|"
         r"nucleotide-excision|base-excision|mismatch repair|replication fork|translesion|DNA synthesis|"
         r"DNA recombination|DNA polymerase|DNA ligase|DNA topological|DNA metabolic|DNA integrity", 2),
        (r"DNA damage|telomere|helicase|nuclease|topoisomerase|site of double-strand", 1),
    ],
    "cell_cycle_mitotic_checkpoint": [
        (r"mitosis|mitotic|^Cell division$|cell division|spindle|kinetochore|centromere|centrosome|centriole|"
         r"chromosome segregation|sister chromatid|cohesin|checkpoint|cyclin|G1/S|G2/M|cytokinesis|meiosis|"
         r"meiotic|chromosome condensation", 2),
        (r"^Cell cycle$|cell cycle|cell cycle arrest", 1),
    ],
    "rna_processing_translation": [
        (r"ribosom|(?<!post-)translation(?!al)|rRNA|tRNA|mRNA (processing|splicing|export|decay|3'-end|stability)|"
         r"spliceosom|RNA splicing|RNA processing|nonsense-mediated|snRNA|snoRNA|RNA modification|"
         r"translation initiation|translation elongation|RNA helicase|exosome \(RNase", 2),
        (r"RNA binding|RNA-binding|^Ribonucleoprotein$|ribonucleoprotein complex|nucleolus|mRNA binding", 1),
    ],
    "lysosomal_vesicular_trafficking": [
        (r"lysosom|vesicle-mediated transport|endocytosis|exocytosis|vesicle (fusion|docking|coating|budding|"
         r"organization|transport)|autophag|Golgi vesicle transport|ER to Golgi|membrane trafficking|SNARE|clathrin|"
         r"COPI|COPII|retromer|intracellular protein transport|protein transport|endosom\w* (transport|organization|"
         r"to)|phagosome|synaptic vesicle (cycle|endocytosis|exocytosis|priming)|vacuol", 2),
        (r"^Endosome$|endosome|^Lysosome$|cytoplasmic vesicle|synaptic vesicle|^Golgi apparatus$|^Protein transport$|"
         r"Golgi membrane|^Transport$", 1),
    ],
    "ion_channel_transporter": [
        (r"channel|ion transport|transmembrane transporter activity|transmembrane transport\b|symporter|antiporter|"
         r"solute carrier|ABC transporter|ATPase-coupled|P-type|membrane transport|voltage-gated|ligand-gated|"
         r"ion homeostasis|^Symport$|^Antiport$|uniporter|exchanger|^Sodium transport$|^Potassium transport$|"
         r"^Calcium transport$", 2),
        (r"^Transport$|transporter|^Chloride$|^Sodium$|^Potassium$", 1),
    ],
    "transcription_chromatin_regulation": [
        (r"^Transcription( regulation)?$|transcription (factor|coactivator|corepressor|cis-regulatory|regulator)|"
         r"DNA-binding transcription|regulation of (DNA-templated )?transcription|RNA polymerase II.*(binding|"
         r"transcription|regulatory)|^Chromatin regulator$|chromatin (organization|remodel|binding)|histone|"
         r"nucleosome|epigenetic|transcription repressor|transcription activator", 2),
        (r"^Activator$|^Repressor$|^DNA-binding$|^Zinc-finger$|^chromatin$|^Developmental protein$|DNA binding|"
         r"sequence-specific.*DNA binding", 1),
    ],
    "structural_cytoskeletal": [
        (r"cytoskelet|actin|microtubule|tubulin|intermediate filament|keratin|myosin|sarcomere|myofibril|"
         r"striated muscle|muscle (contraction|filament)|structural constituent of (?!ribosome)|cilium|ciliary|"
         r"axoneme|dynein|kinesin|Z disc|costamere|desmosome|dystrophin|sarcolemma|lamin|nuclear envelope organization", 2),
        (r"cell adhesion|focal adhesion|cell junction|adherens|tight junction|cell projection|^Cell junction$|"
         r"^Cell adhesion$|^Cytoskeleton$|cytoskeleton|cell cortex|lamellipodium", 1),
    ],
    "extracellular_secreted": [
        (r"^Secreted$|extracellular matrix|extracellular space|blood coagulation|^Hemostasis$|hemostasis|fibrinolysis|"
         r"complement activation|collagen|basement membrane|blood microparticle|proteinaceous extracellular", 2),
        (r"extracellular region|^Signal$|^Glycoprotein$|^Disulfide bond$|^Zymogen$|^Protease$|serine protease", 1),
    ],
    "metabolic_enzyme": [
        (_ENZYME_ACTIVITY, 2),
        (r"metabolic process|biosynthetic process|catabolic process(?<!protein catabolic process)|"
         r"(?<!protein )(?<!lipid )metabolism$|^Lipid metabolism$|biosynthesis$|^Amino-acid|glycoly|gluconeo|glycogen|"
         r"fatty acid|cholesterol (metabolic|biosynthetic|homeostasis)|sphingolipid|glycosaminoglycan|"
         r"carbohydrate|vitamin|cofactor|one-carbon|urea cycle|purine|pyrimidine|nucleotide (metabolic|biosynthetic)|"
         r"heme (biosynthetic|metabolic)|oxidation-reduction process|catalytic activity|detoxification|"
         r"xenobiotic metabolic", 2),
        (r"^Hydrolase$|^Transferase$|^Oxidoreductase$|^Lyase$|^Isomerase$|^Ligase$|^Heme$|^Iron$|"
         r"^Pyridoxal phosphate$|^NAD$|^FAD$|^Flavoprotein$|heme binding", 1),
    ],
    "receptor_signalling": [
        (r"^Receptor$|receptor activity|signaling receptor|receptor (binding|complex|signaling pathway)|"
         r"signal transduction|G protein|GTPase|^Kinase$|kinase activity|-protein kinase$|MAPK|cascade|"
         r"JAK-STAT|NF-kappaB|Wnt|Notch|Ras protein|phosphatidylinositol 3-kinase|^Protein phosphatase$|"
         r"protein phosphorylation|neurotransmitter|hormone activity|growth factor|"
         r"cytokine activity|cytokine-mediated|interleukin|interferon", 2),
        (r"^Cytokine$|^Hormone$|^Growth factor$|^Immunity$|^Innate immunity$|^Adaptive immunity$|immune response|"
         r"inflammatory response|^SH2 domain|SH3 domain|^Proto-oncogene$", 1),
    ],
}

_COMPILED: Dict[str, List[Tuple["re.Pattern[str]", int]]] = {
    cls: [(re.compile(p, re.I), w) for p, w in rules] for cls, rules in RULES.items()
}


def score_terms(terms: Iterable[str]) -> Dict[str, int]:
    """Class -> score (sum over distinct terms of the highest matching weight in that class)."""
    distinct = sorted({t.strip() for t in terms if t and t.strip()})
    scores: Dict[str, int] = {c: 0 for c in ORDER}
    for t in distinct:
        for cls in ORDER:
            best = 0
            for rx, w in _COMPILED[cls]:
                if w > best and rx.search(t):
                    best = w
            scores[cls] += best
    return {c: min(v, MAX_CLASS_SCORE) for c, v in scores.items()}


def matched_terms(terms: Iterable[str], cls: str) -> List[str]:
    """Distinct terms (sorted) that contributed to `cls`."""
    if cls not in _COMPILED:
        return []
    return sorted({t.strip() for t in terms if t and t.strip() and any(rx.search(t) for rx, _ in _COMPILED[cls])})


def classify_general(terms: Sequence[str]) -> Tuple[str, List[str]]:
    """(general mechanism class, matched terms of that class). 'other' when the best score is below MIN_SCORE."""
    scores = score_terms(terms)
    best = max(scores.values()) if scores else 0
    if best < MIN_SCORE:
        return OTHER, []
    for cls in ORDER:  # fixed order = tie-break
        if scores[cls] == best:
            return cls, matched_terms(terms, cls)
    return OTHER, []  # pragma: no cover
