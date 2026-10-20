"""Label-support level of each tabulated drug, from the recorded reading of its US label (`config/pediatric_label_review.json`).

1 explicit: paediatric use stated in the paediatric section, with or without a minimum age.
2 by indication: established by the indication, not by the paediatric section.
3 implicit: paediatric dosing or indications imply children without stating established use or an age.
4 data or safety only: paediatric data or risks, or a paediatric section that only refers to another section (verdict `verify`).
5 not established: no relevant paediatric data. 6 not indicated or contraindicated in children.
`ALLOWED` lists the levels compatible with each verdict of the review; a test checks both against the config.
"""

ALLOWED = {"established": {1, 2, 3}, "verify": {4}, "not_established": {4, 5, 6}}

LEVEL = {**{n: 1 for n in ("BOSUTINIB", "CABOZANTINIB", "ENTRECTINIB", "EVEROLIMUS", "LAROTRECTINIB", "LONAPEGSOMATROPIN", "METFORMIN", "SELPERCATINIB", "SELUMETINIB", "SIROLIMUS", "SOMAPACITAN", "SOMATROGON", "SOMATROPIN", "THIOTEPA", "TRAMETINIB")}, **{n: 2 for n in ("TREOSULFAN", "COLCHICINE", "DAUNORUBICIN")},
         **{n: 2 for n in ("PHENYLBUTANOIC ACID",)}, **{n: 3 for n in ("DOXORUBICIN",)},
         **{n: 4 for n in ("CHLOROQUINE", "VINBLASTINE", "VINCRISTINE", "BUSULFAN", "CABAZITAXEL", "CARBOPLATIN", "CYCLOSPORINE", "DOCETAXEL ANHYDROUS", "PACLITAXEL", "BENDAMUSTINE", "BORTEZOMIB", "CARFILZOMIB", "IXABEPILONE", "LENVATINIB", "OXALIPLATIN", "POMALIDOMIDE", "SUNITINIB", "TEMOZOLOMIDE",
                                    "NINTEDANIB", "VINORELBINE", "MIDOSTAURIN", "EPIRUBICIN", "IFOSFAMIDE", "CISPLATIN", "ETOPOSIDE")},
         **{n: 6 for n in ("PAZOPANIB", "NATALIZUMAB")},
         **{n: 5 for n in ("AVAPRITINIB", "CAPMATINIB", "CARMUSTINE", "CERITINIB", "CHLORAMBUCIL", "DACARBAZINE", "DISOPYRAMIDE", "ERIBULIN", "ETOPOSIDE PHOSPHATE", "GILTERITINIB", "IDARUBICIN",
                                    "IXAZOMIB CITRATE", "LOMUSTINE", "MECHLORETHAMINE", "MELPHALAN", "MITOXANTRONE", "PEXIDARTINIB", "QUIZARTINIB", "REGORAFENIB", "RIPRETINIB", "SELINEXOR", "SORAFENIB", "VALRUBICIN",
                                    "VANDETANIB")}}
