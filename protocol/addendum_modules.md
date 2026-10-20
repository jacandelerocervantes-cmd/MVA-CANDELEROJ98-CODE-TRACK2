# Addendum to the frozen protocol: additional modules

Written on 2026-10-06 **before** running the additional modules on real data (only the main spindle-checkpoint module had been run). It does not modify `frozen_protocol.md`.

- Three modules for **downstream consequences** of BUB1B loss are added, each defined only by verified Reactome pathways and a few core genes: proteotoxic stress (R-HSA-381119, R-HSA-3371453, R-HSA-983168), p53 response and senescence (R-HSA-2559583, R-HSA-69541, R-HSA-6804116; TP53, MDM2, CDKN1A, CDKN2A) and mTORC1/AMPK growth signalling (R-HSA-166208, R-HSA-380972; MTOR, RPTOR, TSC1, TSC2, PRKAA1, PRKAA2).
- The replication demo uses the Fanconi pathway (R-HSA-6783310).
- The same parameters, the same P1 test and the same decision rules of the frozen protocol are applied to each module. No parameter is adjusted after seeing results.
- Observation prior to this addendum, made with the main module: ChEMBL has no curated mechanism for some very well-known drugs (for example metformin and vincristine), so they cannot be scored with this method. This is declared as a coverage limitation and is not compensated for with hand-entered values.
