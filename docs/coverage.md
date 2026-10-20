# Test coverage (offline suite)

Measured on 2026-10-06, Windows 11, Python 3.10.11, pytest-cov (not in the lock file; `pip install pytest-cov`).

```bash
cd HSelith_MVA-Track2
PYTHONPATH=src python -m pytest --cov=track2 --cov-report=term-missing -q
```

Result: 25 tests passed; **overall statement coverage 32 % (1,163 of 1,717 statements not executed)**. The tests are the offline unit tests only; no network, no patient data. Branch coverage was not measured.

| File | Statements | Missed | Cover |
|---|--:|--:|--:|
| `build_candidates.py` | 99 | 99 | 0 % |
| `build_report.py` | 76 | 76 | 0 % |
| `chembl_universe.py` | 57 | 57 | 0 % |
| `direction.py` | 102 | 19 | 81 % |
| `evidence.py` | 85 | 7 | 92 % |
| `gates.py` | 36 | 2 | 94 % |
| `genome_covariates.py` | 85 | 56 | 34 % |
| `http_cache.py` | 49 | 29 | 41 % |
| `labels.py` | 103 | 49 | 52 % |
| `modules.py` | 28 | 28 | 0 % |
| `proximity.py` | 54 | 0 | 100 % |
| `review_page.py` | 92 | 11 | 88 % |
| `run_pipeline.py` | 48 | 48 | 0 % |
| `scoring.py` | 38 | 38 | 0 % |
| `snapshot.py` | 63 | 10 | 84 % |
| `string_net.py` | 69 | 69 | 0 % |
| `toxicity.py` | 29 | 10 | 66 % |
| `validation.py` | 113 | 113 | 0 % |
| `vendor/consequence_calculator.py` | 185 | 167 | 10 % |
| `vendor/genome_wide_io.py` | 122 | 122 | 0 % |
| `vendor/genome_wide_search.py` | 184 | 153 | 17 % |
| **Total** | **1,717** | **1,163** | **32 %** |

`direction.py` and `evidence.py` were added by other work in the same tree while this report was produced; the numbers above are those of the run at the time and will change as files are added. Re-run the command for current figures.

## What is not tested, and why

- **Need network (external services):** `chembl_universe.py`, `string_net.py`, the download paths of `http_cache.py`, the DailyMed/PubMed/openFDA fetches in `labels.py`, `toxicity.py` and `build_candidates.py`. Testing them offline would require recorded fixtures, which do not exist yet.
- **Need the full cached data or hours of compute:** `run_pipeline.py`, `modules.py`, `scoring.py`, `validation.py` (about 30-40 minutes), `build_report.py`. They are exercised only by running the real pipeline; the shipped outputs in `out/` are their evidence, not a unit test.
- **Need patient data (VCF), not distributed:** `genome_covariates.py` and `vendor/*` (the Track 1 modules). Only small pure helpers are touched.
- The tested core is the numerical method (`proximity.py`), the requirement rules and label parsing (partly), the schema, the data manifest and the review page.

The 32 % figure therefore says that most of the pipeline logic is verified only end to end, not by unit tests. Nothing was added to raise the number artificially.
