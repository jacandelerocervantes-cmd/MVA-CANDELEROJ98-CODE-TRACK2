# Note on the freeze

The protocol (`frozen_protocol.md`, SHA-256 in `frozen_protocol.sha256`) was written on 2026-10-06 after downloading the ChEMBL approved-drug universe (`data/universe.json`, 13:02) and downloading the STRING network, and before computing any proximity on real data. Up to that point only the drug count and the tests with synthetic graphs had been seen; no proximity, ranking or validation result.

## Translation note

The protocol was first written in Spanish and hashed on 2026-10-06 (SHA-256 `7396e9478c9d44840baf0ed204425104b6615cd247e483aa997c21bd624a0835`) before any proximity was computed. It was then translated sentence by sentence into English, with every parameter, rule and number unchanged, so that the whole repository is in English; the hash in `frozen_protocol.sha256` is the hash of the English text (`b56795bb30c0c2e9f6ca02706b21c941c4bfa76b3e3861314d794dfc6e1e1131`). The original Spanish hash is kept here for the record.
