"""Standalone embedding worker for engine v4 (phase 1). Run ONLY in isolated mode:

    python -I src/track2/embed_worker.py <model_dir> <mean|cls> <texts.json> <out.npy> [max_len] [batch]

-I ignores PYTHON* variables, the user site and the script directory, so nothing next to the downloaded model files can
be imported. This file imports no project module (stdlib, numpy, torch and transformers only). It reads a JSON list of
strings, writes an (n, hidden) float32 array of L2-normalised vectors (cosine = dot product), and prints one JSON line with
the number of texts and the seconds. Weights are loaded with local_files_only=True; the pickle-based .bin checkpoint of
BiomedBERT is read by torch with weights_only=True (the default of the installed torch).
"""
import json
import sys
import time


def main() -> None:
    if len(sys.argv) < 5:
        raise SystemExit("usage: python -I embed_worker.py model_dir mean|cls texts.json out.npy [max_len] [batch]")
    model_dir, pool, src, dst = sys.argv[1:5]
    max_len = int(sys.argv[5]) if len(sys.argv) > 5 else 128
    batch = int(sys.argv[6]) if len(sys.argv) > 6 else 32
    if pool not in ("mean", "cls"):
        raise SystemExit("pool must be mean or cls")
    import numpy as np
    import torch
    from transformers import AutoModel, AutoTokenizer

    with open(src, encoding="utf-8") as fh:
        texts = json.load(fh)
    t0 = time.time()
    tok = AutoTokenizer.from_pretrained(model_dir, local_files_only=True)
    model = AutoModel.from_pretrained(model_dir, local_files_only=True).eval()
    torch.set_num_threads(max(1, torch.get_num_threads()))
    order = sorted(range(len(texts)), key=lambda i: len(texts[i]))   # length-sorted batches: little padding
    out = np.zeros((len(texts), model.config.hidden_size), dtype=np.float32)
    with torch.no_grad():
        for s in range(0, len(order), batch):
            ids = order[s:s + batch]
            enc = tok([texts[i] for i in ids], padding=True, truncation=True, max_length=max_len, return_tensors="pt")
            h = model(**enc).last_hidden_state
            if pool == "cls":
                v = h[:, 0, :]
            else:
                m = enc["attention_mask"].unsqueeze(-1).to(h.dtype)
                v = (h * m).sum(1) / m.sum(1).clamp(min=1)
            v = torch.nn.functional.normalize(v, dim=-1)
            out[ids] = v.numpy()
            if (s // batch) % 50 == 0:
                print(json.dumps({"done": s, "of": len(texts)}), file=sys.stderr, flush=True)
    np.save(dst, out)
    print(json.dumps({"n_texts": len(texts), "seconds": round(time.time() - t0, 2), "model_dir": model_dir,
                      "pool": pool, "max_len": max_len}), flush=True)


if __name__ == "__main__":
    main()
