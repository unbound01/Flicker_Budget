#!/usr/bin/env python3
"""Print or verify SHA-256 checksums of the regenerated data files (data/*.json).
The digitized data are not distributed (they derive from copyrighted figures); a reader who regenerates them
from the same PDFs with the digitizers can check they match data/CHECKSUMS.sha256.
Usage: python scripts/data_checksums.py [--write]"""
import glob, hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, "..", "data"); OUT = os.path.join(D, "CHECKSUMS.sha256")
def digest(p):   # hash the numbers, not the formatting, so a different json writer gives the same result
    return hashlib.sha256(json.dumps(json.load(open(p)), sort_keys=True, separators=(",", ":")).encode()).hexdigest()
files = sorted(glob.glob(os.path.join(D, "*.json")))
lines = [f"{digest(p)}  {os.path.basename(p)}" for p in files]
if "--write" in sys.argv:
    open(OUT, "w").write("\n".join(lines) + "\n"); print("wrote", OUT)
else:
    ref = dict(l.split()[::-1] for l in open(OUT).read().split("\n") if l.strip()) if os.path.exists(OUT) else {}
    for l in lines:
        h, n = l.split(); print(("OK      " if ref.get(n) == h else "DIFFERS ") + n)
