#!/usr/bin/env python3 -I
"""
Extract MS2 fragment-ion XICs (and MS1 precursor XICs) for SLC25A13 and shared
SLC25A12/SLC25A13 tryptic peptides from DIA mzML files (ThermoRawFileParser output).

Usage:
  python3 -I slc25a13_fragment_xic.py --mzml SM_1.mzML --sample SM_1 --out xic_dir [--ppm 10]

For every MS2 spectrum whose isolation window contains the precursor m/z, the
summed intensity within +/- ppm of each fragment m/z is recorded. Fragments are
generated automatically: singly charged y3..y(n-1) and b3..b6, plus doubly
charged y ions for 3+ precursors, restricted to the acquired MS2 m/z range
(200-1800). From MS1 spectra the precursor monoisotopic and M+1 XICs are
extracted (fragment labels "M0", "M1"). Output: one tidy CSV per sample with
columns sample, peptide, charge, fragment, fragment_mz, rt_min, intensity.
No peak picking and no smoothing are applied here.
"""
import argparse, base64, zlib, csv, os, sys
import numpy as np
from lxml import etree

NS = "{http://psi.hupo.org/ms/mzml}"
MONO = {'G':57.02146,'A':71.03711,'S':87.03203,'P':97.05276,'V':99.06841,'T':101.04768,
        'C':103.00919,'L':113.08406,'I':113.08406,'N':114.04293,'D':115.02694,'Q':128.05858,
        'K':128.09496,'E':129.04259,'M':131.04049,'H':137.05891,'F':147.06841,'R':156.10111,
        'Y':163.06333,'W':186.07931}
H2O = 18.010565; PROTON = 1.007276; C13 = 1.003355
MS2_RANGE = (200.0, 1800.0)

# (peptide, precursor charge as quantified by DIA-NN, class)
TARGETS = [
    ("EEFVLAAQK",          2, "SLC25A13-specific"),
    ("FGLYLPLFKPSVSTSK",   3, "SLC25A13-specific"),
    ("INLPAPNPDHVGGYK",    3, "SLC25A13-specific"),
    ("LAVATFAGIENK",       2, "SLC25A13-specific"),
    ("SSPQFGVTLLTYELLQR",  3, "SLC25A13-specific"),
    ("STGSFVGELMYK",       2, "SLC25A13-specific"),
    ("TVELLSGVVDQTK",      2, "SLC25A13-specific"),
    ("GLLPQLLGVAPEK",      2, "I/L-isobaric with SLC25A12"),
    ("KIYSTLAGTR",         2, "shared SLC25A12/13"),
    ("LQVAGEITTGPR",       2, "shared SLC25A12/13"),
    ("YEGFFGLYR",          2, "shared SLC25A12/13"),
    # SLC25A12 (aralar) counterparts, as paralog-presence controls
    ("LATATFAGIENK",       2, "SLC25A12-specific"),
    ("IDLPPANPDHIGGYR",    3, "SLC25A12-specific"),
    ("EEFAQSAIR",          2, "SLC25A12-specific"),
    ("IVQLLAGVADQTK",      2, "SLC25A12-specific"),
    ("GSGSVVGELMYK",       2, "SLC25A12-specific"),
    ("DLGIFGLYK",          2, "SLC25A12-specific"),
]

def pep_mass(seq): return sum(MONO[a] for a in seq) + H2O
def prec_mz(seq, z): return (pep_mass(seq) + z*PROTON)/z
def frag_mz(seq, label):
    """label like y5, b3, y9++ (doubly charged)."""
    z = 2 if label.endswith("++") else 1
    core = label.rstrip("+")
    kind, n = core[0], int(core[1:])
    if kind == "y": m = sum(MONO[a] for a in seq[-n:]) + H2O
    elif kind == "b": m = sum(MONO[a] for a in seq[:n])
    else: raise ValueError(label)
    return (m + z*PROTON)/z

def fragment_labels(seq, z):
    n = len(seq); labs = []
    labs += [f"y{i}" for i in range(3, n)]
    labs += [f"b{i}" for i in range(3, min(6, n-1)+1)]
    if z >= 3: labs += [f"y{i}++" for i in range(6, n)]
    out = []
    for l in labs:
        mz = frag_mz(seq, l)
        if MS2_RANGE[0] <= mz <= MS2_RANGE[1]: out.append((l, mz))
    return out

def decode(binary_el):
    dtype = None; comp = False
    for cv in binary_el.iter(NS+"cvParam"):
        acc = cv.get("accession")
        if acc == "MS:1000523": dtype = np.float64
        elif acc == "MS:1000521": dtype = np.float32
        elif acc == "MS:1000574": comp = True
    raw = base64.b64decode(binary_el.find(NS+"binary").text or "")
    if comp: raw = zlib.decompress(raw)
    return np.frombuffer(raw, dtype=dtype)

def spectrum_info(spec):
    ms_level = None; rt = None; iso_lo = iso_hi = None; mz = inten = None
    for cv in spec.findall(NS+"cvParam"):
        if cv.get("accession") == "MS:1000511": ms_level = int(cv.get("value"))
    scan = spec.find(NS+"scanList/"+NS+"scan")
    if scan is not None:
        for cv in scan.findall(NS+"cvParam"):
            if cv.get("accession") == "MS:1000016":
                rt = float(cv.get("value"))
                if cv.get("unitName","").startswith("second"): rt /= 60.0
    iw = spec.find(NS+"precursorList/"+NS+"precursor/"+NS+"isolationWindow")
    if iw is not None:
        tgt = lo = hi = None
        for cv in iw.findall(NS+"cvParam"):
            a = cv.get("accession"); v = float(cv.get("value"))
            if a == "MS:1000827": tgt = v
            elif a == "MS:1000828": lo = v
            elif a == "MS:1000829": hi = v
        if tgt is not None and lo is not None and hi is not None:
            iso_lo, iso_hi = tgt - lo, tgt + hi
    for bda in spec.findall(NS+"binaryDataArrayList/"+NS+"binaryDataArray"):
        accs = {cv.get("accession") for cv in bda.findall(NS+"cvParam")}
        arr = decode(bda)
        if "MS:1000514" in accs: mz = arr
        elif "MS:1000515" in accs: inten = arr
    return ms_level, rt, iso_lo, iso_hi, mz, inten

def window_sum(mz, inten, target, ppm):
    tol = target * ppm * 1e-6
    i0 = np.searchsorted(mz, target - tol); i1 = np.searchsorted(mz, target + tol)
    return float(inten[i0:i1].sum()) if i1 > i0 else 0.0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mzml", required=True); ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--ppm", type=float, default=10.0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    targets = []
    for seq, z, cls in TARGETS:
        pmz = prec_mz(seq, z)
        targets.append((seq, z, pmz, fragment_labels(seq, z),
                        [("M0", pmz), ("M1", pmz + C13/z)]))
    rows = []; n_ms1 = n_ms2 = 0; windows = set()
    ctx = etree.iterparse(a.mzml, events=("end",), tag=NS+"spectrum")
    for _, spec in ctx:
        ms_level, rt, lo, hi, mz, inten = spectrum_info(spec)
        if mz is not None and rt is not None:
            if ms_level == 2 and lo is not None:
                n_ms2 += 1; windows.add((round(lo,2), round(hi,2)))
                for seq, z, pmz, frags, _iso in targets:
                    if lo <= pmz <= hi:
                        for f, fmz in frags:
                            rows.append((a.sample, seq, z, f, round(fmz,4), rt, window_sum(mz, inten, fmz, a.ppm)))
            elif ms_level == 1:
                n_ms1 += 1
                for seq, z, pmz, _frags, iso in targets:
                    for f, fmz in iso:
                        rows.append((a.sample, seq, z, f, round(fmz,4), rt, window_sum(mz, inten, fmz, a.ppm)))
        spec.clear()
        while spec.getprevious() is not None: del spec.getparent()[0]
    out = os.path.join(a.out, f"xic_{a.sample}.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample","peptide","charge","fragment","fragment_mz","rt_min","intensity"]); w.writerows(rows)
    with open(os.path.join(a.out, "target_table.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["peptide","charge","precursor_mz","class","fragments"])
        for (seq, z, cls), (_s, _z, pmz, frags, _i) in zip(TARGETS, targets):
            w.writerow([seq, z, round(pmz,4), cls, ";".join(f"{l}={m:.4f}" for l, m in frags)])
    sys.stderr.write(f"[{a.sample}] MS1: {n_ms1}; MS2: {n_ms2}; windows: {len(windows)}; rows: {len(rows)} -> {out}\n")

if __name__ == "__main__":
    main()
