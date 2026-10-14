#!/usr/bin/env python3
"""Fit the 07_rarefaction.sh curves and project the sequencing each sample needs.

Three steps, all from files the pipeline has already written:

1. Library complexity. For each sample, unique mitochondrial reads U are
   modelled against mitochondrial reads sequenced N with the gamma-Poisson
   (negative binomial) saturation curve

       U(N) = C * (1 - (1 + N / (k * C)) ** -k)

   where C is the number of distinct mitochondrial molecules in the library and
   k describes how unevenly they are amplified (k -> infinity is the plain
   Poisson case; small k means PCR jackpotting, so new unique reads appear more
   slowly than the Poisson curve predicts). C caps what any amount of further
   sequencing of the same library can deliver.

2. Calibration. The MitoFinder outcome of 06 (complete: all 13 protein-coding
   genes and both rRNAs) is regressed on log10 unique COX1 depth by logistic
   regression, using only samples whose depth was measured against their own
   mitogenome or COI. The depth at which P(complete) reaches TARGET_P is the
   target.

3. Projection. The target depth is converted to unique reads, inverted
   through the fitted curve to mitochondrial reads, divided by the library's
   mitochondrial fraction to give trimmed reads, and scaled by the sample's
   own fastp retention to raw read pairs, which is the unit sequencing is
   ordered in.

Writes, under OUTDIR (default results/07_rarefaction):
    projection.tsv     one row per sample
    calibration.tsv    logistic fit and target depths
    report_data.json   everything the HTML report plots

Usage, from the project root:
    python lib/rarefaction_model.py [--target-p 0.9] [--outdir DIR]
"""
import argparse
import csv
import json
import math
import os
import sys

import numpy as np

# Below this, the COI that reads were mapped to is too short for its depth to
# be trusted: reads overhang both ends and are clipped or lost.
MIN_COI_BP = 500

# Beyond this multiple of the sequencing a library already has, the projection
# runs far out along a flattening curve, and a new library is the realistic
# answer whatever the model says.
MAX_FOLD = 10


def verdict(d, target):
    """One of: complete, enough depth but failed, sequence more, new library,
    no mitochondrial reads. Mirrored in the HTML report's project()."""
    if d["outcome"] == "complete":
        return "complete"
    if not math.isfinite(d.get("raw_pairs_needed", math.inf)) and d["max_depth"] == 0:
        return "no mitochondrial reads"
    if d["mean_depth"] >= target:
        return "enough depth but failed"
    need = d["raw_pairs_needed"]
    if not math.isfinite(need) or need > MAX_FOLD * d["raw_reads"] / 2:
        return "new library needed"
    return "sequence more"


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


# ----------------------------------------------------------------------------
# Saturation model
# ----------------------------------------------------------------------------

def gp_unique(n, c, k):
    """Expected unique reads after n reads from c molecules, dispersion k."""
    n = np.asarray(n, dtype=float)
    return c * (1.0 - np.power(1.0 + n / (k * c), -k))


def gp_reads_for(u, c, k):
    """Reads needed to reach u unique reads; inf if u >= c."""
    if u >= c:
        return math.inf
    return k * c * ((1.0 - u / c) ** (-1.0 / k) - 1.0)


def fit_saturation(n, u):
    """Least-squares fit of (C, k), weighting each point by 1/max(u, 1).

    Returns (C, k, saturated). saturated is False when the data show too few
    duplicates to bound C, in which case C is a lower bound.
    """
    n = np.asarray(n, dtype=float)
    u = np.asarray(u, dtype=float)
    keep = n > 0
    n, u = n[keep], u[keep]
    if len(n) == 0 or u.max() <= 0:
        return float("nan"), float("nan"), False
    w = 1.0 / np.maximum(u, 1.0)
    umax = u.max()
    cap = umax * 1e3
    best = (math.inf, umax, 1.0)
    for k in np.logspace(-1, 3, 41):
        # C is searched on a log grid, then refined around the best value.
        cs = np.logspace(math.log10(umax), math.log10(cap), 200)
        sse = [float(np.sum(w * (gp_unique(n, c, k) - u) ** 2)) for c in cs]
        j = int(np.argmin(sse))
        lo, hi = cs[max(j - 1, 0)], cs[min(j + 1, len(cs) - 1)]
        for c in np.linspace(lo, hi, 50):
            s = float(np.sum(w * (gp_unique(n, c, k) - u) ** 2))
            if s < best[0]:
                best = (s, c, k)
    _, c, k = best
    saturated = c < cap * 0.5
    return float(c), float(k), saturated


# ----------------------------------------------------------------------------
# Logistic calibration
# ----------------------------------------------------------------------------

def logistic_fit(x, y, iters=100):
    """Newton-Raphson logistic regression of y (0/1) on x. Returns (b0, b1)."""
    X = np.column_stack([np.ones_like(x), x])
    beta = np.zeros(2)
    for _ in range(iters):
        p = 1.0 / (1.0 + np.exp(-X @ beta))
        W = p * (1 - p)
        H = X.T @ (X * W[:, None]) + 1e-9 * np.eye(2)
        step = np.linalg.solve(H, X.T @ (y - p))
        beta += step
        if np.max(np.abs(step)) < 1e-10:
            break
    return float(beta[0]), float(beta[1])


def depth_at_p(b0, b1, p):
    return 10 ** ((math.log(p / (1 - p)) - b0) / b1)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", default="results/07_rarefaction")
    ap.add_argument("--trimdir", default="results/01_trim")
    ap.add_argument("--baitdir", default="results/02_bait")
    ap.add_argument("--mitodir", default="results/06_mitogenome")
    ap.add_argument("--finaldir", default="results/coi_final")
    ap.add_argument("--target-p", type=float, default=0.9,
                    help="probability of a complete mitogenome to plan for")
    args = ap.parse_args()

    rare = read_tsv(os.path.join(args.outdir, "rarefaction.tsv"))
    bait = {r["sample"]: r for r in read_tsv(os.path.join(args.baitdir, "bait_stats.tsv"))}
    mito = {r["sample"]: r for r in read_tsv(os.path.join(args.mitodir, "mitogenome_stats.tsv"))}

    curves = {}
    for r in rare:
        curves.setdefault(r["sample"], []).append(r)

    samples = []
    for s, rows in curves.items():
        rows.sort(key=lambda r: num(r["fraction"]))
        full = rows[-1]
        if num(full["fraction"]) != 1:
            sys.exit(f"{s}: rarefaction has no fraction 1 row")

        # Raw and trimmed read counts from fastp, both counting R1 and R2.
        with open(os.path.join(args.trimdir, f"{s}.fastp.json")) as fh:
            fp = json.load(fh)["summary"]
        raw_reads = fp["before_filtering"]["total_reads"]
        raw_bases = fp["before_filtering"]["total_bases"]
        trimmed = num(bait[s]["total_reads"])

        coi_bp = None
        coi = os.path.join(args.finaldir, f"{s}_COI.fasta")
        if os.path.exists(coi):
            with open(coi) as fh:
                coi_bp = sum(len(l.strip()) for l in fh if not l.startswith(">"))

        g = int(num(mito.get(s, {}).get("genes_annotated", 0)) or 0)
        outcome = "complete" if g == 15 else ("partial" if g > 0 else "none")

        n = [num(r["mito_reads"]) for r in rows]
        u = [num(r["mito_unique"]) for r in rows]
        frac = [num(r["fraction"]) for r in rows]
        c, k, saturated = fit_saturation(n, u)

        mito_reads = n[-1]
        uniq = u[-1]
        mean_d = num(full["cox1_mean_depth"])
        samples.append(dict(
            sample=s,
            ref_source=full["ref_source"],
            outcome=outcome,
            genes_annotated=g,
            coi_bp=coi_bp,
            raw_reads=raw_reads,
            raw_bases=raw_bases,
            trimmed_reads=trimmed,
            retention=trimmed / raw_reads if raw_reads else float("nan"),
            mito_reads=mito_reads,
            mito_unique=uniq,
            mito_fraction=mito_reads / trimmed if trimmed else 0.0,
            dup_rate=1 - uniq / mito_reads if mito_reads else float("nan"),
            mean_depth=mean_d,
            median_depth=num(full["cox1_median_depth"]),
            depth_per_unique=mean_d / uniq if uniq else float("nan"),
            ref_mean_depth=num(bait[s]["cox1_mean_depth"]),
            C=c, k=k, saturated=saturated,
            points=[dict(fraction=f, raw_pairs=f * raw_reads / 2, mito_reads=a, unique=b,
                         mean_depth=num(r["cox1_mean_depth"]))
                    for f, a, b, r in zip(frac, n, u, rows)],
        ))

    # ---- calibration -------------------------------------------------------
    # Only depths measured against the sample's own sequence are unbiased, so
    # only they inform the threshold: its own mitogenome, or its own COI when
    # that is long enough.
    def in_calibration(d):
        if d["mean_depth"] <= 0:
            return False
        if d["ref_source"] == "own_mitogenome":
            return True
        return d["ref_source"] == "own_coi" and (d["coi_bp"] or 0) >= MIN_COI_BP

    for d in samples:
        d["in_calibration"] = in_calibration(d)
    cal = [d for d in samples if d["in_calibration"]]
    x = np.array([math.log10(d["mean_depth"]) for d in cal])
    y = np.array([1.0 if d["outcome"] == "complete" else 0.0 for d in cal])
    b0, b1 = logistic_fit(x, y)
    y_any = np.array([1.0 if d["outcome"] != "none" else 0.0 for d in cal])
    a0, a1 = logistic_fit(x, y_any)
    target = depth_at_p(b0, b1, args.target_p)
    d50 = depth_at_p(b0, b1, 0.5)

    # ---- projection --------------------------------------------------------
    for d in samples:
        d["target_depth"] = target
        r = d["depth_per_unique"]
        if not (d["mito_unique"] > 0 and r > 0 and d["mito_fraction"] > 0):
            d.update(max_depth=0.0, mito_reads_needed=math.inf, raw_pairs_needed=math.inf)
        else:
            d["max_depth"] = d["C"] * r
            n_needed = gp_reads_for(target / r, d["C"], d["k"])
            d["mito_reads_needed"] = n_needed
            d["raw_pairs_needed"] = n_needed / d["mito_fraction"] / d["retention"] / 2
        d["verdict"] = verdict(d, target)
        if d["ref_source"] == "reference":
            d["verdict"] += " (lower-bound depth)"
        elif d["ref_source"] == "own_coi" and (d["coi_bp"] or 0) < MIN_COI_BP:
            d["verdict"] += " (short COI)"

    raw_pairs_now = lambda d: d["raw_reads"] / 2
    for d in samples:
        need = d["raw_pairs_needed"]
        d["raw_pairs_now"] = raw_pairs_now(d)
        d["additional_raw_pairs"] = max(need - d["raw_pairs_now"], 0) if math.isfinite(need) else math.inf
        d["fold"] = need / d["raw_pairs_now"] if math.isfinite(need) else math.inf
        d["bp_per_pair"] = d["raw_bases"] / d["raw_pairs_now"]
        d["additional_raw_gbp"] = (d["additional_raw_pairs"] * d["bp_per_pair"] / 1e9
                                   if math.isfinite(d["additional_raw_pairs"]) else math.inf)
        # Extrapolated curve, in raw read pairs, out to 20x the current effort
        # or the point the target is reached, whichever is further.
        xmax = d["raw_pairs_now"] * 20
        if math.isfinite(need):
            xmax = max(xmax, need * 1.2)
        xs = np.unique(np.concatenate([np.logspace(4, math.log10(max(xmax, 1e5)), 60),
                                       [d["raw_pairs_now"]]]))
        if d["mito_unique"] > 0 and math.isfinite(d["C"]):
            per_pair = 2 * d["retention"] * d["mito_fraction"]
            ys = gp_unique(xs * per_pair, d["C"], d["k"]) * d["depth_per_unique"]
        else:
            ys = np.zeros_like(xs)
        d["curve"] = [dict(raw_pairs=float(a), mean_depth=float(b)) for a, b in zip(xs, ys)]

    # ---- write -------------------------------------------------------------
    order = {"complete": 0, "partial": 1, "none": 2}
    samples.sort(key=lambda d: (order[d["outcome"]], d["sample"]))

    def fmt(v, p=4):
        if isinstance(v, float):
            if math.isinf(v):
                return "Inf"
            if math.isnan(v):
                return "NA"
            return f"{v:.{p}g}"
        return "NA" if v is None else str(v)

    cols = ["sample", "ref_source", "outcome", "genes_annotated", "coi_bp",
            "raw_pairs_now", "mito_fraction", "dup_rate", "mito_unique",
            "mean_depth", "median_depth", "ref_mean_depth", "C", "k", "saturated",
            "max_depth", "target_depth", "raw_pairs_needed", "additional_raw_pairs",
            "additional_raw_gbp", "fold", "verdict"]
    with open(os.path.join(args.outdir, "projection.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for d in samples:
            fh.write("\t".join(fmt(d[c]) for c in cols) + "\n")

    with open(os.path.join(args.outdir, "calibration.tsv"), "w") as fh:
        fh.write("model\tintercept\tslope_per_log10_depth\tn\tn_positive\tdepth_p50\tdepth_target\ttarget_p\n")
        fh.write(f"complete\t{b0:.4f}\t{b1:.4f}\t{len(cal)}\t{int(y.sum())}\t{d50:.2f}\t{target:.2f}\t{args.target_p}\n")
        fh.write(f"any_genes\t{a0:.4f}\t{a1:.4f}\t{len(cal)}\t{int(y_any.sum())}\t"
                 f"{depth_at_p(a0, a1, 0.5):.2f}\t{depth_at_p(a0, a1, args.target_p):.2f}\t{args.target_p}\n")

    def clean(v):
        if isinstance(v, float) and not math.isfinite(v):
            return None
        if isinstance(v, (np.floating,)):
            return clean(float(v))
        if isinstance(v, dict):
            return {kk: clean(vv) for kk, vv in v.items()}
        if isinstance(v, list):
            return [clean(vv) for vv in v]
        if isinstance(v, (np.bool_,)):
            return bool(v)
        return v

    report = dict(
        calibration=dict(complete=dict(b0=b0, b1=b1), any_genes=dict(b0=a0, b1=a1),
                         target_p=args.target_p, target_depth=target, depth_p50=d50,
                         n=len(cal), min_coi_bp=MIN_COI_BP, max_fold=MAX_FOLD),
        samples=samples,
    )
    with open(os.path.join(args.outdir, "report_data.json"), "w") as fh:
        json.dump(clean(report), fh, separators=(",", ":"))

    # ---- console summary ---------------------------------------------------
    print(f"calibration on {len(cal)} samples measured against their own sequence:")
    print(f"  P(complete) = 0.5 at {d50:.1f}x, {args.target_p} at {target:.1f}x unique mean COX1 depth")
    counts = {}
    for d in samples:
        counts[d["verdict"]] = counts.get(d["verdict"], 0) + 1
    for kk in sorted(counts):
        print(f"  {kk}: {counts[kk]}")
    print(f"wrote {args.outdir}/projection.tsv, calibration.tsv, report_data.json")


if __name__ == "__main__":
    main()
