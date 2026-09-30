#!/usr/bin/env python3
"""Compare owned-row dumps (MARS_ROW_DUMP) of two runs by node identity (SFC key).

    compare_rows.py <ref-prefix> <test-prefix> [--rtol 1e-10]

Each run wrote <prefix>.rank<r>.txt with lines "sfc_key diag abs_row_sum rhs". Checks that every
node is owned by exactly one rank in each run, that both runs own the same nodes, and that each
row's diagonal, absolute row sum and RHS agree. Exit code 1 lists the offending nodes.
"""
import argparse
import glob
import sys


def load(prefix):
    files = sorted(glob.glob(prefix + ".rank*.txt"))
    if not files:
        raise SystemExit(f"FAIL: no dump files {prefix}.rank*.txt")
    rows, dup = {}, []
    for fn in files:
        with open(fn) as f:
            for line in f:
                k, d, a, r = line.split()
                k = int(k)
                if k in rows:
                    dup.append(k)
                rows[k] = (float(d), float(a), float(r))
    return rows, dup, len(files)


def compare(ref_prefix, test_prefix, rtol):
    ref, ref_dup, _ = load(ref_prefix)
    test, test_dup, nfiles = load(test_prefix)
    problems = []
    for name, dup in (("reference", ref_dup), ("test", test_dup)):
        if dup:
            problems.append(f"{len(dup)} node(s) owned by more than one rank in the {name} run, e.g. {dup[:5]}")
    missing = sorted(set(ref) - set(test))
    extra = sorted(set(test) - set(ref))
    if missing:
        problems.append(f"{len(missing)} node(s) owned by no rank in the test run, e.g. {missing[:5]}")
    if extra:
        problems.append(f"{len(extra)} node(s) in the test run absent from the reference, e.g. {extra[:5]}")

    # Relative to the row's own size, with a floor so rows that are ~0 in both runs pass.
    scale = max((max(abs(v) for v in t) for t in ref.values()), default=1.0)
    bad = []
    for k in sorted(set(ref) & set(test)):
        for i, what in enumerate(("diag", "abs_row_sum", "rhs")):
            a, b = ref[k][i], test[k][i]
            denom = max(abs(a), abs(b), 1e-12 * scale)
            if abs(a - b) / denom > rtol:
                bad.append((k, what, a, b))
    if bad:
        keys = sorted({k for k, *_ in bad})
        problems.append(f"{len(keys)} node row(s) differ (rtol {rtol:.0e}); first ones:")
        problems += [f"  key {k} {w}: ref {a:.17g} test {b:.17g}" for k, w, a, b in bad[:10]]
    print(f"rows: reference {len(ref)}, test {len(test)} across {nfiles} rank file(s)")
    return problems


def main():
    p = argparse.ArgumentParser()
    p.add_argument("ref")
    p.add_argument("test")
    p.add_argument("--rtol", type=float, default=1e-10)
    a = p.parse_args()
    problems = compare(a.ref, a.test, a.rtol)
    if problems:
        print("FAIL: " + "\n".join(problems))
        sys.exit(1)
    print("PASS: every node owned once, identical rows")


if __name__ == "__main__":
    main()
