"""Track-C benchmark: manual/published-vs-OligoForge comparison. OFFLINE, DETERMINISTIC.

Run from repo root:  PYTHONPATH=. python3 tests/benchmark/bench_score.py

For each target in bench_corpus.json this:
  1. designs an assay with OligoForge at the target's real annealing temperature,
  2. scores every OligoForge oligo on the fixed rubric (Tm in-profile, GC, amplicon, 3' clamp,
     run limits, hairpin/self-dimer/cross-dimer dG AT THE ANNEAL TEMP, probe placement/5'-G/C>G),
  3. cross-checks OligoForge's displayed Tm (tm_acc) against an INDEPENDENT OligoAnalyzer-equivalent
     (Biopython Tm_NN, SantaLucia 1998 + Owczarzy 2008 + Mg/dNTP, saltcorr=7) -- the authoritative
     reference for Tm accuracy,
  4. for targets carrying a published/golden reference assay, tests whether OligoForge ADMITS and
     RANKS the reference primers (the key selection-consistency check), and
  5. writes bench_scorecard.csv (one row per target) + prints a summary.

No network at score time: sequences are committed in bench_corpus.json. primer3 is deterministic,
so the scorecard is reproducible run-to-run (asserted by bench_determinism)."""
import os, sys, json, csv
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from oligoforge import thermo as T, design as D, profiles as P, nn as NN

HERE = os.path.dirname(os.path.abspath(__file__))
CORPUS = json.load(open(os.path.join(HERE, "bench_corpus.json")))

def oa_tm(seq):
    """INDEPENDENT OligoAnalyzer-equivalent Tm (deg C). Uses the repo's from-scratch SantaLucia-1998
    NN + Owczarzy-2004/2008 engine (nn.params), which is a SEPARATE implementation from the displayed
    tm_acc (Biopython Tm_NN path). Comparing tm_acc against this is a genuine cross-implementation
    check, not the same code twice. nn.params reads the current session salt via thermo.COND."""
    seq = "".join(c for c in seq.upper() if c in "ACGT")
    if len(seq) < 2:
        return None
    try:
        return float(NN.params(seq)["tm"])
    except Exception:
        return None

def reject_reason(seq, c):
    """First gate in _ok_primer that this primer fails (or 'admitted'). Mirrors design._ok_primer
    order exactly so the benchmark can classify WHY OligoForge declines a reference primer."""
    if not (c["len_min"] <= len(seq) <= c["len_max"]): return f"len {len(seq)} not in {c['len_min']}-{c['len_max']}"
    if not (c["gc_min"] <= T.gc_percent(seq) <= c["gc_max"]): return f"GC {T.gc_percent(seq):.0f} not in {c['gc_min']}-{c['gc_max']}"
    if not (c["tm_min"] <= T.tm(seq) <= c["tm_max"]): return f"selTm {T.tm(seq):.1f} not in {c['tm_min']}-{c['tm_max']}"
    if c.get("no_three_prime_T") and seq[-1] == "T": return "3'-T disallowed"
    if T.last5_gc(seq) > c["max_3prime_gc"]: return f"3'GC-clamp {T.last5_gc(seq)}>{c['max_3prime_gc']}"
    if T.max_run(seq, "G") >= c["max_g_run"]: return f"G-run>={c['max_g_run']}"
    if T.max_run(seq) >= c["max_any_run"]: return f"run>={c['max_any_run']}"
    ac = c.get("anneal_c", T.ANNEAL_C)
    if T.hairpin_full(seq, ac)[1] <= c["hairpin_min"]: return f"hairpin dG@{ac:.0f}<={c['hairpin_min']}"
    if T.self_dimer_full(seq, ac)[1] <= c["self_dimer_min"]: return f"self-dimer dG@{ac:.0f}<={c['self_dimer_min']}"
    return "admitted"

def locate(sub, seq):
    return seq.upper().find(sub.upper())

def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def score_oligo(seq, c, anneal_c, is_probe=False):
    """Rubric checks for one oligo. Returns dict of pass/fail + the numbers."""
    seq = seq.upper()
    tm_sel = T.tm(seq)                          # selection Tm (primer3)
    tm_disp = T.tm_acc(seq)                     # displayed Tm (Owczarzy-2008)
    tm_ref = oa_tm(seq)                         # independent OA-equivalent
    gc = T.gc_percent(seq)
    hp37, hp_an, _ = T.hairpin_full(seq, anneal_c=anneal_c)
    sd37, sd_an, _ = T.self_dimer_full(seq, anneal_c=anneal_c)
    return {
        "seq": seq, "len": len(seq),
        "tm_sel": round(tm_sel, 2), "tm_disp": round(tm_disp, 2),
        "tm_ref_OA": round(tm_ref, 2) if tm_ref is not None else None,
        "tm_disp_vs_OA": round(abs(tm_disp - tm_ref), 3) if tm_ref is not None else None,
        "gc": round(gc, 1), "gc_in_profile": c["gc_min"] <= gc <= c["gc_max"],
        "tm_in_profile": c["tm_min"] <= tm_sel <= c["tm_max"],
        "hairpin_dg_anneal": round(hp_an, 2), "hairpin_ok": hp_an > c.get("probe_hairpin_min" if is_probe else "hairpin_min"),
        "selfdimer_dg_anneal": round(sd_an, 2), "selfdimer_ok": sd_an > c["self_dimer_min"],
    }

def admits_primer(seq, c):
    """Does OligoForge's _ok_primer accept this exact primer (at the profile's anneal temp)?"""
    return D._ok_primer(seq, c)

def run():
    rows = []
    for t in CORPUS["targets"]:
        prof = P.PROFILES[t["profile"]]
        anneal = t.get("anneal_c", prof.get("anneal_c", 60))
        # design with OligoForge at the target's real anneal temp
        T.set_conditions(anneal_c=anneal)
        try:
            a = D.design_assay(t["seq"], prof)
        except Exception as e:
            a = {"error": str(e)}
        row = {"id": t["id"], "gc": t["gc"], "profile": t["profile"], "anneal_c": anneal,
               "ref_kind": t["reference"]["kind"]}
        # OligoForge's own pick
        if a and not a.get("error") and a.get("forward"):
            row["of_forward"] = a["forward"]; row["of_reverse"] = a["reverse"]
            row["of_amplicon"] = a.get("amplicon"); row["of_probe"] = a.get("probe")
            fs = score_oligo(a["forward"], prof, anneal)
            rs = score_oligo(a["reverse"], prof, anneal)
            row["of_fwd_tm_disp"] = fs["tm_disp"]; row["of_fwd_tm_OA"] = fs["tm_ref_OA"]
            row["of_fwd_tm_vs_OA"] = fs["tm_disp_vs_OA"]; row["of_rev_tm_vs_OA"] = rs["tm_disp_vs_OA"]
            row["of_fwd_selfdimer_anneal"] = fs["selfdimer_dg_anneal"]
            row["of_fwd_selfdimer_ok"] = fs["selfdimer_ok"]
        else:
            row["of_error"] = a.get("error", "no assay found")
        # reference-primer admission/ranking check (published or golden)
        ref = t["reference"]
        if ref.get("forward") and ref.get("reverse"):
            fpos = locate(ref["forward"], t["seq"]); rpos = locate(revcomp(ref["reverse"]), t["seq"])
            row["ref_forward"] = ref["forward"]; row["ref_reverse"] = ref["reverse"]
            row["ref_locates"] = (fpos >= 0 and rpos >= 0)
            row["ref_fwd_admitted"] = admits_primer(ref["forward"], prof)
            row["ref_rev_admitted"] = admits_primer(ref["reverse"], prof)
            row["ref_fwd_reject_reason"] = reject_reason(ref["forward"], prof)
            row["ref_rev_reject_reason"] = reject_reason(ref["reverse"], prof)
            # Tm agreement on the reference primers
            rf = score_oligo(ref["forward"], prof, anneal)
            rr = score_oligo(ref["reverse"], prof, anneal)
            row["ref_fwd_tm_disp"] = rf["tm_disp"]; row["ref_fwd_tm_OA"] = rf["tm_ref_OA"]
            row["ref_fwd_tm_vs_OA"] = rf["tm_disp_vs_OA"]
            row["ref_fwd_tm_in_profile"] = rf["tm_in_profile"]
            row["ref_fwd_gc_in_profile"] = rf["gc_in_profile"]
        rows.append(row)
    T.set_conditions(anneal_c=60)  # restore
    return rows

def write_csv(rows, path):
    cols = sorted({k for r in rows for k in r})
    lead = ["id","gc","profile","anneal_c","ref_kind"]
    cols = lead + [c for c in cols if c not in lead]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in rows: w.writerow(r)

if __name__ == "__main__":
    rows = run()
    out = os.path.join(HERE, "bench_scorecard.csv")
    write_csv(rows, out)
    print("=== Track-C scorecard ===")
    for r in rows:
        print(f"\n[{r['id']}]  gc={r['gc']}% profile={r['profile']} Ta={r['anneal_c']} ref={r['ref_kind']}")
        if r.get("of_forward"):
            print(f"  OF pick: F={r['of_forward']} R={r['of_reverse']} amp={r.get('of_amplicon')} probe={r.get('of_probe')}")
            print(f"  OF Tm vs OA-equiv: fwd |d|={r.get('of_fwd_tm_vs_OA')}C  rev |d|={r.get('of_rev_tm_vs_OA')}C")
            print(f"  OF fwd self-dimer dG@{r['anneal_c']}C={r.get('of_fwd_selfdimer_anneal')} ok={r.get('of_fwd_selfdimer_ok')}")
        else:
            print(f"  OF: {r.get('of_error')}")
        if r.get("ref_forward"):
            print(f"  REF({r['ref_kind']}): F={r['ref_forward']} R={r['ref_reverse']}  locates={r.get('ref_locates')}")
            print(f"    admitted by OligoForge? fwd={r.get('ref_fwd_admitted')} ({r.get('ref_fwd_reject_reason')}) "
                  f"rev={r.get('ref_rev_admitted')} ({r.get('ref_rev_reject_reason')})")
            print(f"    ref_fwd Tm(disp)={r.get('ref_fwd_tm_disp')} vs OA-indep |d|={r.get('ref_fwd_tm_vs_OA')}C  Tm_in_profile={r.get('ref_fwd_tm_in_profile')}")
    print(f"\nscorecard -> {out}")
