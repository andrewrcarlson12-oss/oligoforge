"""Vendor / chemistry constraint profiles.

Each profile is the rule-set a design must satisfy for that vendor+chemistry.
Stage 2 expands this into a full per-vendor linter (IDT PrimeTime / Affinity Plus,
Thermo TaqMan, Bio-Rad, generic SYBR). For now it carries the IDT TaqMan profile
the engine designs against.
"""

IDT_TAQMAN = dict(
    name="IDT PrimeTime (hydrolysis probe)",
    # primers
    len_min=18, len_max=24, gc_min=35.0, gc_max=65.0,
    tm_min=59.0, tm_max=64.5, tm_opt=62.0,
    no_three_prime_T=True, max_3prime_gc=3, min_3prime_gc=1, max_g_run=4, max_any_run=5,
    hairpin_min=-2.0, self_dimer_min=-6.0, pair_dimer_min=-6.0,
    pair_tm_gap_max=2.0,
    # amplicon
    amp_min=70, amp_max=150, min_probe_gap=24,
    # probe
    probe_len_min=18, probe_len_max=30,
    probe_offset_min=5.0, probe_offset_max=10.5,
    probe_hairpin_min=-1.5,
)


# ---- additional chemistry/vendor profiles (guideline defaults; edit to your kit) ----
IDT_AFFINITY = dict(
    name="IDT Affinity Plus (LNA probe)",
    len_min=18, len_max=24, gc_min=35.0, gc_max=65.0,
    tm_min=59.0, tm_max=64.5, tm_opt=62.0,
    no_three_prime_T=True, max_3prime_gc=3, min_3prime_gc=1, max_g_run=4, max_any_run=5,
    hairpin_min=-2.0, self_dimer_min=-6.0, pair_dimer_min=-6.0, pair_tm_gap_max=2.0,
    amp_min=60, amp_max=150, min_probe_gap=20,
    probe_len_min=10, probe_len_max=25,        # short LNA core
    probe_offset_min=8.0, probe_offset_max=14.0,  # LNA raises effective Tm
    probe_hairpin_min=-1.5,
    notes="LNA probe: <=6 LNA bases, <=4 sequential; 5'FAM/3'IBFQ, no ZEN; amplicon <=150.",
)
THERMO_TAQMAN = dict(
    name="Thermo Fisher TaqMan (MGB)",
    len_min=15, len_max=30, gc_min=30.0, gc_max=80.0,
    tm_min=58.0, tm_max=60.5, tm_opt=59.0,       # ABI designs primers cooler
    no_three_prime_T=False, max_3prime_gc=3, min_3prime_gc=1, max_g_run=4, max_any_run=5,
    hairpin_min=-2.0, self_dimer_min=-6.0, pair_dimer_min=-6.0, pair_tm_gap_max=2.0,
    amp_min=50, amp_max=150, min_probe_gap=10,
    probe_len_min=13, probe_len_max=25,          # MGB probes are short
    probe_offset_min=6.0, probe_offset_max=12.0,
    probe_hairpin_min=-1.5,
    notes="MGB probe raises effective Tm; keep amplicon <100 where possible; 5' no G; more C than G.",
)
BIORAD_PROBE = dict(
    name="Bio-Rad PrimePCR / iTaq (probe)",
    len_min=18, len_max=25, gc_min=40.0, gc_max=60.0,
    tm_min=60.0, tm_max=64.0, tm_opt=62.0,
    no_three_prime_T=True, max_3prime_gc=3, min_3prime_gc=1, max_g_run=4, max_any_run=5,
    hairpin_min=-2.0, self_dimer_min=-6.0, pair_dimer_min=-6.0, pair_tm_gap_max=2.0,
    amp_min=70, amp_max=200, min_probe_gap=20,
    probe_len_min=18, probe_len_max=30, probe_offset_min=5.0, probe_offset_max=10.0,
    probe_hairpin_min=-1.5,
)
SYBR_GENERIC = dict(
    name="Generic SYBR Green (no probe)",
    len_min=18, len_max=24, gc_min=40.0, gc_max=60.0,
    tm_min=58.0, tm_max=62.0, tm_opt=60.0,
    no_three_prime_T=True, max_3prime_gc=3, min_3prime_gc=1, max_g_run=4, max_any_run=5,
    hairpin_min=-2.0, self_dimer_min=-6.0, pair_dimer_min=-6.0, pair_tm_gap_max=1.5,
    amp_min=70, amp_max=200, min_probe_gap=0,
    probe_len_min=0, probe_len_max=0, probe_offset_min=0, probe_offset_max=0,
    probe_hairpin_min=-99,
    no_probe=True,
    notes="No probe. Primer-dimers are SYBR-critical: F/R, F-self, R-self must all clear the dimer floor.",
)

PARASITE_MTDNA = dict(
    name="AT-rich parasite mtDNA (low-Tm TaqMan, ~54C)",
    len_min=18, len_max=28, gc_min=25.0, gc_max=60.0,
    tm_min=52.0, tm_max=58.0, tm_opt=55.0,
    # AT-rich: allow 6-mer poly-A/T runs (G-run still capped at 4)
    no_three_prime_T=False, max_3prime_gc=3, max_g_run=4, max_any_run=7,
    hairpin_min=-2.0, self_dimer_min=-6.0, pair_dimer_min=-6.0, pair_tm_gap_max=2.5,
    amp_min=70, amp_max=150, min_probe_gap=10,
    probe_len_min=12, probe_len_max=28, probe_offset_min=6.0, probe_offset_max=12.0,
    probe_hairpin_min=-1.5,
    notes="AT-rich mitochondrial targets (e.g. apicomplexan cytb), ~54C anneal. "
          "An all-DNA probe often can't reach Tm here; order the probe as Affinity Plus (LNA).",
)
PARASITE_SYBR = dict(
    name="AT-rich parasite mtDNA (low-Tm SYBR, ~54C)",
    len_min=18, len_max=28, gc_min=25.0, gc_max=60.0,
    tm_min=52.0, tm_max=58.0, tm_opt=55.0,
    # AT-rich: allow 6-mer poly-A/T runs (G-run still capped at 4)
    no_three_prime_T=False, max_3prime_gc=3, max_g_run=4, max_any_run=7,
    hairpin_min=-2.0, self_dimer_min=-6.0, pair_dimer_min=-6.0, pair_tm_gap_max=2.5,
    amp_min=70, amp_max=250, min_probe_gap=0,
    probe_len_min=0, probe_len_max=0, probe_offset_min=0, probe_offset_max=0,
    probe_hairpin_min=-99, no_probe=True,
    notes="No probe; AT-rich mtDNA primers only (e.g. apicomplexan cytb), ~54C anneal.",
)

GC_RICH = dict(
    name="GC-rich target (high-Tm TaqMan)",
    len_min=16, len_max=24, gc_min=50.0, gc_max=85.0,
    tm_min=60.0, tm_max=68.0, tm_opt=64.0,
    # GC-rich: short primers avoid runaway Tm; tolerate higher GC% and short G/C runs
    no_three_prime_T=True, max_3prime_gc=4, min_3prime_gc=1, max_g_run=5, max_any_run=6,
    hairpin_min=-3.0, self_dimer_min=-8.0, pair_dimer_min=-8.0, pair_tm_gap_max=2.5,
    amp_min=70, amp_max=150, min_probe_gap=10,
    probe_len_min=15, probe_len_max=30, probe_offset_min=5.0, probe_offset_max=10.5,
    probe_hairpin_min=-2.5,
    notes="GC-rich targets (e.g. high-GC Actinobacteria like Mycobacterium/Streptomyces, "
          "GC-rich viral genomes). Short primers keep Tm in range; higher GC% and short G/C "
          "runs are tolerated. These templates fold strongly \u2014 check the amplicon structure "
          "and consider DMSO or 7-deaza-dGTP in the reaction.",
)
GC_RICH_SYBR = dict(
    name="GC-rich target (high-Tm SYBR, no probe)",
    len_min=16, len_max=24, gc_min=50.0, gc_max=85.0,
    tm_min=60.0, tm_max=68.0, tm_opt=64.0,
    no_three_prime_T=True, max_3prime_gc=4, min_3prime_gc=1, max_g_run=5, max_any_run=6,
    hairpin_min=-3.0, self_dimer_min=-8.0, pair_dimer_min=-8.0, pair_tm_gap_max=2.5,
    amp_min=70, amp_max=200, min_probe_gap=0,
    probe_len_min=0, probe_len_max=0, probe_offset_min=0, probe_offset_max=0,
    probe_hairpin_min=-99, no_probe=True,
    notes="No probe; GC-rich primers only. Strong template structure \u2014 consider DMSO / 7-deaza-dGTP.",
)

PROFILES = {
    "idt_taqman": IDT_TAQMAN, "idt_affinity": IDT_AFFINITY,
    "thermo_taqman": THERMO_TAQMAN, "biorad_probe": BIORAD_PROBE,
    "sybr_generic": SYBR_GENERIC,
    "parasite_mtdna": PARASITE_MTDNA, "parasite_sybr": PARASITE_SYBR,
    "gc_rich": GC_RICH, "gc_rich_sybr": GC_RICH_SYBR,
}


def lint_oligo(seq, role, profile):
    """Return [(rule, status, detail)] for one oligo against a profile.
    role: 'forward' | 'reverse' | 'primer' | 'probe'."""
    from . import thermo as T
    seq = seq.upper()
    out = []
    def chk(name, ok, detail):
        out.append((name, "PASS" if ok else "FAIL", detail))
    def warn(name, ok, detail):
        out.append((name, "PASS" if ok else "WARN", detail))

    _deg = T.has_degenerate(seq)
    warn("no degenerate bases", not _deg,
         "contains IUPAC degeneracy (N/R/Y…); Tm & structure values are approximate" if _deg else "none")

    if role == "probe":
        chk("length", profile["probe_len_min"] <= len(seq) <= profile["probe_len_max"],
            f"{len(seq)} nt (allowed {profile['probe_len_min']}-{profile['probe_len_max']})")
        chk("5' not G", seq[0] != "G", f"5'={seq[0]} (5' G quenches FAM)")
        warn("no G in first 3", "G" not in seq[:3], f"first3={seq[:3]}")
        chk("more C than G", seq.count("C") >= seq.count("G"),
            f"C/G = {seq.count('C')}/{seq.count('G')}")
        _ptm = T.tm(seq)
        _lna_mgb = any(x in profile.get("name", "") for x in ("LNA", "Affinity", "MGB"))
        if _lna_mgb:
            out.append(("probe Tm (as DNA)", "PASS",
                        f"{_ptm:.1f} C — {profile.get('name','')} chemistry raises the effective Tm; confirm in OligoAnalyzer"))
        else:
            _floor = profile.get("tm_min", 59.0) + profile.get("probe_offset_min", 6.0)
            warn("probe Tm above primers", _ptm >= _floor,
                 f"{_ptm:.1f} C (want >= {_floor:.0f} C, ~{profile.get('probe_offset_min',6):.0f}+ over the {profile.get('tm_min',59):.0f} C primer floor)")
        hdg, htm = T.hairpin(seq)
        warn("hairpin", hdg > profile["probe_hairpin_min"], f"dG={hdg:.2f} (Tm {htm:.0f})")
        warn("self-dimer", T.self_dimer(seq) > profile["self_dimer_min"], f"dG={T.self_dimer(seq):.2f}")
    else:
        chk("length", profile["len_min"] <= len(seq) <= profile["len_max"],
            f"{len(seq)} nt (allowed {profile['len_min']}-{profile['len_max']})")
        chk("GC%", profile["gc_min"] <= T.gc_percent(seq) <= profile["gc_max"],
            f"{T.gc_percent(seq):.0f}% (allowed {profile['gc_min']:.0f}-{profile['gc_max']:.0f})")
        chk("Tm", profile["tm_min"] <= T.tm(seq) <= profile["tm_max"],
            f"{T.tm(seq):.1f} C (target {profile['tm_min']}-{profile['tm_max']})")
        if profile.get("no_three_prime_T"):
            warn("3' not T", seq[-1] != "T", f"3'={seq[-1]}")
        chk("3' clamp <=3 G/C in last 5", T.last5_gc(seq) <= profile["max_3prime_gc"],
            f"{T.last5_gc(seq)} G/C in last 5")
        if profile.get("min_3prime_gc"):
            warn("3' GC clamp present", T.last5_gc(seq) >= profile["min_3prime_gc"],
                 f"{T.last5_gc(seq)} G/C in last 5 (>= {profile['min_3prime_gc']} helps polymerase extension)")
        chk("no G-run >=%d" % profile["max_g_run"], T.max_run(seq, "G") < profile["max_g_run"],
            f"max G-run {T.max_run(seq,'G')}")
        chk("no homopolymer >=%d" % profile["max_any_run"], T.max_run(seq) < profile["max_any_run"],
            f"max run {T.max_run(seq)}")
        warn("hairpin", T.hairpin(seq)[0] > profile["hairpin_min"], f"dG={T.hairpin(seq)[0]:.2f}")
        warn("self-dimer", T.self_dimer(seq) > profile["self_dimer_min"], f"dG={T.self_dimer(seq):.2f}")
    return out


def lint_pair(fseq, rseq, amplicon_len, profile):
    from . import thermo as T
    out = []
    gap = abs(T.tm(fseq) - T.tm(rseq))
    out.append(("pair Tm match", "PASS" if gap <= profile["pair_tm_gap_max"] else "FAIL",
                f"gap {gap:.1f} C (<= {profile['pair_tm_gap_max']})"))
    fxr = T.hetero_dimer(fseq, rseq)
    out.append(("F x R dimer", "PASS" if fxr > profile["pair_dimer_min"] else "FAIL",
                f"dG={fxr:.2f}" + ("  [SYBR-critical]" if profile.get("no_probe") else "")))
    out.append(("amplicon size", "PASS" if profile["amp_min"] <= amplicon_len <= profile["amp_max"] else "FAIL",
                f"{amplicon_len} bp (allowed {profile['amp_min']}-{profile['amp_max']})"))
    return out
