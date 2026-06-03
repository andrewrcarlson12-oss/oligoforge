"""Relative gene-expression from qPCR Cq values: the back half of an expression study,
so a finished run can be analysed in OligoForge without leaving for a stats package.

Two standard methods, both normalised to one or more reference genes:

  Livak   -- the 2^-ddCq method (Livak & Schmittgen 2001, Methods 25:402-408). Assumes every
             assay amplifies at ~100% efficiency (E=2). dCq = Cq(target) - reference; ddCq is
             that minus the control group's mean dCq; fold change = 2^-ddCq.
  Pfaffl  -- the efficiency-corrected ratio (Pfaffl 2001, NAR 29:e45), using the per-assay
             amplification efficiencies the user supplies. With every efficiency at 100% it
             reduces exactly to Livak.

Multiple reference genes are combined as the geometric mean of their relative quantities
(Vandesompele et al. 2002, Genome Biol 3:research0034) -- which, on the E=2 Livak scale, is
the arithmetic mean of the reference Cqs. Pure Python (no numpy/scipy) to match the engine.
Nothing here is a substitute for proper experimental design, technical replicates, or a
validated reference-gene panel; it computes what the Cq values imply, with the assumptions stated.
"""
import math
import statistics


def _mean(xs):
    return sum(xs) / len(xs)


def _geomean(xs):
    # all xs must be > 0 (they are: amplification factors raised to a power, or 2^-Cq)
    return math.exp(sum(math.log(x) for x in xs) / len(xs))


def _E(eff_pct):
    """Amplification factor from percent efficiency. 100% -> E=2.0 (a perfect doubling)."""
    if eff_pct is None:
        return 2.0
    return 1.0 + float(eff_pct) / 100.0


def _cq(sample, gene):
    """A finite float Cq for `gene` in `sample`, or None (missing / undetermined)."""
    v = (sample.get("cq") or {}).get(gene)
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def parse_table(text):
    """Parse a pasted CSV/TSV table into the sample list `analyze` expects.

    Header must be:  sample, group, <gene1>, <gene2>, ...
    (first column an id, second the group/condition, the rest gene Cq columns). Blank cells and
    NA / NaN / undetermined are treated as missing for that gene.
    """
    text = (text or "").strip()
    if not text:
        raise ValueError("no data pasted")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    delim = "\t" if "\t" in lines[0] else ","
    header = [h.strip() for h in lines[0].split(delim)]
    hlow = [h.lower() for h in header]
    if len(header) < 3 or hlow[0] not in ("sample", "id", "name") or \
       hlow[1] not in ("group", "condition", "treatment"):
        raise ValueError("header must be: sample, group, <gene1>, <gene2>, ...")
    genes = header[2:]
    _NA = ("", "na", "nan", "-", "undetermined", "und", "none")
    samples = []
    for ln in lines[1:]:
        cells = [c.strip() for c in ln.split(delim)]
        if len(cells) < 2 or not cells[0]:
            continue
        cq = {}
        for j, g in enumerate(genes, start=2):
            if j < len(cells) and cells[j].lower() not in _NA:
                try:
                    cq[g] = float(cells[j])
                except ValueError:
                    pass
        samples.append(dict(sample=cells[0], group=cells[1], cq=cq))
    if not samples:
        raise ValueError("no data rows")
    return samples


def analyze(samples, reference_genes, control_group, efficiencies=None):
    """Relative expression of every target gene per group, vs the control group.

    samples: [{sample, group, cq:{gene: Cq, ...}}, ...]
    reference_genes: list of gene names used for normalisation (>=1).
    control_group: the calibrator group every fold change is expressed relative to.
    efficiencies: optional {gene: percent} (e.g. 98.5). Absent -> 100% (Livak only).
    """
    if not samples:
        raise ValueError("no samples")
    refs = [r for r in (reference_genes or []) if r]
    if not refs:
        raise ValueError("at least one reference gene is required")
    if not control_group:
        raise ValueError("a control / calibrator group is required")
    eff = {}
    for k, v in (efficiencies or {}).items():
        try:
            eff[k] = float(v)
        except (TypeError, ValueError):
            pass

    genes = []
    for s in samples:
        for g in (s.get("cq") or {}):
            if g not in genes:
                genes.append(g)
    targets = [g for g in genes if g not in refs]
    if not targets:
        raise ValueError("no target genes found (every gene in the data is a reference gene)")
    missing = [r for r in refs if r not in genes]
    if missing:
        raise ValueError("reference gene(s) not in the data: " + ", ".join(missing))

    groups = []
    for s in samples:
        g = s.get("group")
        if g and g not in groups:
            groups.append(g)
    if control_group not in groups:
        raise ValueError("control group '%s' is not among the sample groups" % control_group)

    warnings = []
    dcq_by = {}        # (target, group) -> [dCq]      (Livak)
    norm_by = {}       # (target, group) -> [Q_norm]   (Pfaffl, efficiency-corrected)
    per_sample = []

    for s in samples:
        grp = s.get("group")
        ref_cqs = [_cq(s, r) for r in refs]
        if any(v is None for v in ref_cqs):
            warnings.append("sample %s is missing a reference Cq and was skipped" % s.get("sample", "?"))
            continue
        ref_mean = _mean(ref_cqs)                                            # Livak normaliser (E=2)
        nf = _geomean([_E(eff.get(r)) ** (-c) for r, c in zip(refs, ref_cqs)])  # Pfaffl norm. factor
        for t in targets:
            ct = _cq(s, t)
            if ct is None:
                continue
            dcq = ct - ref_mean
            qn = (_E(eff.get(t)) ** (-ct)) / nf
            dcq_by.setdefault((t, grp), []).append(dcq)
            norm_by.setdefault((t, grp), []).append(qn)
            per_sample.append(dict(sample=s.get("sample"), group=grp, target=t,
                                   dcq=round(dcq, 4), rel_norm=qn))

    results = []
    for t in targets:
        ctrl_dcq = dcq_by.get((t, control_group))
        if not ctrl_dcq:
            warnings.append("no control-group measurements for %s; its fold changes are undefined" % t)
            continue
        mean_ctrl = _mean(ctrl_dcq)
        ctrl_norm = norm_by.get((t, control_group)) or []
        ctrl_center = _geomean(ctrl_norm) if ctrl_norm else None
        for grp in groups:
            ds = dcq_by.get((t, grp))
            if not ds:
                continue
            mean_dcq = _mean(ds)
            sd_dcq = statistics.stdev(ds) if len(ds) > 1 else 0.0
            ddcq = mean_dcq - mean_ctrl
            fold = 2.0 ** (-ddcq)
            row = dict(target=t, group=grp, n=len(ds),
                       mean_dcq=round(mean_dcq, 3), sd_dcq=round(sd_dcq, 3),
                       ddcq=round(ddcq, 3),
                       fold_livak=round(fold, 4),
                       fold_low=round(2.0 ** (-(ddcq + sd_dcq)), 4),
                       fold_high=round(2.0 ** (-(ddcq - sd_dcq)), 4),
                       log2_fold=round(math.log2(fold), 4) if fold > 0 else None)
            ns = norm_by.get((t, grp))
            if ns and ctrl_center:
                row["fold_pfaffl"] = round(_geomean(ns) / ctrl_center, 4)
            if len(ds) < 3:
                row["low_n"] = True
            results.append(row)

    used_eff = bool(eff)
    note = ("Livak 2^-ddCq" + ("" if used_eff else " (all efficiencies assumed 100% / E=2)") + ". "
            + ("Pfaffl efficiency-corrected ratios use the supplied per-assay efficiencies; "
               "with every efficiency at 100% Pfaffl equals Livak. " if used_eff
               else "Supply per-assay efficiencies (from your standard curves) to also get "
                    "efficiency-corrected Pfaffl ratios. ")
            + ("Reference normalisation uses the geometric mean of %d reference genes (Vandesompele). "
               % len(refs) if len(refs) > 1 else "Normalised to the single reference gene %s. " % refs[0])
            + "Fold range is 2^-(ddCq +/- SD of dCq). All fold changes are relative to group '%s'."
              % control_group)
    return dict(targets=targets, reference_genes=refs, control_group=control_group,
                groups=groups, n_samples=len(samples), efficiencies=(eff or None),
                results=results, per_sample=per_sample, method_note=note, warnings=warnings)
