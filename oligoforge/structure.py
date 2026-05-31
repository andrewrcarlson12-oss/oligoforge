"""Template secondary-structure check via ViennaRNA (RNAfold).

The gap Beacon Designer fills: it folds the template and avoids regions with stable
secondary structure, because a stem over a primer's 3' end or the probe footprint
impedes priming / hybridisation. This computes the amplicon's MFE fold and the fraction
of each oligo's binding site that is base-paired in that fold.

ViennaRNA folds RNA; DNA nearest-neighbour parameters (Mathews 2004) are loaded when the
installed build exposes them, otherwise the RNA model is used as a structure-tendency
proxy. `dna_params` in the result says which was used, so the UI can label it honestly.
"""
try:
    import RNA
    _HAVE = True
except Exception:
    _HAVE = False

_DNA = False
if _HAVE:
    for _fn in ("params_load_DNA_Mathews2004", "params_load_DNA_Mathews1999"):
        try:
            getattr(RNA, _fn)()
            _DNA = True
            break
        except Exception:
            pass


def available():
    return _HAVE


def fold(seq):
    if not _HAVE or not seq:
        return None
    s = "".join(c for c in seq.upper() if c in "ACGTUN").replace("T", "U") if not _DNA else \
        "".join(c for c in seq.upper() if c in "ACGTUN")
    if len(s) < 8:
        return None
    try:
        fc = RNA.fold_compound(s)
        struct, mfe = fc.mfe()
    except Exception:
        return None
    paired = {i for i, ch in enumerate(struct) if ch in "()"}
    return dict(mfe=round(float(mfe), 1), mfe_per_nt=round(float(mfe) / max(1, len(s)), 3),
                structure=struct, paired=paired, dna_params=_DNA, n=len(s))


def site_paired_fraction(paired, start, end):
    if paired is None or end <= start:
        return None
    return round(sum(1 for i in range(start, end) if i in paired) / (end - start), 2)
