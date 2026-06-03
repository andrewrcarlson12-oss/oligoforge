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

import threading
_LOCK = threading.Lock()      # ViennaRNA uses non-thread-safe global state; serialize folds

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
    if len(s) < 8 or len(s) > 1000:   # RNAfold is O(n^3); amplicons are short — never fold a long template
        return None
    try:
        with _LOCK:
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


def fold_ensemble(seq, anneal_c=None):
    """Richer fold for post-ranking annotation. Returns fold()'s keys PLUS:
      - paired_prob: per-position ensemble base-pairing probability from the partition
        function (fc.pf/bpp). A base can be paired in the single MFE structure yet largely
        unpaired across the Boltzmann ensemble; the ensemble probability is the more honest
        "is this site accessible" signal than the MFE alone.
      - paired_anneal / mfe_anneal: the MFE structure at the annealing temperature (most
        37 C structure has melted by ~60 C, so the 37 C fold overstates structure at priming).
    Heavier than fold() (adds pf + a second fold), so used only on the handful of ranked
    candidates, never in the enumeration loop. MFE keys are identical to fold()."""
    if not _HAVE or not seq:
        return None
    s = "".join(c for c in seq.upper() if c in "ACGTUN").replace("T", "U") if not _DNA else \
        "".join(c for c in seq.upper() if c in "ACGTUN")
    if len(s) < 8 or len(s) > 1000:
        return None
    n = len(s)
    try:
        with _LOCK:
            fc = RNA.fold_compound(s)
            struct, mfe = fc.mfe()
            fc.pf()
            bpp = fc.bpp()                          # 1-indexed upper triangle
            probs = [0.0] * n
            for i in range(1, n + 1):
                row = bpp[i] if i < len(bpp) else None
                if not row:
                    continue
                for j in range(i + 1, n + 1):
                    try:
                        pij = row[j]
                    except Exception:
                        pij = 0.0
                    if pij > 0:
                        probs[i - 1] += pij
                        probs[j - 1] += pij
            paired_anneal = mfe_anneal = None
            if anneal_c is not None:
                _old = RNA.cvar.temperature
                try:
                    RNA.cvar.temperature = float(anneal_c)
                    sa, ma = RNA.fold_compound(s).mfe()
                    paired_anneal = {k for k, ch in enumerate(sa) if ch in "()"}
                    mfe_anneal = round(float(ma), 1)
                finally:
                    RNA.cvar.temperature = _old      # always restore the 37 C default
    except Exception:
        return None
    paired = {k for k, ch in enumerate(struct) if ch in "()"}
    return dict(mfe=round(float(mfe), 1), mfe_per_nt=round(float(mfe) / max(1, n), 3),
                structure=struct, paired=paired,
                paired_prob=[round(min(p, 1.0), 3) for p in probs],
                paired_anneal=paired_anneal, mfe_anneal=mfe_anneal,
                anneal_c=(float(anneal_c) if anneal_c is not None else None),
                dna_params=_DNA, n=n)


def site_paired_prob(probs, start, end):
    """Mean ensemble paired probability over a binding site [start, end)."""
    if not probs or end <= start or start < 0 or end > len(probs):
        return None
    return round(sum(probs[start:end]) / (end - start), 2)
