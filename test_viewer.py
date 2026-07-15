"""Regression pins for the SnapGene-viewer backend (oligoforge.design.design_candidates / probe_span).

Verifies the CRITICAL coordinate invariants used to highlight oligos on the base-level view --
template[f_xy]==forward, template[r_xy]==revcomp(reverse), and template[probe_xy]==probe (or its
revcomp for a minus-strand probe) -- plus Tm-window and GC-window enforcement (the multiplex rules).
Offline + deterministic: primer3 Tm is deterministic and the template is a real Florida Scrub-Jay
RPL13 transcript (NCBI XM_069027304.1) embedded as a constant.
"""
import sys; sys.path.insert(0, ".")
import oligoforge.design as D, oligoforge.thermo as T, oligoforge.profiles as P
rc = T.revcomp

SEQ = "CCTCGCGCGCGGCCCGTTGGGCGCCGCCGCGGCACTTCCGGCGCCGCCGCGCCAGCGCTTCCTTTCCGGTGGCCATTGTCGTGGCGGGAGGCCGCAGCCATGGCGCCCAGCCGCAATGGGATGATCCTGAAGCCCCACTTCCACAAGGACTGGCAGCGCCGGGTCGCCACCTGGTTCAACCAACCCGCCCGCAAGCTCCGCAGGAGGAAGGCTCGCCAGGCCAAGGCTCGCCGCATCGCCCCTCGCCCCGTGGCTGGGCCCATCCGGCCCATCGTGAGGTGCCCTACGATCAGATACCACAAAAAAGTTCGTGCTGGTAGAGGCTTCAGCCTGGAGGAGCTTAAACTCGCTGGCATCAACAAGAAGTTTGCCCGGACTATCGGAATCTCCGTGGATCCCCGGAGACGGAACAAGTCCACCGAGTCCCTGCAGGCCAACGTGCAGAGGCTGAAGGAGTACCACTCCAAGCTCATCCTCTTCCCGAGGAAGCCAGCCATGCCCAAGAAGGGAGACAGCTCTCCAGAGGAACTCAAGATGGCCACTCAGCTCACAGGACCCGTTATGCCCATCAAGAACGTTTTCAAGCGGGAGAAGGCGCGTGTCATCTCGGAAGACGAGAAGAACTTCAAGGCCTTTGCCAGCCTTCGCATGGCCCGGGCCAATGCCCGCCTCTTTGGCATCCGCGCCAAGCGCGCCAAGGAAGCAGCGGAGCAGGACGTGGAGAAAAAGAAATGAACTGTTCTCCCCAGAACTGTCAATAAAAAGCCGTAGAGA"

def cfg(tmlo, tmhi, probe=True):
    c = dict(P.PROFILES["idt_taqman"])
    c.update(tm_min=tmlo, tm_max=tmhi, tm_opt=(tmlo + tmhi) / 2.0, gc_min=30.0, gc_max=70.0, amp_min=70, amp_max=200)
    if probe:
        c.update(no_probe=False, probe_offset_min=5.0, probe_offset_max=12.0)
    else:
        c["no_probe"] = True
    return c

def invariants(cands, tmlo, tmhi, probe):
    assert cands, "no candidates found"
    for a in cands:
        f0, f1 = a["f_xy"]; r0, r1 = a["r_xy"]
        assert SEQ[f0:f1] == a["forward"], ("forward coordinate mismatch", a["f_xy"])
        assert SEQ[r0:r1] == rc(a["reverse"]), ("reverse coordinate mismatch", a["r_xy"])
        assert a["amplicon"] == r1 - f0, ("amplicon size", a["amplicon"])
        assert a["amplicon_xy"] == [f0, r1], "amplicon_xy"
        assert tmlo <= a["f_tm"] <= tmhi and tmlo <= a["r_tm"] <= tmhi, ("Tm window", a["f_tm"], a["r_tm"])
        assert 30.0 <= T.gc_percent(a["forward"]) <= 70.0 and 30.0 <= T.gc_percent(a["reverse"]) <= 70.0, "GC window"
        if probe and a.get("probe"):
            px = a.get("probe_xy"); assert px, "probe_xy missing"
            sub = SEQ[px[0]:px[1]]
            assert sub == a["probe"] or sub == rc(a["probe"]), ("probe coordinate mismatch", a.get("probe_strand"), px)

cT = D.design_candidates(SEQ, cfg(58.0, 64.5, True), n=5)
invariants(cT, 58.0, 64.5, True)
assert any(a.get("probe") and a.get("probe_xy") for a in cT), "expected a probe candidate with coordinates"
cF = D.design_candidates(SEQ, cfg(58.0, 64.5, False), n=5)
invariants(cF, 58.0, 64.5, False)
# tight multiplex window: every returned primer must fall inside 61-62 C
cTight = D.design_candidates(SEQ, cfg(61.0, 62.0, True), n=5)
invariants(cTight, 61.0, 62.0, True)
# amp_max above the default 400 bp window must be honored (the window adapts) -- a >400 bp amplicon is
# impossible under the old fixed window; coordinates must still be exact
cBig = dict(P.PROFILES["idt_taqman"]); cBig.update(amp_min=350, amp_max=600, no_probe=True,
            gc_min=30.0, gc_max=70.0, tm_min=58.0, tm_max=65.0, tm_opt=61.5)
cb = D.design_candidates(SEQ, cBig, n=5)
assert cb and max(a["amplicon"] for a in cb) > 400, "amp_max above 400 must be reachable (window must adapt)"
invariants(cb, 58.0, 65.0, False)

# degenerate (IUPAC) bases must not crash design, and coordinates stay exact against the input
iup = SEQ[:300] + "RYSWKM" + SEQ[300:]
ciup = dict(P.PROFILES["idt_taqman"]); ciup.update(gc_min=30.0, gc_max=70.0, tm_min=57.0, tm_max=65.0, amp_min=70, amp_max=160)
for a in D.design_candidates(iup, ciup, n=3):
    assert iup[a["f_xy"][0]:a["f_xy"][1]] == a["forward"], "IUPAC F coord"
    assert iup[a["r_xy"][0]:a["r_xy"][1]] == rc(a["reverse"]), "IUPAC R coord"

# probe_span returns None cleanly for a probe-less assay dict
assert D.probe_span(SEQ, {"forward": "A", "reverse": "A", "f_xy": [0, 1], "r_xy": [2, 3]}) is None
print("ALL VIEWER ASSERTS PASS")
