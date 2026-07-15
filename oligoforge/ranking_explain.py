"""Human-readable explanations and explicit uncertainty for structured ranks.

The authoritative order remains in :mod:`oligoforge.ranking`.  This module does
not invent a calibrated probability that rank 1 will win at the bench.  It reports
whether the modeled difference is decisive, modest, near-equivalent, or impossible
to interpret because important evidence was not supplied.
"""

UNCERTAINTY_SCHEMA_VERSION = "1.1.0"


def evidence_completeness(evidence):
    e = evidence or {}
    ev = e.get("evaluations") or {}
    objective = str(e.get("objective") or "balanced")
    required = {
        "target_epcr": True,
        "condition_robustness": True,
        # Any claim about exclusivity is corpus-dependent.  Absence of an
        # off-target set is therefore visible even for inclusive/screening runs.
        "offtarget_epcr": objective not in {"broad_inclusivity", "screening"},
        "panel": objective == "multiplex",
        "junction": objective == "transcript_specific",
    }
    missing = [name for name, needed in required.items() if needed and not ev.get(name)]
    optional_missing = [name for name in ("offtarget_epcr", "panel", "junction")
                        if not required.get(name) and not ev.get(name)]
    state = "complete for declared computational objective" if not missing else "incomplete computational evidence"
    return {
        "schema_version": UNCERTAINTY_SCHEMA_VERSION,
        "state": state,
        "required_evaluations": required,
        "missing_required": missing,
        "missing_optional": optional_missing,
        "empirical_validation_present": bool(ev.get("empirical_feedback")),
        "scope": "computational evidence only; no probability of wet-lab superiority",
    }


def _offtarget_tuple(e):
    o = (e or {}).get("offtarget") or {}
    return int(o.get("signal_subjects", 0)), int(o.get("product_subjects", 0))


def preference_strength(item, competitor=None):
    if not competitor:
        return {
            "schema_version": UNCERTAINTY_SCHEMA_VERSION,
            "state": "insufficient evidence to compare",
            "basis": "no adjacent fully evaluated candidate",
            "confidence": "not assessed",
            "decisive_margin": None,
            "missing_evidence": evidence_completeness((item or {}).get("evidence") or {}).get("missing_required", []),
        }
    e, c = item.get("evidence") or {}, competitor.get("evidence") or {}
    ce, cc = evidence_completeness(e), evidence_completeness(c)
    missing = sorted(set(ce["missing_required"] + cc["missing_required"]))

    if bool(e.get("hard_valid")) != bool(c.get("hard_valid")):
        return dict(schema_version=UNCERTAINTY_SCHEMA_VERSION, state="strong preference",
                    basis="hard validity differs", confidence="high within the declared computational rules",
                    decisive_margin=1, missing_evidence=missing)

    eo, co = _offtarget_tuple(e), _offtarget_tuple(c)
    if eo != co:
        if not ((e.get("evaluations") or {}).get("offtarget_epcr") and
                (c.get("evaluations") or {}).get("offtarget_epcr")):
            return dict(schema_version=UNCERTAINTY_SCHEMA_VERSION,
                        state="insufficient evidence to distinguish",
                        basis="off-target burden differs but the compared off-target corpus was not fully evaluated",
                        confidence="low", decisive_margin=None,
                        missing_evidence=sorted(set(missing + ["offtarget_epcr"])))
        return dict(schema_version=UNCERTAINTY_SCHEMA_VERSION, state="strong preference",
                    basis="predicted off-target product burden differs",
                    confidence="high for the supplied off-target corpus",
                    decisive_margin=[co[0] - eo[0], co[1] - eo[1]], missing_evidence=missing)

    dcov = float(e.get("target_coverage", 0)) - float(c.get("target_coverage", 0))
    if abs(dcov) >= 0.05:
        return dict(schema_version=UNCERTAINTY_SCHEMA_VERSION, state="strong preference",
                    basis="coherent target coverage differs by at least 5 percentage points",
                    confidence="high for the supplied target corpus", decisive_margin=round(dcov, 4),
                    missing_evidence=missing)
    if abs(dcov) >= 0.01:
        return dict(schema_version=UNCERTAINTY_SCHEMA_VERSION, state="moderate preference",
                    basis="coherent target coverage differs",
                    confidence="moderate for the supplied target corpus", decisive_margin=round(dcov, 4),
                    missing_evidence=missing)

    if int(e.get("panel_risk", 0)) != int(c.get("panel_risk", 0)):
        if not ((e.get("evaluations") or {}).get("panel") and (c.get("evaluations") or {}).get("panel")):
            return dict(schema_version=UNCERTAINTY_SCHEMA_VERSION,
                        state="insufficient evidence to distinguish",
                        basis="panel-risk values are not backed by the same evaluated panel",
                        confidence="low", decisive_margin=None,
                        missing_evidence=sorted(set(missing + ["panel"])))
        return dict(schema_version=UNCERTAINTY_SCHEMA_VERSION, state="moderate preference",
                    basis="multiplex-panel interaction burden differs",
                    confidence="moderate for the supplied panel", decisive_margin=int(c.get("panel_risk", 0))-int(e.get("panel_risk", 0)),
                    missing_evidence=missing)

    er = float((e.get("condition_robustness") or {}).get("valid_fraction", 0))
    cr = float((c.get("condition_robustness") or {}).get("valid_fraction", 0))
    dr = er - cr
    if abs(dr) >= 0.34:
        return dict(schema_version=UNCERTAINTY_SCHEMA_VERSION, state="moderate preference",
                    basis="condition-envelope stability differs",
                    confidence="moderate within the three declared scenarios", decisive_margin=round(dr, 4),
                    missing_evidence=missing)

    dp = float(c.get("triplet_penalty", 0)) - float(e.get("triplet_penalty", 0))
    practical = float(c.get("practical_penalty", 0)) - float(e.get("practical_penalty", 0))
    if abs(dp) <= 0.75 and abs(practical) <= 0.75:
        state = "insufficient evidence to distinguish" if missing else "near-equivalent alternatives"
        basis = ("important objective evidence is missing and the remaining modeled differences are small"
                 if missing else "no decisive modeled difference under the selected objective")
        return dict(schema_version=UNCERTAINTY_SCHEMA_VERSION, state=state, basis=basis,
                    confidence="low", decisive_margin={"triplet_penalty": round(dp, 4), "practical_penalty": round(practical, 4)},
                    missing_evidence=missing)

    if missing:
        return dict(schema_version=UNCERTAINTY_SCHEMA_VERSION,
                    state="insufficient evidence to distinguish",
                    basis="ordering depends on soft thermodynamic or practical evidence while required evidence is missing",
                    confidence="low", decisive_margin={"triplet_penalty": round(dp, 4)},
                    missing_evidence=missing)
    return dict(schema_version=UNCERTAINTY_SCHEMA_VERSION, state="moderate preference",
                basis="lower complete-triplet or practical penalty",
                confidence="moderate computational support", decisive_margin={"triplet_penalty": round(dp, 4), "practical_penalty": round(practical, 4)},
                missing_evidence=missing)



def rank_reversal_scenarios(item, competitor=None):
    """Return concrete, evidence-linked conditions that could reverse an ordering.

    This is sensitivity analysis, not a probability model.  It only names changes
    that are supported by a measurable trade-off or by an evaluation that was not
    performed.  Generic boilerplate is deliberately avoided.
    """
    if not competitor:
        return [{"trigger": "another fully evaluated candidate becomes available",
                 "reason": "no adjacent comparator was supplied", "evidence_level": "informational"}]
    e, c = (item or {}).get("evidence") or {}, (competitor or {}).get("evidence") or {}
    ev, cv = e.get("evaluations") or {}, c.get("evaluations") or {}
    out = []

    eo, co = _offtarget_tuple(e), _offtarget_tuple(c)
    ecov, ccov = float(e.get("target_coverage", 0)), float(c.get("target_coverage", 0))
    er = float((e.get("condition_robustness") or {}).get("valid_fraction", 0))
    cr = float((c.get("condition_robustness") or {}).get("valid_fraction", 0))
    ep = float(e.get("practical_penalty", 0)); cp = float(c.get("practical_penalty", 0))
    et = float(e.get("triplet_penalty", 0)); ct = float(c.get("triplet_penalty", 0))

    if ccov > ecov + 1e-9:
        out.append({"trigger": "the objective shifts toward broad inclusivity or screening sensitivity",
                    "reason": "the competing assay has higher coherent target coverage",
                    "current": round(ecov, 4), "competitor": round(ccov, 4),
                    "evidence_level": "computed for the supplied target corpus"})
    if co < eo:
        out.append({"trigger": "the objective shifts toward confirmatory exclusivity",
                    "reason": "the competing assay has fewer predicted off-target products/signals",
                    "current": list(eo), "competitor": list(co),
                    "evidence_level": "computed only for the supplied off-target corpus"})
    if cr > er + 1e-9:
        out.append({"trigger": "reaction conditions vary outside the nominal setting",
                    "reason": "the competing assay remains valid in more tested condition scenarios",
                    "current": round(er, 4), "competitor": round(cr, 4),
                    "evidence_level": "three-scenario computational envelope"})
    if cp + 0.5 < ep:
        out.append({"trigger": "short-product recovery, synthesis burden, or practical simplicity is prioritized",
                    "reason": "the competing assay has a materially lower practical penalty",
                    "current": round(ep, 4), "competitor": round(cp, 4),
                    "evidence_level": "computational preference"})
    if ct + 0.75 < et:
        out.append({"trigger": "thermodynamic triplet quality is given greater priority",
                    "reason": "the competing assay has a materially lower complete-triplet penalty",
                    "current": round(et, 4), "competitor": round(ct, 4),
                    "evidence_level": "model-based"})

    for key, label in (("offtarget_epcr", "a representative off-target corpus is supplied"),
                       ("panel", "the intended multiplex panel is supplied"),
                       ("junction", "the transcript/junction requirement is verified")):
        if not ev.get(key) or not cv.get(key):
            out.append({"trigger": label,
                        "reason": "%s evidence is incomplete for one or both candidates" % key,
                        "evidence_level": "missing evidence could change the order"})

    if (item.get("equivalence_group") and competitor.get("equivalence_group") and
            item.get("equivalence_group", {}).get("group_id") == competitor.get("equivalence_group", {}).get("group_id")):
        out.append({"trigger": "small changes in sequence corpus, reaction conditions, or empirical results",
                    "reason": "the candidates are already in the same modeled equivalence band",
                    "evidence_level": "near-equivalent computational evidence"})
    if not out:
        out.append({"trigger": "new empirical efficiency, specificity, or inhibition evidence",
                    "reason": "the current order is based on modeled evidence that does not cover all wet-lab effects",
                    "evidence_level": "unmeasured biological performance"})
    return out

def explain(item, competitor=None):
    e = item.get("evidence") or {}
    a = item.get("assay") or {}
    strengths, weaknesses = [], []
    cov = e.get("target_coverage", 0.0)
    if cov >= 0.999:
        strengths.append("coherent amplification%s across every supplied target" % (" with probe binding" if a.get("probe") else ""))
    elif cov:
        weaknesses.append("coherent target coverage is %.1f%%" % (100 * cov))
    off = e.get("offtarget") or {}
    if (e.get("evaluations") or {}).get("offtarget_epcr"):
        if off.get("signal_subjects", 0) == 0:
            strengths.append("no signal-generating product in the supplied off-target set")
        else:
            weaknesses.append("%d supplied off-target sequence(s) can yield a signal-generating product" % off.get("signal_subjects", 0))
    else:
        weaknesses.append("no supplied off-target corpus was evaluated; exclusivity is unresolved")
    rf = ((e.get("condition_robustness") or {}).get("valid_fraction"))
    if rf == 1:
        strengths.append("retains acceptable Tm, hairpin, self-dimer, cross-dimer, and primer–probe interaction states across the declared condition envelope")
    elif rf is not None:
        weaknesses.append("passes the complete thermodynamic and interaction screen in %d%% of tested condition scenarios" % round(100 * rf))
    if e.get("hard_failures"):
        weaknesses.extend(e["hard_failures"])
    strongest = strengths[0] if strengths else "no single decisive advantage was identified"
    weakest = weaknesses[0] if weaknesses else "no modeled hard failure; wet-lab performance remains uncertain"
    why = "Ranked by hard validity first, then the %s objective evidence hierarchy." % e.get("objective", "selected")
    above = None
    if competitor:
        pref_preview = preference_strength(item, competitor)
        above = "%s: %s." % (pref_preview["state"].capitalize(), pref_preview["basis"])
    limitations = []
    evals = e.get("evaluations") or {}
    if not evals.get("offtarget_epcr"):
        limitations.append("no supplied off-target set was evaluated")
    if not evals.get("panel"):
        limitations.append("multiplex-panel interactions were not evaluated")
    limitations.append("software does not predict all wet-lab efficiency, inhibition, fluorescence, or synthesis effects")
    pref = preference_strength(item, competitor)
    completeness = evidence_completeness(e)
    comp_summary = None
    if competitor:
        ca = competitor.get("assay") or {}
        comp_summary = dict(rank=competitor.get("rank"), forward=ca.get("forward"),
                            reverse=ca.get("reverse"), probe=ca.get("probe"),
                            display_score=competitor.get("display_score"))
    return dict(
        schema_version="oligoforge-ranking-explanation/v2",
        why_ranked=why,
        strongest_feature=strongest,
        weakest_feature=weakest,
        strengths=strengths,
        weaknesses=weaknesses,
        hard_constraints_satisfied=not bool(e.get("hard_failures")),
        hard_failures=e.get("hard_failures") or [],
        comparison_to_next=above,
        closest_competitor=comp_summary,
        preference_state=pref["state"],
        preference_basis=pref["basis"],
        preference_confidence=pref.get("confidence"),
        decisive_margin=pref.get("decisive_margin"),
        evidence_completeness=completeness,
        equivalence_group=item.get("equivalence_group"),
        rank_reversal_scenarios=rank_reversal_scenarios(item, competitor),
        ranking_may_reverse_if="; ".join(x["trigger"] for x in rank_reversal_scenarios(item, competitor)),
        evaluations_not_performed=[k for k, v in evals.items() if not v],
        limitations=limitations,
        optimum_status="heuristic bounded search; strongest candidate among the retained and fully evaluated pool",
    )
