# OligoForge Ranking Audit — v1.35.0

Version 1.35.0 does not change the authoritative ranking order, ranking schema, objective priorities, or ranker version (2.2.0). Search version 2.1.1 evaluates a deterministic minimum of up to three spread-ordered target windows before its soft runtime cutoff, preventing cache warmth from changing the minimum candidate corpus; search remains heuristic-bounded. Staged automatic design invokes the same refactored scientific functions as the retained synchronous compatibility endpoint; Validation Studio and Assurance consume completed assay evidence without silently changing rank.

## Executive finding

The legacy pipeline often produced usable assays, but it did not establish that rank 1 was the strongest supported complete assay. Target-wide windows existed, yet early primer, pair, probe, and triplet beams could discard the eventual best assay before target coverage, paired specificity, robustness, junction, or panel evidence was available. A preliminary scalar then determined which survivors received full annotation.

Version 1.34.0 retains the v1.32 complete-triplet repair and adds evidence-completeness, explicit indifference, and reproducible provenance. It replaces that failure pattern with a bounded, auditable search over complete forward-primer/reverse-primer/probe triplets; diversity-preserving retention; target/off-target-aware probe augmentation; complete annotation of the retained pool; hard constraints; objective-specific lexicographic ranking; Pareto annotation; truthful finalist categories; and a machine-readable attrition ledger. The search remains heuristic-bounded, not an exhaustive proof of a global wet-lab optimum.

## Pipeline traced

1. Input normalization and reaction-condition snapshot.
2. Target/off-target corpus validation and reference selection.
3. Target-spanning window scheduling.
4. Primer enumeration and inexpensive hard screens.
5. Primer-pair geometry, Tm compatibility, and dimer screens.
6. Pair-diversity beam across regions and product-length bands.
7. Multiple probe candidates per retained pair.
8. Probe hard screens and probe-diversity beam.
9. Complete-triplet construction and regional/near-duplicate retention.
10. Discrimination-specialist pair augmentation when off-targets are supplied.
11. Target/off-target-aware probe augmentation.
12. Objective-aware full-annotation retention.
13. Coherent target products, all-site off-target products, conservation, 3′ mismatch, junction, panel, degeneracy, and condition-envelope analysis.
14. Hard constraints, objective lexicographic ordering, Pareto fronts, deterministic tie-breakers, and subordinate display score.
15. Truthful finalist selection, explanations, Workbench transfer, reports, RDML, and order export.

## Legacy ranking architecture

The legacy cross-window order was approximately:

`preliminary conservation/discrimination benefit − sequence penalties − long-product penalty − selected dimer penalties`

This scalar was not calibrated to experimental outcome. More importantly, it acted before several decisive metrics existed. A good nominal Tm or product length could therefore preserve a candidate with an off-target product, weak isolate coverage, or fragile condition behavior while a stronger candidate was discarded upstream.

The v1.34.0 authoritative order is a structured evidence vector:

`hard validity → objective-specific primary evidence → condition/variant robustness → complete-triplet quality → practical preferences → deterministic tie-breakers`

The 0–100 display score cannot override that order.

## Current bounded-search limits

Default automatic multi-sequence design records these values in each run:

- Up to 18 target-spanning windows, 420 nt each, 140 nt nominal step.
- 35-second base search budget.
- Up to 12 primer pairs per window after pair diversity retention.
- Up to 3 local probe candidates per retained pair.
- Up to 30 triplets per window.
- Objective probe augmentation over up to `max(8, requested × 2)` retained pairs, scanning up to 24 probes and adding at most 2 alternatives per pair.
- Final expensive-annotation pool of `max(20, min(28, requested × 5))` complete assays; five requested results therefore use 25.
- Exact duplicate suppression, near-duplicate family limits, regional round-robin retention, and finalist preference for distinct primer pairs and alternate regions.

These are engineering budgets, not biological truths. Every truncation is labeled heuristic and reversible in the ledger.

## Confirmed ranking defects and repairs

| ID | Confirmed defect | Reproducible effect | Repair |
|---|---|---|---|
| RF-01 | Early-window bias | A dense 5′ region consumed the time/candidate budget. | Endpoint/midpoint target-spanning window order and regional retention. |
| RF-02 | Early primer truncation | A mediocre isolated primer could form the strongest complete triplet. | Retain diverse pairs and rank complete assays. |
| RF-03 | Greedy probe attachment | The first acceptable probe hid a cleaner probe on the same pair. | Multiple probes per pair plus objective-aware probe augmentation. |
| RF-04 | One-slot strand-order bug | With one probe slot, `+` strand alphabetical order beat the actual best probe. | Order strand representatives by preliminary evidence. |
| RF-05 | Specialist-pair probe monopoly | Early discrimination pairs consumed all specialist slots with nearby probes. | Preserve unique primer pairs first, then add probe alternatives round-robin. |
| RF-06 | Corpus-superior probe outside local beam | Better isolate conservation or discrimination never reached final rank. | Bounded target/off-target-aware second-stage probe scan. |
| RF-07 | Dense-region monopoly | Near-identical assays consumed the retained budget. | Exact/near-duplicate suppression and region/product-band round-robin retention. |
| RF-08 | Plain top-N annotation cap | Later regions were lost immediately before expensive evidence. | Objective-aware diversity retention before full annotation. |
| RF-09 | Specificity appended too late | Preliminary rank 1 produced a supplied off-target amplicon. | Paired all-site specificity participates in hard/lexicographic rank. |
| RF-10 | Conservation appended too late | A narrow assay outranked coherent multi-isolate amplification. | Coherent target-product and probe coverage precede soft preferences. |
| RF-11 | Weighted compensation | Good Tm could mask a true hard failure. | Non-compensable hard requirements. |
| RF-12 | Degenerate probe disappearance | Primer degeneracy could be evaluated without the unchanged probe. | Effective triplet always preserves every unchanged component. |
| RF-13 | Robustness label exceeded calculation | UI implied broad robustness while only Tm windows were rechecked. | Recalculate Tm, hairpin, self-dimer, cross-dimer, and primer–probe interactions in each scenario. |
| RF-14 | Finalist-category mislabeling | Rank 2 could be called “best specificity” merely because rank 1 was already displayed. | Category labels attach only to the true category winner. |
| RF-15 | Manual/automatic divergence | Viewer and batch paths could use first-hit design while auto used structured rank. | Automatic, direct, batch, viewer, manual, and rescue paths share the authoritative ranker. |
| RF-16 | Ambiguous placement hidden | Repeated motifs highlighted the first sequence match, not the selected coordinates. | Carry exact coordinates and return every manual exact/near placement. |
| RF-17 | Locked component vanished | A lock failed when the component was outside the generic beam. | Inject locked components and optimize only unlocked components. |
| RF-18 | Panel risk underweighted | Small local probe differences outranked severe multiplex interaction. | Multiplex objective places panel evidence after hard specificity. |
| RF-19 | False rank precision | Decimal scores implied biological calibration. | Evidence vectors, rounded subordinate score, preference-strength categories. |
| RF-20 | Shared/stale reaction state | Concurrent requests could alter another run’s thermodynamics. | Request-local immutable condition snapshots. |
| RF-21 | Native-call/runtime opacity | Very large uninterrupted Primer3 annotation batches could stall. | Auditable, diverse 20–28 candidate full-annotation pool and bounded tests. |

Every confirmed defect above has a deterministic regression or is exercised by the versioned biological fixture tests.

## Evidence-completeness and provenance audit

Version 1.32 still had three truthfulness gaps even when its authoritative order was correct:

1. deterministic tie-breakers could make close candidates look more distinguishable than the modeled evidence supported;
2. missing off-target, panel, or junction evaluations were listed as limitations but did not always suppress comparison confidence;
3. run manifests did not record the full application/native-library/model/database provenance needed to reconstruct a later rank.

Version 1.33 repairs those gaps. Consecutive candidates receive a shared equivalence group when no decisive modeled difference exists. Missing required evidence produces `insufficient evidence to distinguish` rather than a soft-confidence preference. Every automatic/manual/rescue run receives a deterministic self-hashing manifest with application, ranker, schema, search/retention/manual/rescue versions, scientific model identifiers, installed scientific-library versions, reaction conditions, candidate limits, constraints, sequence/input hashes, declared external-database metadata, warnings, and fallbacks. Wall-clock time, host names, and local paths are excluded from the run identity.

New confirmed defects and repairs:

| ID | Confirmed defect | Reproducible effect | Repair |
|---|---|---|---|
| RF-22 | Deterministic tie-break portrayed as evidence | Near-identical candidates appeared strictly ordered despite negligible modeled separation. | Explicit rank bands/equivalence groups; deterministic order retained only for reproducibility. |
| RF-23 | Missing corpus evidence did not fully suppress confidence | A soft thermodynamic difference could look authoritative without a supplied off-target set. | Evidence-completeness gate and insufficient-evidence comparison state. |
| RF-24 | Incomplete run provenance | The same visible assay could not be traced to exact model/library/condition/database state. | Deterministic self-hashing provenance manifest propagated through all design workflows. |
| RF-25 | Feedback rows accepted without dataset audit | Duplicate, conflicting, impossible, or target-leaking records could be treated as model-ready. | Validated import, deduplication, conflict/range checks, completeness report, and target-group-isolated split. |
| RF-26 | Benchmark split integrity assumed | A malformed corpus could leak a target group or name a nonexistent expected candidate. | Machine-enforced unique IDs, split completeness/disjointness, target-group isolation, and expectation validation. |
| RF-27 | Selection/export dropped rank provenance | A fully evaluated candidate could enter the Workbench as bare sequences, then appear in HTML/CSV/RDML without its run ID, objective, evidence state, or manifest digest. | Preserve manifest, objective, rank, evidence, trace, attrition and workflow in the Workbench; verify and surface provenance in reports and RDML, including altered-manifest warnings. |
| RF-28 | Report thermodynamics used mutable session conditions | An assay designed at one Mg/salt/oligo/annealing condition could be exported later with Tm and structure recomputed under another user's current settings. | Recompute each assay under its self-hashed manifest snapshot; label legacy/current-session fallback explicitly and export the condition basis in CSV/HTML. |
| RF-29 | Direct manual analysis was a dead-end for selection | The exact entered assay could be fully mapped and ranked but could not enter the Workbench unless the user first ran redesign or rescue, encouraging an unnecessary sequence change and breaking the intended direct-entry workflow. | Treat the exact analyzed assay as a first-class candidate and transfer its sequences, evidence, objective, trace, condition snapshot, workflow, and manifest directly into the Workbench. |

## Candidate attrition ledger

Every design run can report:

- stage name, schema/version, unit of analysis, entered/retained/rejected counts;
- exact hard-screen and heuristic-retention reasons;
- reversibility and candidate decisions;
- windows attempted/completed and runtime-budget status;
- duplicate, regional, pair, probe, and annotation-budget decisions;
- retained candidate identities and exact coordinates;
- objective, ranker/profile/model versions, conditions, input hashes, and candidate limits;
- evaluations performed and fallbacks used.

A downstream rank cannot recover an earlier discarded candidate, so the ledger explicitly identifies every such irreversible fact and every reversible engineering truncation.

## Before-and-after trace: frozen Plasmodium/Haemoproteus fixture

Input: 11 *Plasmodium* cytb targets, 12 *Haemoproteus* off-targets, discrimination objective, `parasite_mtdna` profile.

Current trace:

- 18/18 windows evaluated.
- 157,608 primer candidates entered hard screening; 33,076 survived.
- 15,193,998 possible pair relationships entered geometry/Tm screening; 1,404,393 survived.
- Dimer/native-call cap retained 5,400 pairs; pair-diversity beams retained 216.
- 297,688 probe candidates entered hard screening; 7,874 survived; probe-diversity beams retained 414.
- 328 complete triplets survived window beams; base diversity retention kept 22.
- 15 discrimination-specialist triplets and 19 target/off-target-aware probe alternatives were added without bypassing final hard gates.
- 56 candidates entered the final annotation-retention stage; 25 diverse complete assays received full annotation.
- Five finalists were returned in 30.94 seconds on the audit machine.
- Rank 1 modeled 100% coherent target/probe coverage, zero supplied off-target products, and full validity across all three declared condition scenarios.

The complete machine-readable trace is `tests/benchmark/plasmodium_ranking_trace.json`.

## Remaining limitations

- Search is bounded and can miss an unretained triplet.
- External database completeness, accession versions, and user-supplied labeling constrain specificity evidence.
- Three condition scenarios are a robustness screen, not exhaustive chemistry validation.
- Degenerate order pools, modified probes, amplicon structure, fluorescence, inhibition, reagent competition, synthesis behavior, and biological matrix effects remain imperfectly modeled.
- The synthetic ranking benchmark validates software decision behavior, not wet-lab superiority.
- No computational rank establishes a universal “best assay.”

## v1.34.0 decision-analysis addendum

### Confirmed residual defects

29. **Failed 3′ placements were hidden in Manual Design.** The default primer-site scan anchored on an exact 3′ base and therefore omitted near-match placements that explained why an existing primer would not extend. The manual inventory now scans all allowed near-matches, reports exact/ambiguous/mismatching 3′ states, and separately marks extension eligibility. Product construction still accepts only eligible placements.
30. **Rank sensitivity was generic rather than evidence-linked.** Explanations could say that a rank might reverse without naming the measurable trade-off. The explanation layer now emits concrete triggers tied to target coverage, supplied off-target products, condition robustness, practical burden, and missing objective-required evidence.
31. **Two design runs could not be forensically compared.** A changed winner could not be separated into a legitimate input/objective/model change versus unexplained non-reproducibility. `run_compare` now verifies manifests, normalizes candidate identity, compares contexts, measures top-k overlap/Spearman/pair reversals, and classifies reproducibility.
32. **Bench feedback lacked an assay-local evidence summary.** Records could be stored but not consolidated without inviting over-generalization. The feedback layer now summarizes each assay only within its declared contexts and creates pairwise preferences only for unanimous success-versus-failure comparisons under the same target-group/condition context.
33. **Batch Design bypassed the authoritative run object.** It called the convenience single-assay function and returned bare oligos, losing objective, full retained-pool evidence, attrition, uncertainty, and manifest provenance. Batch Design now calls `design_from_sequences`, retains the winner's evidence and run manifest, and declares a bounded per-template search budget.
34. **Manual sequence edits lacked complete before/after diagnosis.** Users could re-run analysis but could not see which hard failures were resolved or introduced, which metrics improved or worsened, or whether the edited assay still mapped uniquely. `compare_edits` now recalculates both versions under an identical context and returns an auditable delta.

### Before-and-after traces

- **3′ mismatch:** before, a one-base terminal mismatch produced no manual placement row; after, the row is present with `three_prime_status=mismatch` and `extension_eligible=false`.
- **Batch winner:** before, `/api/batch_design` returned the first convenience-path winner with no run manifest; after, the result reports `pipeline=authoritative_structured_ranker`, hard validity, explanation, candidate attrition, alternatives evaluated, and a self-hashing manifest.
- **Manual edit:** before, changing one nucleotide required manually comparing two unrelated result blocks; after, the edit comparator reports sequence operations, resolved/new failures, objective evidence deltas, and the conditional preferred assay.
- **Run reproducibility:** before, identical manifests with reordered candidates were not diagnosed; after, the comparison classifies this as `critical_non_reproducibility`.

The new functions do not prove biological superiority. They make candidate loss, manual edits, rank changes, and evidence gaps reconstructable.
