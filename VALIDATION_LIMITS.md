# Validation limits — OligoForge 1.37.0

OligoForge produces computational candidates, structured evidence, uncertainty labels, provenance manifests, and assay-readiness records. It does not biologically validate an assay.

Before describing an assay as validated or publication-ready, obtain representative target and near-neighbor panels; replicate standard curves; efficiency and linearity; adequately replicated LOD/LOQ; product identity; inclusivity/exclusivity; inhibition and matrix tests; probe-signal performance; multiplex interaction tests; synthesis/vendor review for modified oligos; and locked analysis criteria.

The ranking benchmark proves that the software obeys frozen adversarial decision preferences and avoids known compensation and candidate-loss failures. It does **not** prove that rank 1 will outperform ranks 2–20 at the bench. The repository does not contain a sufficiently large target-grouped held-out experimental preference dataset.

Rank confidence is conditional on supplied evidence. Missing off-target, panel, junction, sequence-version, or database-version information can produce an explicit insufficient-evidence state. “Near-equivalent” means the implemented computational evidence does not support a decisive distinction; it is not an equivalence proof.

Search is heuristic-bounded. A stronger triplet can exist outside the retained pool. Every run records candidate limits, attrition, conditions, input hashes, model versions, external-database state, warnings, and fallbacks so that this uncertainty remains visible.

The canonical design contract verifies which shared policy ran, whether objective prerequisites were supplied, whether the winner cleared recorded computational hard gates, and whether two results are configuration-comparable. It does not certify database completeness, synthesis, amplification, fluorescence, analytical performance, or wet-lab reproducibility.

Experimental-feedback records remain local evidence. The learned-reranker gate is only a minimum eligibility screen and does not establish that a model is valid or deployable.

Manual edit comparison is also computational. A reported improvement means that the edited assay has stronger modeled evidence under the supplied context; it does not establish improved amplification efficiency, fluorescence, inhibition tolerance, synthesis quality, or multiplex behavior. Batch designs use a shorter declared search budget than a full interactive design and can therefore miss candidates retained by a broader run; compare the manifests and rerun critical targets interactively before ordering.

Validation Studio chooses among user-supplied cases that distinguish supplied candidates under the declared complete-product model. Its selection is bounded and deterministic, not proof of global experimental optimality. Its interpretation is valid only for the declared samples, controls, reaction conditions, observations, and acceptance criteria. Invalid controls suppress a candidate conclusion.

AssaySBOM and sequence-snapshot hashes establish record integrity, not scientific truth. DriftGuard reports sequence-level modeled concerns over supplied snapshots; it does not estimate the probability of assay failure, future evolution, clinical impact, or regulatory acceptability. Snapshot sampling and metadata can be incomplete or unrepresentative. See `ASSURANCE_VALIDATION_LIMITS.md`.
