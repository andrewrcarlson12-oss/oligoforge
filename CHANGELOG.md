# Changelog

## 1.31.1 — scientific-correctness and release-hardening audit

### Design and specificity
- Search the full target in a bounded spread of windows, collect candidates globally, and rank before trimming.
- Preserve exact primer/probe/amplicon coordinates through windowed search.
- Reject unresolved template bases in oligo windows instead of resolving them arbitrarily.
- Evaluate ambiguous subject bases conservatively in offline specificity and scan all viable primer sites.
- Rank multi-isolate assays using coherent complete-amplicon coverage and penalize impractical degeneracy.

### Thermodynamics and multiplexing
- Clarify model/concentration conventions and remove claims of exact vendor equivalence.
- Preserve explicit modified-probe order notation while separating DNA-backbone and LNA-model estimates.
- Strictly reject malformed multiplex oligos; compare assays by identity; report dimer ΔG at 37 °C and anneal, dimer Tm, and 3′ engagement.
- Make orthogonal-panel reaction-condition overrides request-local.
- Use exact branch-and-bound or rigorous clique-cover bounds for model certificates; theta remains diagnostic only.

### Quantitative analysis and reporting
- Require sustained amplification and baseline detrending for raw Cq screening.
- Regress standard curves on dilution-level means rather than treating technical replicates as independent levels.
- Scope reference-gene output to geNorm-style M/V plus a Cq-SD screen.
- Prevent computational-only records from passing MIQE-aligned readiness checks.
- Escape report fields, neutralize spreadsheet formulas, validate order notation, and reject unresolved gBlock bases.
- Upgrade assay-definition export to RDML 1.3 without guessing target roles or reporter dyes.

### Hosted reliability and security
- Disable shared persistence and shared condition mutation by default in hosted mode.
- Block client-selected local BLAST paths on hosted deployments.
- Bound NCBI timeouts/retries and request sizes.
- Sanitize validation and exception responses; add browser security headers.
- Pin runtime dependencies and run the container as a non-root user.

### Compatibility
- Version synchronized to 1.31.1.
- Existing deterministic scientific fixtures are retained unless a test encoded a superseded or scientifically invalid assumption.
