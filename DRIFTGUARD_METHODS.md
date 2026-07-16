# DriftGuard methods

DriftGuard compares an immutable baseline with a follow-up snapshot. For each registered complete assay and unique supplied sequence, it uses the existing OligoForge isolate engine to reconstruct convergent primer products and probe recognition. It reports new target product loss, new probe-recognition loss and new signal-capable off-target products with exact record hashes and affected component identifiers.

States are reason-coded: Stable, Stable with new variation, Watch, Redundancy degraded, Possible target dropout, Possible signal-generating off-target, Interpretation ambiguity, Action review recommended, Evidence insufficient or Scan incomplete. The initial vertical slice emits only states justified by its implemented evidence; it does not manufacture a numeric risk score.
