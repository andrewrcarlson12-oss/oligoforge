# Paper-readiness assessment — OligoForge 1.36.0

## 1.36.0 impact

The first-class Validation Studio and Assurance browser workspaces improve usability, evidence visibility, and workflow reproducibility. They do not add a sufficiently large, target-grouped, held-out wet-lab dataset and therefore do not change the paper's biological-validity ceiling. Any manuscript must describe candidate selection and DriftGuard as bounded computational evidence and keep experimental performance claims separate.

The software methods, candidate provenance, attrition accounting, explicit rank uncertainty, benchmark-integrity checks, and adversarial computational validation are suitable for a methods preprint after independent code review. Claims must remain limited to computational behavior.

A paper claiming superior assay performance requires a preregistered, target-group-separated wet-lab comparison of rank 1 against diverse alternatives and external design tools under matched inputs. Preserve failed assays, adjudicate conflicting outcomes, freeze development/tuning/held-out/final-test splits, report all exclusions and candidate limits, and retain versioned run manifests.

Do not use the current synthetic/adversarial benchmark or the minimum feedback-data gate as evidence of biological superiority or a validated learned reranker.

## 1.34.0 impact

Decision traceability is stronger: manual edits, batch winners, run-to-run changes, and local experimental feedback are now auditable. This improves methods reporting and reproducibility, but the central publication blocker remains unchanged: no sufficiently large, target-group-separated wet-lab ranking dataset demonstrates that rank 1 outperforms lower-ranked candidates across independent assay families.
