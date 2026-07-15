# Contributing

1. Create a focused branch and keep scientific and UI changes reviewable.
2. Add a regression that fails before the repair and passes afterward.
3. Run `python run_tests.py` before submitting changes.
4. Do not update locked assay fixtures merely to make a changed algorithm pass; explain and review any intentional scientific-output change.
5. Keep computational predictions, empirical observations, and validation claims explicitly separated.
6. Reject malformed or ambiguous inputs rather than silently rewriting them into different oligos.
7. Never expose credentials, local paths, raw exception text, or user sequence payloads in hosted responses or logs committed to the repository.
8. Synchronize the release version in `oligoforge/__init__.py`, `app.py`, `launcher.py`, and `static/index.html`.
