# OligoForge tests — run all from the repo root

    OLIGOFORGE_EMAIL=you@x.com python3 tests/test_regression.py   # offline science + wiring asserts
    node tests/ui_handlers.js     # 32 tab-handler renders (stubbed DOM)
    node tests/ui_workbench.js    # workbench card flow: add -> edit -> check -> attach std -> report -> multiplex -> remove
    node tests/ui_conditions.js   # reaction conditions init/apply + Load-example(FSJ) seed
    node tests/ui_projects.js     # server-side project save/list/load/delete

The UI harnesses read static/index.html relative to the current directory, so run them from the repo root.
Re-run ALL of these after any edit, then rebuild the zip from this folder before shipping.
