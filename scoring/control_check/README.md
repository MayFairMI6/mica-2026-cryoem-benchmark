# Recorded control checks

In the archived `metrics.json` files and the failed combined-control `.log`,
`<WORKSPACE>/mica_benchmark` replaces the original machine-specific temporary
workspace prefix. It is a path placeholder, not a literal directory or a claim
that these checks were rerun. Relative paths and filenames are retained.

This normalization changes only the path prefix. Scores, input SHA-256 hashes,
control definitions, commands, timestamps, and diagnostic messages are retained.
The original records remain available in Git history.
