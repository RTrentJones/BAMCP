"""Pre-registered visual-vs-text eval: do pileup renders improve model analysis over text?

Design and analysis plan: ``docs/eval/visual-preregistration.md``. Modules:

- :mod:`.cases` / :mod:`.fixtures` — case types, labels and the deterministic ``vt_eval.bam``
- :mod:`.stimuli` — opaque per-(item, condition) stimulus directories + ``runs.jsonl``
- :mod:`.leaklint` — fails a run set whose stimuli could reveal labels
- :mod:`.answers` — answer schema / validation; :mod:`.grade` + :mod:`.stats` — metrics
- :mod:`.charts` — PNG charts; :mod:`.selftest` — oracle / random grader self-tests
- :mod:`.cli` — ``python -m bamcp.eval.visual {fixtures,build,lint,selftest,grade,dry}``
"""
