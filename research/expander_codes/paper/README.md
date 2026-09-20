# Unified manuscript

The manuscript presents one dependency-ordered argument:

1. exact enumerators for the original binary EA and EC ensembles;
2. regular expander ensembles and their finite-parameter improvements;
3. the labeled finite-field generalization;
4. certified finite parameters and the degree--memory frontier.

Material in `../notes/previous_draft/` can be reused only after it has been
placed into this unified structure and checked against the new notation.

Build from this directory with:

```powershell
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

`CWC_REVIEW.md` records the section-by-section controlled-writing pass.  The
submission source contains no visible editorial placeholders.  From the parent
directory, `python scripts/verify_manifest.py --strict-versions --list`
checks the complete claim-to-artifact manifest.  Use `--run headline` or
`--run all` to execute rigorous entries serially.
