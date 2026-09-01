# New SPIN manuscript

The new manuscript entry point is `main.tex`. Build from this directory:

```text
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The pre-restart root-level TeX draft is preserved in
`../old/paper_draft_pre_spin_2026-08-31/`. Its `ARCHIVE_MANIFEST.md` records
the moved files, byte lengths, and SHA-256 hashes.

Red `TODO` paragraphs mark unresolved definitions, integrations, and theorem
inputs. They are intentionally visible in the draft and must not be read as
claims.

The draft was last built successfully with TeX Live 2026 on 2026-09-01. The
verified 34-page PDF is written to
../output/pdf/spin_codes_draft.pdf.
