# Proposed paper update: formal verification of distance and rate

Review date: 2026-09-28. **Applied after user approval.** Both manuscript
variants build and the edited pages pass visual review. The theorem statement
is unchanged; see `scripts/map_data/paper_revision_verification.json`.
The review covers the scalable theorem and its proof, abstract/introduction
framing, certificate discussion, artifact statement, and AI disclosure.
It is not a fresh review of the finite BCH or application-security claims.

## Recommendation

Keep the theorem, construction, 11% threshold, and rate unchanged. Use the
paper's existing terminology: say that the theorem's distance and rate claims
are formally verified in Lean. Do not introduce “native-length” or
“exact-rate” labels; the theorem already specifies the lengths and rate.

Explain the scope once in the formalization appendix. Put differences in
intermediate bounds there too. Keep verification commands and replay status
in the artifact documentation. No new mathematical error was found in the
reviewed argument.

## 1. Abstract: make the formal result visible

At [abstract.tex](../paper/abstract.tex), after the sentence ending
“structured support spreading and recursive state mixing” (line 11), add:

```tex
We formally verify these distance and rate guarantees in Lean.
```

This sentence attaches verification to the mathematical clauses that were
checked. It does not extend the claim to encoding work, finite BCH parameters,
or application security. Preserve the current theorem and performance claims.

## 2. Introduction: add a short contribution statement

In [introduction.tex](../paper/introduction.tex), insert the following after
the asymptotic encoding-work sentence (line 36), before “For finite lengths”:

```tex
The distance and rate claims are formally verified in Lean, including
the construction's random setup, the analytic estimates, and the required
numerical certificates. Appendix~\ref{app:structured-lean} describes the
formalization.
```

This placement makes the reference to distance and rate unambiguous before
the paragraph moves to the separate finite-length results. No additional
contribution heading or repetition of the theorem is needed.

## 3. Theorem: retain its statement and identify the checked clauses

Immediately after Theorem `thm:structured-spin-scalable` in
[structured_proof.tex](../paper/structured_proof.tex) (line 91), insert:

```tex
The distance and rate claims of Theorem~\ref{thm:structured-spin-scalable}
are formally verified in Lean; see Appendix~\ref{app:structured-lean}.
```

Keep the theorem, rate/work proof, and requested-length argument in place.
There is no need to explain Lean's zero-based index here: it is merely a
shift of the paper's positive index and does not change the construction.

Also clarify the first sentence of the final proof (line 446), which currently
says “Condition on any outer satisfying” the selection event. Replace it with:

```tex
Fix a shared outer realization in $\mathcal G_b$. The route--inner bounds
hold uniformly over this realization. Averaging under the outer law
conditioned on $\mathcal G_b$ gives the conditional first-moment bounds below.
```

The later union bound already charges selection failure separately. This edit
distinguishes a fixed realization from conditioning on an event; it does not
repair an independence error.

## 4. Appendix: distinguish the formal proof from numerical replay

Extend the existing “Completion and Certificate Boundary” discussion at
[structured_imt_appendix.tex](../paper/structured_imt_appendix.tex), lines
466–485. Retain the distinction between high-precision numerical replay and
proof of the analytic steps. Add a paragraph/subsubsection with label
`app:structured-lean` and the following content:

```tex
\subsubsection{Formal Verification}
\label{app:structured-lean}
The Lean development defines the binary encoder and setup distribution
used in Theorem~\ref{thm:structured-spin-scalable}. It proves the distance
statement \eqref{eq:structured-spin-distance} and rate $1/2$ for every setup
realization. The proof includes the shared outer spectrum, the routed IMT
process, the three occupation regimes, and the final first-moment argument.
Lean checks both the analytic arguments and the required numerical
inequalities. The encoding-work claim, subsequent padding corollary, and
separate finite-length BCH results are outside the final checked statement.
```

Then explain the differences without suggesting that the sharper paper
bounds were disproved:

- The dense assembly uses margin `4×10^-7`, smaller than the reported
  `4.10335×10^-7`. All 1,023 indexed boxes are covered.
- For each fixed positive occupation `Q`, the formalization obtains the
  following bounds for all sufficiently large family indices. The threshold
  index may depend on `Q`:

  | Occupation | Bound on the conditional first moment |
  | --- | --- |
  | `Q=1` | `1000 exp(-b/50)` |
  | `Q=2` | `1000 exp(-b/20)` |
  | Fixed `Q≥3` | `(1600/3) exp(-Qb/2000)` |

- For `Q=1,2`, Lean uses different tilts. For larger fixed occupations, it
  reserves a fixed positive margin in the continuum transfer. These bounds
  suffice for the same limit but do not establish every sharper intermediate
  estimate displayed in the paper. Only `1≤Q<4096` is summed using fixed-Q
  convergence; the growing sparse range has its own uniform bound.

A short concluding paragraph should identify `SpinCodes.Native` and
`SpinCodes.NativePin`, state that the recursive axiom audit reports only
`propext`, `Classical.choice`, and `Quot.sound`, and direct readers to the
artifact's reproduction guide. Keep the account of reused dependencies and
the separate full-source replay in that guide. The paper need not recount
the build process or make a claim about a fresh replay.

**Recommended treatment of constants:** retain the sharper paper estimates
with this explicit correspondence note. Replacing them throughout would be
an alternate proof presentation, requiring coordinated changes to the tilt
choices, matrix estimates, and supporting certificate references. It is not
necessary to state the verified final result accurately.

## 5. Artifact description: make the proof reproducible

The introduction currently advertises only the core encoder artifact
(lines 149–154). Expand both the anonymous and public branches to include
Lean sources, generated certificates, and reproduction instructions when
those files are included in the distributed artifact. Use `\repositorycite`
or the existing submission conditional to preserve anonymity.

Identify a versioned source snapshot for publication. The new Lean files
are currently untracked in this working tree; the existing public repository
citation alone does not establish their availability. Packaging/publication
is a separate next action, not something this review has performed.

Keep toolchain pins, commands, provenance categories, and fresh-replay status
in [FINAL_REPRODUCTION.md](FINAL_REPRODUCTION.md), rather than putting local
paths, compiler job counts, or mutable progress numbers in the manuscript.
A versioned artifact report can record the precise verification run.

## 6. AI disclosure: include formal proof development

[ai_disclosure.tex](../paper/ai_disclosure.tex), lines 4–12, describes the
manuscript and implementation work but omits the Lean development. Add:

```tex
AI tools also assisted the Lean proof development, numerical-certificate
generation, and verification orchestration. The Lean development included
Claude-assisted work and GPT-assisted completion.
```

The conversation establishes Claude's participation, but not its exact
version. Preserve the authors' responsibility statement. Keep the technical
trust explanation in the formalization appendix rather than duplicating it
here. The abstract footnote currently lists only GPT models; simplify it to:

```tex
Generative AI tools were used extensively in this work. See
Section~\ref{sec:ai-disclosure} for details.
```

## Integration checks after editing

Check the submission and public builds, all new labels/references, and each
use of “verified” or “formalized” for scope. Keep the mathematical theorem
unchanged. Recheck numerical constants against the correspondence note.
Changing paper sources intentionally changes the paper hashes preserved by
the earlier Lean closure audit: keep that audit as a historical snapshot and
record the new manuscript version separately; do not relabel its old hashes.

Next: apply this focused edit set, compile both manuscript variants, and
review the resulting pages. Full source replay and publication of a versioned
Lean artifact remain separate verification/release tasks.

The short conditioning rewrite is a clarity improvement. The artifact and AI
disclosure changes update the record of the completed work. Neither changes
the theorem or introduces new terminology for its existing assumptions.
