# Paper and research

Start with the [SPIN manuscript](paper/README.md). The production library is
separate, in [`../spin/`](../spin/README.md), and does not depend on this tree.

| Directory | Contents |
|---|---|
| [paper/](paper/README.md) | Current manuscript, figures, and supporting paper scripts. |
| [artifact/](artifact/README.md) | Paper artifact packaging and author-side reproduction tools. |
| [workstreams/](workstreams/) | Proof development, experimental implementations, and selected results. |
| `scripts/`, `constructions/`, `explorations/`, `bch_spectrum_work/` | Supporting research and historical experiments. |
| `BA_paper/`, `enumerator_paper/`, `expander_codes/` | Related manuscripts. |
| `collaboration/` and top-level Markdown notes | Historical plans, audits, and handoffs. |

## Running historical commands

The research tree was moved together, preserving its internal relative paths
and the bytes of frozen sources and manifests. Run commands in its historical
notes from this directory, not the repository root:

```sh
cd research
python -B -m unittest discover -s artifact -p 'test_*.py'
```

Old absolute checkout paths and generated build directories may need to be
reconfigured. Historical GitHub links using paths at the former repository root
can be read at commit `557e40119f58323ad4e10238532d131be0a17f02`; current paths
have the `research/` prefix. Frozen manifests retain their original relative
names and hashes.

See the [publication policy](GITHUB_PUBLISH_POLICY.md) before adding research data.
Do not commit experiment dumps or build products.
