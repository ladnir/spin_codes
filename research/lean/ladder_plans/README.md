# Turbo s22 ladder plans (regeneration inputs)

Per rung `g{e}` (K = 2^e): `config.json` (geometry, q1 tilts, q1/q2 bits from
`scripts/ladder22.py`) and `tail_plan.json` (the tail pieces: exact pieces and marker
ranges with tilt, alpha, t, q range) from the untrusted planners
(`scripts/ladder_exact_plan22.py`, `scripts/ladder_marker_plan22.py`).

The generated Lean (`SpinCodes/Turbo/Generated/`, `SpinCodes/Turbo/PhysicalProfile*G*Tilt*`,
`PhysicalSmallBirth22G*Tilt*`) is NOT committed; regenerate it from these plans:

    cp ladder_plans/g{e}/*.json .lake/turbo/ladder/g{e}/
    python scripts/ladder22.py {e} config families q1 q2
    TURBO_QCHUNK_EXACT=4 python scripts/ladder_tail22.py {e} families pairs exact marker emit endpoint

then build `SpinCodes/Turbo/Generated/G{e}/Endpoint` on Peach (`pbuild.py`). Scoped build
limits beyond `-M 14000` (G24 TailP15MQ0; G26 TailP18-21MQ0) are listed in
`REVISED_SPIN_PROOF_PLAN.md`. Use `TURBO_QCHUNK_EXACT=6` for G22, `30` (default) below.

Certified endpoints (`Spin.Turbo.G{e}Final.g{e}_endpoint`, standard axioms):
2^14 64.25152, 2^16 62.34042, 2^18 60.36266, 2^20 58.36785, 2^22 56.3691,
2^24 54.36978, 2^26 52.32161 bits.
