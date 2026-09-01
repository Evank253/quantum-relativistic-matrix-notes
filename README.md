# Quantum-relativistic matrix notes + Quantara 0.1

Working notes on wormholes, chronology protection, holography, JT Page curves, Gao-Jafferis-Wall, and a small Python lab that computes the formulas from the session.

These are study notes, not a claim that vacuum energy, dark matter, or gravity have been solved. Not a SaaS AGI stack.

## Repo

https://github.com/Evank253/quantum-relativistic-matrix-notes

Default branch: `Python-3`

## Quantara commands (from the package in this chat)

```bash
cd quantara
PYTHONPATH=. python -m quantara.cli qed
PYTHONPATH=. python -m quantara.cli page --S 100
PYTHONPATH=. python -m quantara.cli cylinder
PYTHONPATH=. python -m quantara.cli cylramp --t 20 --beta 2
PYTHONPATH=. python -m quantara.cli chords --q 0 --nmax 6
python -m pytest tests -q
```

### Checks that should hold

- `t_Page = 600` for `S=100`, `c=1`, `beta=2π`
- `Z_{0,2}` numeric = analytic to ~1e-15
- `Z_{0,2}(β+it,β-it) ~ t/(4πβ)` at large t (ramp)
- chord moments at `q=0` are Catalan numbers; at `q=1` they are double factorials
