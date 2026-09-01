"""Crossing-weighted chord moments for double-scaled SYK.

[n]_q = (1 - q^n) / (1 - q)
H |n> = sqrt([n]_q) |n-1> + sqrt([n+1]_q) |n+1>
m_{2n}(q) = <0| H^{2n} |0>

q -> 0: Catalan numbers (noncrossing).
q -> 1: (2n-1)!! (all matchings).
"""
from __future__ import annotations

from functools import lru_cache
from math import comb


def q_number(n, q):
    if n == 0:
        return 0.0
    if abs(q - 1.0) < 1e-15:
        return float(n)
    return (1.0 - q**n) / (1.0 - q)


def double_factorial_odd(n):
    acc = 1
    for k in range(n):
        acc *= 2 * k + 1
    return acc


def catalan(n):
    return comb(2 * n, n) // (n + 1)


@lru_cache(maxsize=None)
def _amp(steps, level, q):
    if level < 0 or level > steps:
        return 0.0
    if steps == 0:
        return 1.0 if level == 0 else 0.0
    down = q_number(level + 1, q) ** 0.5 * _amp(steps - 1, level + 1, q)
    up = q_number(level, q) ** 0.5 * _amp(steps - 1, level - 1, q)
    return down + up


def moment(n_even, q):
    if n_even < 0:
        raise ValueError("n >= 0")
    if n_even % 2 == 1:
        return 0.0
    return _amp(n_even, 0, q)


def moment_table(n_max, q):
    rows = []
    for n in range(0, n_max + 1, 2):
        m = moment(n, q)
        rows.append({
            "n": n,
            "m_n": m,
            "matchings": double_factorial_odd(n // 2) if n else 1,
            "noncrossing": catalan(n // 2) if n else 1,
        })
    return {"q": q, "moments": rows}


def crossing_weight_limits(n):
    if n % 2:
        raise ValueError("even n")
    return {
        "n": n,
        "q0_catalan": moment(n, 0.0),
        "catalan": float(catalan(n // 2)),
        "q1_all_matchings": moment(n, 1.0),
        "double_factorial": float(double_factorial_odd(n // 2)),
    }
