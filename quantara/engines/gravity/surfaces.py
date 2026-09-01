"""Disk and trumpet integrals for JT gravity (SSS conventions used in-session).

Trumpet:
    Z_tr(β, b) = exp(-b² / 4β) / sqrt(4π β)

Cylinder (double trumpet):
    Z_{0,2}(β1, β2) = ∫_0^∞ b db  Z_tr(β1,b) Z_tr(β2,b)
                    = sqrt(β1 β2) / (2π (β1 + β2))

Continued to β ± it the cylinder grows ~ t, which is the SFF ramp.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.integrate import quad


def trumpet(beta, b):
    if b < 0:
        raise ValueError("geodesic length b >= 0")
    return np.exp(-(b * b) / (4.0 * beta)) / np.sqrt(4.0 * np.pi * beta)


def disk_schwarzian(beta, C=1.0, S0=0.0):
    if beta <= 0:
        raise ValueError("beta > 0")
    return (math.pi / beta) ** 1.5 * math.exp(S0 + 2.0 * math.pi**2 * C / beta)


def cylinder_analytic(beta1, beta2):
    return np.sqrt(beta1 * beta2) / (2.0 * np.pi * (beta1 + beta2))


def cylinder_numeric(beta1, beta2, bmax=80.0):
    if beta1 <= 0 or beta2 <= 0:
        raise ValueError("use real positive beta for the numeric check")

    def integrand(b):
        return b * float(np.real(trumpet(beta1, b) * trumpet(beta2, b)))

    val, err = quad(integrand, 0.0, bmax, epsabs=1e-10, limit=200)
    closed = complex(cylinder_analytic(beta1, beta2)).real
    return {
        "numeric": float(val),
        "analytic": float(closed),
        "abs_err_quad": float(err),
        "rel_err_vs_closed": float(abs(val - closed) / max(abs(closed), 1e-16)),
        "beta1": beta1,
        "beta2": beta2,
    }


def cylinder_sff_ramp(t, beta):
    if beta <= 0:
        raise ValueError("beta > 0")
    b1 = complex(beta, t)
    b2 = complex(beta, -t)
    z = cylinder_analytic(b1, b2)
    return {
        "Z02": {"real": z.real, "imag": z.imag},
        "abs": abs(z),
        "linear_asymptote": abs(t) / (4.0 * math.pi * beta),
        "exact_abs": math.sqrt(beta * beta + t * t) / (4.0 * math.pi * beta),
        "t": t,
        "beta": beta,
        "shape": "ramp" if abs(t) > beta else "near-zero",
    }
