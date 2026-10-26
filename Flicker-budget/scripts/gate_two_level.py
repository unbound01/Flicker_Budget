#!/usr/bin/env python3
"""Independent check of the gate criterion, Eq. (8) (Gong et al. JSSC 2023, Eqs. 1-2), and of the 2/3 state-average conversion.

Model: a resonant pi rotation at Rabi frequency f_R driven by an LO with frequency noise df(t). In the frame of the LO the
qubit sees H(t) = (Omega/2) sigma_x + (delta(t)/2) sigma_z, with Omega = 2 pi f_R and delta = 2 pi df.

(1) first_order_weight(u): to first order the error is a rotation exp(-i(a sigma_y + b sigma_z)) with
    a = int (delta/2) cos(Omega t) dt, b = int (delta/2) sin(Omega t) dt over the pulse T = pi/Omega.
    Its worst-placed input state has infidelity a^2 + b^2. For a frequency tone of variance s^2 at f = u f_R, averaged over
    phase, E[a^2 + b^2] / (s^2 / f_R^2) is the weight on S_f / f_R^2; it should equal |H_LO(u)|^2 of Eq. (8).
(2) avg_over_worst(...): exact two-level evolution (piecewise-constant expm), average gate fidelity
    (2 + |Tr(U_id^dag U)|^2) / 6 versus the worst input state found on a Bloch-sphere grid. For any single-qubit
    unitary error the ratio of (1 - F_avg) to (1 - F_worst) is 2/3.

Run: python scripts/gate_two_level.py
"""
import numpy as np
from scipy.linalg import expm

SX = np.array([[0, 1], [1, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def gong_weight(u):
    """|H_LO(u)|^2 of Eq. (8)."""
    return (1 + np.cos(np.pi * u)) * (1 + u**2) / (2 * (1 - u**2) ** 2)


def first_order_weight(u, n_t=4001, n_phase=64):
    fr = 1.0
    w = 2 * np.pi * fr
    t = np.linspace(0, np.pi / w, n_t)
    out = []
    for th in np.linspace(0, 2 * np.pi, n_phase, endpoint=False):
        d = 2 * np.pi * np.cos(2 * np.pi * u * fr * t + th)          # unit-amplitude frequency tone, variance 1/2
        a = np.trapezoid(d / 2 * np.cos(w * t), t)
        b = np.trapezoid(d / 2 * np.sin(w * t), t)
        out.append(a * a + b * b)
    return np.mean(out) / (0.5 / fr**2)


def evolve(u, amp, th, n=400, fr=1.0):
    w = 2 * np.pi * fr
    T = np.pi / w
    dt = T / n
    U = np.eye(2, dtype=complex)
    for i in range(n):
        tm = (i + 0.5) * dt
        d = 2 * np.pi * amp * np.cos(2 * np.pi * u * fr * tm + th)
        U = expm(-1j * dt * (w / 2 * SX + d / 2 * SZ)) @ U
    return U


def avg_over_worst(U, grid=31):
    E = expm(-1j * np.pi / 2 * SX).conj().T @ U
    f_avg = (2 + abs(np.trace(E)) ** 2) / 6
    worst = 1.0
    for a in np.linspace(0, np.pi, grid):
        for p in np.linspace(0, 2 * np.pi, grid):
            s = np.array([np.cos(a / 2), np.exp(1j * p) * np.sin(a / 2)])
            worst = min(worst, abs(s.conj() @ E @ s) ** 2)
    return (1 - f_avg) / (1 - worst)


def main():
    print(" u        |H_LO|^2 (Eq. 8)   first-order worst-state weight")
    for u in (0.01, 0.3, 0.7, 1.0001, 1.5, 2.2, 3.0, 7.3):
        print(f"{u:7.4f}   {gong_weight(u):12.5g}   {first_order_weight(u):12.5g}")
    rng = np.random.default_rng(1)
    r = [avg_over_worst(evolve(rng.uniform(0.1, 4), rng.uniform(0.01, 0.05), rng.uniform(0, 2 * np.pi))) for _ in range(200)]
    print(f"exact two-level, 200 random tones: (1-F_avg)/(1-F_worst) min {min(r):.4f} max {max(r):.4f}")


if __name__ == "__main__":
    main()
