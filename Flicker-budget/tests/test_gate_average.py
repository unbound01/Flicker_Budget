"""Eq. (8) is the first-order worst-state infidelity of a pi rotation with LO frequency noise, and the state average is 2/3 of it."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import gate_two_level as G  # noqa: E402


def test_first_order_weight_equals_gong_filter():
    for u in (0.01, 0.3, 0.7, 1.0001, 1.5, 2.2, 7.3):
        assert abs(G.first_order_weight(u) / G.gong_weight(u) - 1) < 1e-3
    assert G.first_order_weight(3.0) < 1e-12 and G.gong_weight(3.0) < 1e-12      # notch at odd u >= 3


def test_state_average_is_two_thirds_of_worst_case():
    rng = np.random.default_rng(7)
    for _ in range(12):
        U = G.evolve(rng.uniform(0.1, 4), rng.uniform(0.01, 0.05), rng.uniform(0, 2 * np.pi), n=200)
        assert abs(G.avg_over_worst(U, grid=25) - 2 / 3) < 3e-3
