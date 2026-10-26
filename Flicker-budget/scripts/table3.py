#!/usr/bin/env python3
"""Every number in Table III and Sections V-C, VI and VII of the manuscript, from the digitized data.

Writes
  <out>/table3_rows.tex   the body rows of Table III (\\input by main.tex)
  <out>/numbers.json      all values quoted in the text, keyed by name

Run:  python scripts/table3.py [--out paper/generated]
Needs the digitized data in data/ (regenerate with the digitize_* scripts; verify with data_checksums.py).

Conventions (see docs/filter_criteria.md): S_f = 2 f^2 L; gate = Gong Eq. 1-2, pi rotation, worst-case |H_LO|^2;
idle at tau = 120 us: Ramsey var = (2 pi tau)^2 int S_f sinc^2, echo var = int S_phi 16 sin^4(pi f tau/2);
idle infidelity average (1-e^{-v/2})/3, worst case (1-e^{-v/2})/2. Integration band 1 Hz - 100 MHz unless stated.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
import budget as B  # noqa: E402
import filter_criteria as F  # noqa: E402
from a0_gong_measured import h_lo2  # noqa: E402

TAU = 1.2e-4
TARGET = 1e-5
FRS = (0.5e6, 1e6, 2e6, 4e6, 10e6)
SIGMA_Q = 1.9e3
D = lambda n: json.load(open(os.path.join(ROOT, "data", n + ".json")))


def points():
    """Nine 4-K operating points: (key, cite key, own carrier GHz, trace, dB offset to own carrier)."""
    out = [("G7-10.0", "gong2023", 10.0, np.array(D("gong_fig13b")["trace"]), F.DB_DIV4)]
    for n, p in D("gong_thesis_fig721")["panels"].items():
        out.append((f"G10-{n[:-3]}", "gongthesis", float(n[:-3]), np.array(p["trace"]), 20 * np.log10(8)))
    out.insert(0, ("X12-8.6", "xiao2025", 8.6, np.array(D("xiao_rfic2025_fig7")["trace"]), 0.0))
    for n, p in D("peng_rfit2024_fig6")["panels"].items():
        out.append((f"P13-{n[:-3]}", "peng2024", float(n[:-3]), np.array(p["trace"]), 0.0))
    return out


def spec(trace, off, f0, at10=False, hi="slope"):
    return F.trace_spectrum(trace, off + (20 * np.log10(10.0 / f0) if at10 else 0.0), "A", hi)[0]


def idle_all(L):
    vr, ve = F.idle_var(L, TAU), F.echo_var(L, TAU)
    return dict(ramsey_var=vr, ramsey_avg=F.idle_infidelity(vr), ramsey_worst=F.idle_infidelity(vr, True),
                echo_var=ve, echo_avg=F.idle_infidelity(ve), echo_worst=F.idle_infidelity(ve, True))


def metrics(L):
    m = dict(brick=F.brick(L), gate={f"{fr/1e6:g}": F.gate(L, fr) for fr in FRS})
    m.update(idle_all(L))
    return m


def tail_share(L, f_last, fr=1e6, k=3.0):
    g = 2 * F.GRID**2 * L * h_lo2(F.GRID, fr) / fr**2
    m_all, m_tail = F.GRID <= f_last, (F.GRID >= k * fr) & (F.GRID <= f_last)
    return np.trapezoid((g * F.GRID)[m_tail], F.LNF[m_tail]) / np.trapezoid((g * F.GRID)[m_all], F.LNF[m_all])


def spearman(x, y):
    rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def discordant(x, y):
    n = len(x)
    return sum(1 for i in range(n) for j in range(i + 1, n) if (x[i] - x[j]) * (y[i] - y[j]) < 0), n * (n - 1) // 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(ROOT, "paper", "generated"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    N = {"conventions": {"tau_s": TAU, "target": TARGET, "band_Hz": [1.0, 1e8]}}

    # ---- Table III: own carrier and referred to 10 GHz --------------------------------------------------
    F.set_band(1.0, 1e8)
    P = points()
    rows = {}
    for key, cite, f0, tr, off in P:
        own, ten = spec(tr, off, f0), spec(tr, off, f0, at10=True)
        rows[key] = dict(cite=cite, f0_GHz=f0, f_last_Hz=float(tr[-1, 0]), own=metrics(own), at10=metrics(ten),
                         tail_share_3fR=tail_share(own, tr[-1, 0]))
    N["points"] = rows

    # ---- Section V-C: band sensitivity (flat floor at the last measured level) -------------------------
    sens = {}
    for fmax in (1e9, 1e10):
        F.set_band(1.0, fmax, n=800001)
        sens[f"{fmax:.0e}"] = {k: F.gate(spec(tr, off, f0, hi="flat"), 1e6) for k, _, f0, tr, off in P}
    F.set_band(1.0, 1e8)
    N["band_sensitivity_gate_fR1MHz"] = {
        fm: {k: dict(gate=v, ratio=v / rows[k]["own"]["gate"]["1"]) for k, v in d.items()} for fm, d in sens.items()}

    # ---- Section VI: correlation, gate average conversion ----------------------------------------------
    b = [rows[k]["own"]["brick"] for k in rows]
    g = [rows[k]["own"]["gate"]["1"] for k in rows]
    nd, npair = discordant(b, g)
    N["spearman_brick_gate1"] = spearman(b, g)
    N["discordant_pairs"] = [nd, npair]
    N["gate_avg_over_worst"] = 2.0 / 3.0
    for where in ("own", "at10"):
        N[f"pass_gate1_{where}"] = {"worst": [k for k in rows if rows[k][where]["gate"]["1"] <= TARGET],
                                    "average": [k for k in rows if 2 / 3 * rows[k][where]["gate"]["1"] <= TARGET]}
        for idle in ("ramsey_avg", "ramsey_worst", "echo_avg", "echo_worst"):
            N[f"fail_{idle}_{where}"] = [k for k in rows if rows[k][where][idle] > TARGET]
    N["echo_over_ramsey_locked"] = {k: rows[k]["own"]["echo_var"] / rows[k]["own"]["ramsey_var"] for k in rows}
    g300 = F.gong_spectrum(os.path.join(ROOT, "data", "gong_fig13a.json"), "A")[0]
    N["gong300K_gate1"] = F.gate(g300, 1e6)

    # ---- free-running oscillators -------------------------------------------------------------------
    from pll_comparison import fit_primary
    al, fc = fit_primary(np.array(D("patra_digitized")["fig22_PN_dBc_per_Hz_vs_offset_Hz"]["4K"]))
    qvar = (2 * np.pi * SIGMA_Q * TAU) ** 2
    fr = {}
    for name, (alpha, corner) in {"patra_6.3": (al, fc), "gong2022_4.65": (0.041, 1.3e6),
                                  "gong2022_6.3": (0.041 * (6.3 / 4.65) ** 2, 1.3e6)}.items():
        L = F.patra_free(alpha, corner)
        d = dict(alpha=alpha, fc=corner, brick=F.brick(L), **idle_all(L))
        d["sigma_eff_Hz"] = np.sqrt(d["ramsey_var"]) / (2 * np.pi * TAU)
        d["ramsey_var_over_qubit"] = d["ramsey_var"] / qvar
        d["echo_var_ratio"] = d["ramsey_var"] / d["echo_var"]
        d["echo_infid_ratio"] = d["ramsey_avg"] / d["echo_avg"]
        # frequency tracking that removes the qubit's quasi-static offset also removes the LO's own slow drift:
        # Ramsey with the LO noise below f_a = 8.3 kHz (about 1/T2*) removed
        vt = F.idle_var(np.where(F.GRID >= 8.3e3, L, 0.0), TAU)
        d["ramsey_tracked_var"], d["ramsey_tracked_avg"] = vt, F.idle_infidelity(vt)
        fr[name] = d
    N["free_running"] = fr
    N["qubit_own"] = dict(var=qvar, avg=F.idle_infidelity(qvar))
    N["sigma_eff_for_target_Hz"] = np.sqrt(6 * TARGET) / (2 * np.pi * TAU)
    N["fc_max_gong2022_6.3"] = {f"{fp:g}": B.fc_max(fr["gong2022_6.3"]["alpha"], fpll=fp) for fp in (3e5, 1e6)}

    # ---- Section VII: plateau jitter limits (locked LO, var ~ 2 sigma^2 Ramsey, 6 sigma^2 echo) ---------
    to_fs = lambda s, f0=10e9: s / (2 * np.pi * f0) * 1e15
    sr, se = np.sqrt(3 * TARGET), np.sqrt(TARGET)
    s60 = 2 * np.pi * 10e9 * 60e-15
    N["jitter_spec"] = dict(ramsey_mrad=sr * 1e3, ramsey_fs_10GHz=to_fs(sr), echo_mrad=se * 1e3, echo_fs_10GHz=to_fs(se),
                            gong60fs_sigma2=s60**2, gong60fs_ramsey=s60**2 / 3, gong60fs_echo=s60**2)

    # ---- Section VI: Patra's oscillator in the Gong-matched loop (type-II, zeta = 1.6), four input floors -------
    import pll
    arch, zeta = pll.GONG_MATCHED
    loop = {}
    for bw in (1e6, 4e6):
        h = pll.make_h_lp(arch, bw, zeta)(F.GRID)
        h2, hp2 = np.abs(h) ** 2, np.abs(1 - h) ** 2          # input low-pass, oscillator high-pass
        losc = F.patra_free(al, fc)
        floors = {"white_FM": al / bw**2, "flicker_aware": float(pll.l_osc(np.array([bw]), al, fc)[0]),
                  "flat_-122": 10 ** (-12.2), "flat_-131": 10 ** (-13.1)}
        loop[f"{bw:g}"] = {}
        for nm, lin in floors.items():
            L = h2 * lin + hp2 * losc                                       # Eq. (7)
            v = F.idle_var(L, TAU)
            loop[f"{bw:g}"][nm] = dict(brick=F.brick(L), gate1=F.gate(L, 1e6), ramsey_var=v, ramsey_avg=F.idle_infidelity(v))
    N["patra_in_loop"] = loop

    # ---- Section IV: where each weight looks, Gong [7] at 10 GHz --------------------------------------------
    g4 = rows["G7-10.0"]
    _, _, f0g, trg, offg = [t for t in P if t[0] == "G7-10.0"][0]
    L = spec(trg, offg, f0g)
    def share(y, lo, hi):
        m = (F.GRID >= lo) & (F.GRID <= hi)
        return float(np.trapezoid((y * F.GRID)[m], F.LNF[m]))
    sf = 2 * F.GRID**2 * L
    bt = share(sf, 8.3e3, 1e7)
    gw = sf * h_lo2(F.GRID, 1e6) / 1e12
    gt = F.integ(gw)
    N["weights_gong7"] = dict(brick_4_10MHz=share(sf, 4e6, 1e7) / bt, gate1_1_4MHz=share(gw, 1e6, 4e6) / gt,
                              gate1_4_10MHz=share(gw, 4e6, 1e7) / gt, fidelity_fR1=1 - g4["own"]["gate"]["1"],
                              fidelity_fR10=1 - g4["own"]["gate"]["10"])

    json.dump(N, open(os.path.join(a.out, "numbers.json"), "w"), indent=1, default=float)

    # ---- LaTeX rows ---------------------------------------------------------------------------------
    e6 = lambda v: (f"{v*1e6:.2f}" if v < 1e-6 else f"{v*1e6:.1f}")
    bold = lambda v, s: f"\\textbf{{{s}}}" if v > TARGET else s
    lines, last = [], None
    for key in rows:
        r = rows[key]
        o, t = r["own"], r["at10"]
        cite = f"\\cite{{{r['cite']}}}" if r["cite"] != last else ""
        last = r["cite"]
        cells = [cite, (f"{r['f0_GHz']:.1f}" if r["f0_GHz"] < 10 else f"{r['f0_GHz']:.1f}".rstrip("0").rstrip(".")), f"{o['brick']/1e3:.2f}"]
        cells += [bold(o["gate"][f"{fr_/1e6:g}"], e6(o["gate"][f"{fr_/1e6:g}"])) for fr_ in FRS]
        cells += [bold(o[k], e6(o[k])) for k in ("ramsey_avg", "echo_avg")]
        cells += [f"{t['brick']/1e3:.2f}", bold(t["gate"]["1"], e6(t["gate"]["1"]))]
        cells += [bold(t[k], e6(t[k])) for k in ("ramsey_avg", "echo_avg")]
        lines.append(" & ".join(cells) + "\\\\")
    # own-carrier-only variant (10 columns: ref, carrier, brick, gate x5, Ramsey, echo), as in the manuscript's Table II
    own_lines, last = [], None
    for key in rows:
        r, o = rows[key], rows[key]["own"]
        cite = f"\\cite{{{r['cite']}}}" if r["cite"] != last else ""
        last = r["cite"]
        f0s = f"{r['f0_GHz']:.1f}" if r["f0_GHz"] < 10 else f"{r['f0_GHz']:g}"
        cells = [cite, f0s, f"{o['brick']/1e3:.2f}"]
        cells += [bold(o["gate"][f"{fr_/1e6:g}"], e6(o["gate"][f"{fr_/1e6:g}"])) for fr_ in FRS]
        cells += [bold(o[k], e6(o[k])) for k in ("ramsey_avg", "echo_avg")]
        own_lines.append(" & ".join(cells) + "\\\\")
    open(os.path.join(a.out, "table_own_carrier_rows.tex"), "w").write(
        "% generated by scripts/table3.py -- own-carrier columns only\n" + "\n".join(own_lines) + "\n")
    open(os.path.join(a.out, "table3_rows.tex"), "w").write(
        "% generated by scripts/table3.py -- do not edit by hand\n" + "\n".join(lines) + "\n")
    print(open(os.path.join(a.out, "table3_rows.tex")).read())
    print(json.dumps({k: N[k] for k in N if k not in ("points", "band_sensitivity_gate_fR1MHz")}, indent=1, default=float))


if __name__ == "__main__":
    main()
