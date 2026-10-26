"""Figures for the TCAS-I manuscript (Figs. 1-7) and fig_numbers.json.

Run from the repo root:  python paper/make_figs.py [repo_root] [out_dir]
Defaults: repo_root = this repo, out_dir = paper/figs. Table III and the numbers quoted in the text: scripts/table3.py."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
R = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, "..")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(_HERE, "figs")
os.makedirs(OUT, exist_ok=True)
sys.path.insert(0, os.path.join(R, "src")); sys.path.insert(0, os.path.join(R, "scripts"))
import budget as B
import pll
import filter_criteria as F
from a0_gong_measured import h_lo2

plt.rcParams.update({"font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"], "mathtext.fontset": "stix",
                     "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "axes.linewidth": 0.6, "grid.linewidth": 0.4, "grid.color": "#d9d9d6", "lines.linewidth": 1.2,
                     "axes.edgecolor": "#52514e", "xtick.color": "#52514e", "ytick.color": "#52514e",
                     "axes.labelcolor": "#0b0b0b", "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.02})
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]          # validated categorical slots 1-4
INK2 = "#52514e"
D = lambda f: json.load(open(os.path.join(R, "data", f + ".json")))


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".pdf")); fig.savefig(os.path.join(OUT, name + ".png"), dpi=200); plt.close(fig)


# ---------- spectra ----------
S = {}
S["Gong [7], 10 GHz"], _ = F.gong_spectrum(os.path.join(R, "data", "gong_fig13b.json"), "A")
th = D("gong_thesis_fig721")["panels"]
for n, p in th.items():
    S[f"Gong [10], {n[:-3]} GHz"], _ = F.trace_spectrum(np.array(p["trace"]), 20 * np.log10(8), "A")
S["Xiao [12], 8.6 GHz"], _ = F.trace_spectrum(np.array(D("xiao_rfic2025_fig7")["trace"]), 0.0, "A")
for n, p in D("peng_rfit2024_fig6")["panels"].items():
    S[f"Peng [13], {n[:-3]} GHz"], _ = F.trace_spectrum(np.array(p["trace"]), 0.0, "A")
S300, _ = F.gong_spectrum(os.path.join(R, "data", "gong_fig13a.json"), "A")
A = np.array(D("patra_digitized")["fig22_PN_dBc_per_Hz_vs_offset_Hz"]["4K"])
from pll_comparison import fit_primary
al, fc = fit_primary(A)
PFR = F.patra_free(al, fc); PL1 = F.patra_in_loop(al, fc, 1e6); PL4 = F.patra_in_loop(al, fc, 4e6)
chip = lambda k: 0 if k.startswith("Gong [7]") else 1 if k.startswith("Gong [10]") else 2 if k.startswith("Xiao") else 3
MK = ["o", "s", "D", "^"]
rows = {k: (F.brick(L), F.gate(L, 1e6), F.gate(L, 4e6), (2 * np.pi * 1.2e-4) ** 2 * F.idle_sigma_eff(L, 1.2e-4) ** 2) for k, L in S.items()}
json.dump({k: list(map(float, v)) for k, v in rows.items()}, open(os.path.join(OUT, "fig_numbers.json"), "w"), indent=1)

# ---------- Fig. 1: corner tolerance ----------
fig, ax = plt.subplots(figsize=(3.45, 2.35))
a = np.linspace(0.02, 0.1842, 400)
for fp, col, lab in ((3e5, C[0], "300 kHz"), (1e6, C[1], "1 MHz")):
    for cont, ls in ((False, "-"), (True, "--")):
        y = np.array([B.fc_max(x, fpll=fp, continuous=cont) for x in a])
        m = y > 0
        ax.plot(a[m], y[m] / 1e6, ls, color=col, label=f"$f_{{PLL}}$ = {lab}, " + ("continuous" if cont else "fixed") + " floor")
ax.plot([al], [fc / 1e6], "k*", ms=8, mec="#fcfcfb", mew=0.5); ax.annotate("Patra 4 K oscillator [4]", (al, fc / 1e6), xytext=(0.128, 12.0), fontsize=7, arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2))
AG = 0.041 * (6.3 / 4.65) ** 2
ax.plot([AG], [1.3], "kP", ms=6, mec="#fcfcfb", mew=0.5); ax.annotate("Gong osc. [6] at 6.3 GHz", (AG, 1.3), xytext=(0.03, 0.68), fontsize=7, arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2))
ax.set_yscale("log"); ax.set_ylim(0.05, 40); ax.set_xlim(0.02, 0.2)
ax.set_xlabel(r"White-FM coefficient $\alpha$ (Hz$^2$/Hz)"); ax.set_ylabel(r"Tolerable corner $f_{c,\max}$ (MHz)")
ax.grid(True, which="major"); ax.legend(loc="lower left", handlelength=1.8, fontsize=6.5)
save(fig, "fig1_corner")

# ---------- Fig. 2: spectra and filters ----------
G = F.GRID
fig, (a1, a2) = plt.subplots(2, 1, figsize=(3.45, 4.1), sharex=True, gridspec_kw={"hspace": 0.12})
rep = ["Gong [7], 10 GHz", "Gong [10], 10.6 GHz", "Xiao [12], 8.6 GHz", "Peng [13], 6.0 GHz"]
for i, k in enumerate(rep):
    m = (G >= 1e4) & (G <= 3e7)
    a1.semilogx(G[m], 10 * np.log10(S[k][m]), color=C[i], label=k)
a1.set_ylabel("$\\mathcal{L}(f)$ at the PLL output (dBc/Hz)"); a1.set_ylim(-160, -95); a1.grid(True); a1.legend(loc="lower left", ncol=1)
a1.text(0.98, 0.94, "(a)", transform=a1.transAxes, ha="right", va="top")
f = np.logspace(3, 8, 4000)
a2.semilogx(f, np.where((f >= 8.3e3) & (f <= 1e7), 1.0, 0.0), color="#0b0b0b", lw=1.0, label="brick, 8.3 kHz–10 MHz [4]")
a2.semilogx(f, h_lo2(f, 1e6), color=C[1], label="gate, $\\pi$ rotation, $f_R$ = 1 MHz [7]")
a2.semilogx(f, np.sinc(f * 1e-5) ** 2, color=C[2], label="idle, $\\tau$ = 10 µs")
a2.set_ylabel("Weight on $S_f(f)$"); a2.set_xlabel("Offset frequency $f$ (Hz)"); a2.set_ylim(-0.03, 1.32); a2.grid(True)
a2.legend(loc="center left", bbox_to_anchor=(0.0, 0.42)); a2.set_xlim(1e3, 1e8)
a2.text(0.98, 0.94, "(b)", transform=a2.transAxes, ha="right", va="top")
save(fig, "fig2_spectra_filters")

# ---------- Fig. 3: brick vs gate ----------
fig, ax = plt.subplots(figsize=(3.45, 3.0))
ax.axvspan(1.0, 1.9, color="#eef3fb", lw=0); ax.axhspan(1e-8, 1e-5, color="#eef7f1", lw=0, zorder=0)
ax.axvline(1.9, color=INK2, lw=0.8, ls="--"); ax.axhline(1e-5, color=INK2, lw=0.8, ls="--")
labs = ["Gong [7]", "Gong [10]", "Xiao [12]", "Peng [13]"]
done = set()
for k, (b, g1, g4, iv) in rows.items():
    i = chip(k)
    ax.plot(b / 1e3, g1, MK[i], color=C[i], ms=6, mec="#fcfcfb", mew=0.6, label=labs[i] if i not in done else None); done.add(i)
pts = [("Gong [7], 300 K", F.brick(S300), F.gate(S300, 1e6), "x"), ("Patra osc., free running", F.brick(PFR), F.gate(PFR, 1e6), "v"),
       ("Patra osc., 1-MHz loop", F.brick(PL1), F.gate(PL1, 1e6), "<"), ("Patra osc., 4-MHz loop", F.brick(PL4), F.gate(PL4, 1e6), ">")]
for lab, b, g, mk in pts:
    ax.plot(b / 1e3, g, mk, color=INK2, ms=5, mfc="none", mew=0.8, label=lab)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1.0, 25); ax.set_ylim(2e-7, 4e-5)
ax.set_xlabel("Brick-wall FN, 8.3 kHz–10 MHz (kHz rms)"); ax.set_ylabel("Gate infidelity $1-F$, $f_R$ = 1 MHz")
ax.text(1.05, 2.4e-7, "passes 1.9 kHz", fontsize=6.5, color=INK2); ax.text(24, 1.1e-5, "fails 99.999%", fontsize=6.5, color=INK2, ha="right")
ax.text(24, 8.4e-6, "meets 99.999%", fontsize=6.5, color=INK2, ha="right", va="top")
ax.set_xticks([1, 2, 5, 10, 20]); ax.set_xticklabels(["1", "2", "5", "10", "20"])
ax.grid(True, which="major"); ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), fontsize=6.3, ncol=2, columnspacing=0.8, handletextpad=0.3)
save(fig, "fig3_brick_vs_gate")

# ---------- Fig. 4: idle phase variance vs tau ----------
fig, ax = plt.subplots(figsize=(3.45, 2.35))
taus = np.logspace(-7, -3, 41)
v = lambda L: np.array([(2 * np.pi * t) ** 2 * F.idle_sigma_eff(L, t) ** 2 for t in taus])
ax.loglog(taus * 1e6, (2 * np.pi * 1.9e3 * taus) ** 2, color="#0b0b0b", lw=1.0, ls="--", label="qubit's own (1.9 kHz, $T_2^*$ = 118 µs)")
ax.loglog(taus * 1e6, v(PFR), color=INK2, lw=1.2, label="Patra osc., free running [4]")
for i, k in enumerate(rep):
    ax.loglog(taus * 1e6, v(S[k]), color=C[i], label=k)
ax.set_xlabel("Idle duration $\\tau$ (µs)"); ax.set_ylabel("Phase variance (rad$^2$)"); ax.set_ylim(1e-7, 1e2); ax.grid(True, which="major")
ax.legend(loc="upper left", fontsize=6.3)
save(fig, "fig4_idle")
print("figures written;", len(rows), "measured 4 K points")

# ---------- Fig. 5: FN versus loop bandwidth (Patra 4-K oscillator, white-FM floor) ----------
bw = np.logspace(5, np.log10(9.9e6), 61)
fig, ax = plt.subplots(figsize=(3.45, 2.5))
ax.semilogx(bw / 1e6, [B.fn(al, fc, fpll=b) / 1e3 for b in bw], color="#0b0b0b", lw=1.0, label="brick wall, fixed floor")
ax.semilogx(bw / 1e6, [B.fn(al, fc, fpll=b, continuous=True) / 1e3 for b in bw], color="#0b0b0b", lw=1.0, ls="--", label="brick wall, continuous floor")
fams = [("first-order", None, C[0], "first order"), ("type-II (PI zero)", 1.0, C[1], "type II, $\\zeta$ = 1"),
        ("2nd-order low-pass", 1.6, C[2], "second-order low-pass, $\\zeta$ = 1.6"), ("type-II (PI zero)", 1.6, C[3], "type II, $\\zeta$ = 1.6 (as in [7])")]
for arch, z, col, lab in fams:
    y = [pll.fn_rms(pll.make_h_lp(arch, b, z), al, fc, al / b ** 2) / 1e3 for b in bw]
    ax.semilogx(bw / 1e6, y, color=col, label=lab)
ax.axhline(1.9, color=INK2, lw=0.8, ls=":"); ax.text(0.105, 1.93, "budget 1.9 kHz", fontsize=6.5, color=INK2, va="bottom")
ax.axhline(np.sqrt(B.fn2_free_running(al, fc)) / 1e3, color=INK2, lw=0.8, ls="-."); ax.text(0.105, 3.31, "free running", fontsize=6.5, color=INK2, va="bottom")
ax.set_xlim(0.1, 10); ax.set_ylim(1.0, 4.0); ax.set_xlabel("Loop bandwidth (MHz)"); ax.set_ylabel("FN, 8.3 kHz–10 MHz (kHz rms)")
ax.set_xticks([0.1, 0.3, 1, 3, 10]); ax.set_xticklabels(["0.1", "0.3", "1", "3", "10"])
ax.grid(True, which="major"); ax.legend(loc="upper right", fontsize=6.0)
save(fig, "fig5_pll_bandwidth")

# ---------- Fig. 6: fit identifiability over 16 windows ----------
fits = [(lo, hi, *fit_primary(A, lo, hi)) for lo in (1.5e5, 2e5, 3e5, 5e5) for hi in (5e6, 1e7, 2e7, 2.9e7)]
fig, ax = plt.subplots(figsize=(3.45, 2.3))
Ks = [a * f for *_, a, f in fits]
xx = np.linspace(0.05, 0.3, 100)
for Kv, ls in ((min(Ks), "--"), (max(Ks), "--")):
    ax.plot(xx, Kv / xx / 1e6, color=INK2, lw=0.8, ls=ls)
ax.text(0.205, max(Ks) / 0.205 / 1e6 * 1.12, f"$K$ = {max(Ks)/1e5:.1f}$\\times10^5$ Hz$^2$", fontsize=6.5, color=INK2)
ax.text(0.19, 1.62, f"$K$ = {min(Ks)/1e5:.1f}$\\times10^5$ Hz$^2$", fontsize=6.5, color=INK2)
los = sorted({f[0] for f in fits})
for i, lo in enumerate(los):
    pts = [(a, f) for l, h, a, f in fits if l == lo]
    ax.plot([p[0] for p in pts], [p[1] / 1e6 for p in pts], MK[i], color=C[i], ms=5, mec="#fcfcfb", mew=0.5, label=f"lower edge {lo/1e3:.0f} kHz")
ax.plot([al], [fc / 1e6], "k*", ms=8, mec="#fcfcfb", mew=0.5, label="primary fit (200 kHz–10 MHz)")
ax.axhline(5.7, color=INK2, lw=0.6, ls=":"); ax.text(0.225, 6.0, "[4, Fig. 24]", fontsize=6.5, color=INK2)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.06, 0.3); ax.set_ylim(1.5, 15)
ax.set_xticks([0.07, 0.1, 0.15, 0.2, 0.3]); ax.set_xticklabels(["0.07", "0.1", "0.15", "0.2", "0.3"])
ax.set_yticks([2, 3, 5, 10]); ax.set_yticklabels(["2", "3", "5", "10"])
from matplotlib.ticker import NullFormatter
ax.xaxis.set_minor_formatter(NullFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())
ax.set_xlabel(r"Fitted $\alpha$ (Hz$^2$/Hz)"); ax.set_ylabel(r"Fitted $f_c$ (MHz)"); ax.grid(True, which="major")
ax.legend(loc="lower left", fontsize=6.0)
save(fig, "fig6_identifiability")
print("fig5, fig6 written")


# ---------- Fig. 7: gate infidelity versus rotation speed (native carriers) ----------
frs = np.logspace(np.log10(3e5), np.log10(2e7), 31)
fig, ax = plt.subplots(figsize=(3.45, 2.5))
for thr, lab in ((1e-4, "$10^{-4}$"), (1e-5, "$10^{-5}$ [7]"), (1e-6, "$10^{-6}$")):
    ax.axhline(thr, color=INK2, lw=0.7, ls=":"); ax.text(17, thr * 1.12, lab, fontsize=6.5, color=INK2, ha="right")
done = set()
for k, L in S.items():
    i = chip(k)
    ax.loglog(frs / 1e6, [F.gate(L, x) for x in frs], color=C[i], lw=1.0, ls=["-", "-", "-", "--"][i], label=labs[i] if i not in done else None); done.add(i)
ax.set_xlim(0.3, 20); ax.set_ylim(5e-8, 3e-4)
ax.set_xticks([0.3, 1, 3, 10]); ax.set_xticklabels(["0.3", "1", "3", "10"])
ax.set_xlabel("Rotation speed $f_R$ (MHz)"); ax.set_ylabel("Gate infidelity $1-F$"); ax.grid(True, which="major")
ax.legend(loc="lower left", fontsize=6.3)
save(fig, "fig7_gate_vs_fr")
print("fig7 written")
