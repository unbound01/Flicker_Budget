"""Linear phase-domain PLL noise model (Phase 3).

Output phase noise (single-sideband L, 1/Hz, same convention as budget.py):

    L_out(f) = |H_lp(f)|^2 * L_in(f)  +  |1 - H_lp(f)|^2 * L_osc(f)

  L_osc : free-running oscillator, alpha/f^2 * (1 + fc/f)
  L_in  : everything entering at the phase detector, referred to the OUTPUT carrier
          (reference noise * N^2, detector, divider, input buffer); modelled as a flat in-band floor
  H_lp  : closed-loop transfer from the (output-referred) input to the output; the oscillator sees 1 - H_lp.

Because L_in is output-referred, the division ratio N is already inside L_in. Conventions are
single-sideband throughout; FN^2 = 2 * int f^2 L_out df, exactly as in budget.py / Patra Eq. 8.

"Bandwidth" for the generic loops is the closed-loop -3 dB frequency of H_lp. For the type-II family this is the
same definition Gong uses: BW = sqrt(2z^2 + 1 + sqrt((2z^2 + 1)^2 + 1)) * wn / (2 pi)  (thesis, Ch. 5 footnote).
"""
import numpy as np
from scipy.optimize import brentq

from budget import FA, FB

_trapz = getattr(np, "trapezoid", None) or np.trapz      # NumPy < 2.0 compatibility


def l_osc(f, alpha, fc):
    return alpha / f**2 * (1.0 + fc / f)


# ---- closed-loop low-pass transfer functions in units of a corner fn (complex) --------------------
def _first(r, zeta):
    return 1.0 / (1.0 + 1j * r)


def _type2(r, zeta):                     # charge-pump PLL with PI filter: H = (2 z s + 1)/(s^2 + 2 z s + 1)
    s = 1j * r
    return (2 * zeta * s + 1) / (s**2 + 2 * zeta * s + 1)


def _lp2(r, zeta):                       # type-I loop with an extra pole: H = 1/(s^2 + 2 z s + 1)
    s = 1j * r
    return 1.0 / (s**2 + 2 * zeta * s + 1)


def _butter4(r, zeta):                   # steep stable loop: closed-loop 4th-order Butterworth (extreme case)
    s = 1j * r
    return 1.0 / (s**4 + 2.613126 * s**3 + 3.414214 * s**2 + 2.613126 * s + 1)


ARCH = {
    "first-order": (_first, 1.0),        # (function, default damping)
    "type-II (PI zero)": (_type2, 1.0),
    "2nd-order low-pass": (_lp2, 1.6),   # generic family; NOT Gong's loop (see GONG_MATCHED below)
    "4th-order Butterworth (steep)": (_butter4, 1.0),   # damping unused; shows the cost of a steep roll-off
}


# Gong et al. JSSC 2023 DAPLL: V/I stage driving (R_Z + 1/sC_Z) || 1/sC_P (thesis Fig. 6.6(a)) -> type-II loop with
# a zero; zeta ~ 1.6 and BW = 4 MHz (paper pp. 1363, 1369). The extra C_P pole and the zero-order hold of K_PD(s)
# are NOT included here; this is the closest generic family, not a model of the circuit.
GONG_MATCHED = ("type-II (PI zero)", 1.6)


def ratio_3db(arch, zeta):
    """f_3dB / fn for the given architecture and damping (|H_lp|^2 = 1/2)."""
    fun = ARCH[arch][0]
    g = lambda r: abs(fun(r, zeta)) ** 2 - 0.5
    return brentq(g, 1e-3, 1e3) if g(1e-3) * g(1e3) < 0 else np.nan


def make_h_lp(arch, f3, zeta=None):
    """Return H_lp(f) for a loop whose closed-loop -3 dB bandwidth is f3."""
    fun, z0 = ARCH[arch]
    z = z0 if zeta is None else zeta
    fn = f3 / ratio_3db(arch, z)
    return lambda f: fun(np.asarray(f) / fn, z)


def brick_wall_h_lp(fp):
    return lambda f: (np.asarray(f) <= fp).astype(complex)


# ---- integrated frequency noise --------------------------------------------------------------------
def _grid(fa, fb, n):
    return np.logspace(np.log10(fa), np.log10(fb), n)


def l_osc_from_data(f_data, y_dbc, alpha, fc):
    """Oscillator spectrum taken directly from digitized points: log-log interpolation in dBc/Hz above the first
    data point, the fitted two-term model below it. Returns a callable L(f) in 1/Hz."""
    f_data, y_dbc = np.asarray(f_data, float), np.asarray(y_dbc, float)

    def fun(f):
        f = np.asarray(f, float)
        out = np.empty_like(f)
        m = f >= f_data[0]
        out[m] = 10 ** (np.interp(np.log(f[m]), np.log(f_data), y_dbc) / 10.0)
        out[~m] = l_osc(f[~m], alpha, fc)
        return out
    return fun


def fn_components(h_lp, alpha, fc, l_in, fa=FA, fb=FB, n=8001, osc=None):
    """Return (osc_part, in_part) of FN^2 in Hz^2. l_in: float (flat, 1/Hz) or callable of f.
    osc: optional callable L_osc(f) (e.g. from l_osc_from_data); default is the two-term model alpha/f^2 (1 + fc/f)."""
    f = _grid(fa, fb, n)
    H = h_lp(f)
    lin = l_in(f) if callable(l_in) else l_in
    l_o = osc(f) if osc is not None else l_osc(f, alpha, fc)
    o_term = np.abs(1.0 - H) ** 2 * l_o
    inb = np.abs(H) ** 2 * lin
    w = 2.0 * f**3                                   # 2 f^2 * df, df = f dln f
    lnf = np.log(f)
    return _trapz(w * o_term, lnf), _trapz(w * inb, lnf)


def fn_rms(h_lp, alpha, fc, l_in, **kw):
    o, i = fn_components(h_lp, alpha, fc, l_in, **kw)
    return np.sqrt(o + i)


def dbc_to_lin(x):
    return 10 ** (np.asarray(x) / 10.0)
