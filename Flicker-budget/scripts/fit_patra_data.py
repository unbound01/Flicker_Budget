"""Fit (alpha, fc) to the digitized Patra Figs 22/24/25; test the 3.4 kHz and '<= 1 MHz PLL' claims.
Run from the repo root:  python scripts/fit_patra_data.py   (writes figures/patra_fit.png, which is git-ignored)"""
import json, sys, os
import numpy as np
from scipy.optimize import least_squares, brentq
here=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(here,'..','src'))
import budget as c8
DATA_PATH=os.path.join(here,'..','data','patra_digitized.json')
if not os.path.exists(DATA_PATH):
    sys.exit("data/patra_digitized.json not found. It is not distributed with this repo (derived from figures of a\n"
             "copyrighted paper). Generate it from your own copy of the paper:\n"
             "  python scripts/digitize_patra_figs.py --pdf PATH/TO/CryoCMOS_BPatra.pdf --out data/patra_digitized.json")
D=json.load(open(DATA_PATH))
def Ldb(f,al,fc): return 10*np.log10(al/f**2*(1+fc/f))
def fit(f,y,p0=(0.3,1e6)):
    r=least_squares(lambda p:Ldb(f,*np.exp(p))-y,np.log(p0)); return (*np.exp(r.x),np.sqrt(np.mean(r.fun**2)))
def curve(T,lo,hi):
    a=np.array(D['fig22_PN_dBc_per_Hz_vs_offset_Hz'][T]); m=(a[:,0]>=lo)&(a[:,0]<=hi*1.001); return a[m,0],a[m,1]
FNp=lambda al,fc,fp,cont:c8.fn(al,fc,fpll=fp,continuous=cont)
print("== Fig 22 fits (range 200k-10M) ==")
fits={}
for T in ('4K','165K','300K'):
    al,fc,rms=fit(*curve(T,2e5,1e7)); fits[T]=(al,fc)
    print(f"{T:5s} alpha={al:.3f} fc={fc/1e6:.3f} MHz K={al*fc:.3g} Hz^2 rms={rms:.2f} dB  free-run FN={np.sqrt(c8.fn2_free_running(al,fc))/1e3:.2f} kHz")
print("== 4 K sensitivity (fit range, +-0.5 dB y-shift) ==")
for lo,hi,sh in ((2e5,1e7,0),(1.5e5,2.9e7,0),(3e5,1e7,0),(2e5,1e7,.5),(2e5,1e7,-.5)):
    f,y=curve('4K',lo,hi); al,fc,rms=fit(f,y+sh)
    print(f" {lo/1e3:4.0f}k-{hi/1e6:4.1f}M {sh:+.1f}dB: a={al:.3f} fc={fc/1e6:.2f}M free={np.sqrt(c8.fn2_free_running(al,fc))/1e3:.2f} "
          f"FN@300k={FNp(al,fc,3e5,False)/1e3:.2f}/{FNp(al,fc,3e5,True)/1e3:.2f} FN@1M={FNp(al,fc,1e6,False)/1e3:.2f}/{FNp(al,fc,1e6,True)/1e3:.2f} kHz (fixed/cont)")
al,fc=fits['4K']
print(f"fc_max at alpha={al:.3f}: fPLL=300k {c8.fc_max(al,fpll=3e5)/1e6:.2f} MHz, fPLL=1M {c8.fc_max(al,fpll=1e6)/1e6:.2f} MHz (fitted fc={fc/1e6:.2f} MHz)")
for cont in (False,True):
    fp=brentq(lambda x:FNp(al,fc,x,cont)-1.9e3,1e5,9e6); print(f"PLL bandwidth for 1.9 kHz ({'cont' if cont else 'fixed'}): {fp/1e6:.2f} MHz")
print("== Fig 25 exact 2-point solve (1 & 10 MHz) ==")
P=D['fig25_PN_dBc_per_Hz_vs_GHz']; tab=[]
for i in range(8):
    g=P['1MHz'][i]['f_GHz']; p1=P['1MHz'][i]['PN']; p10=P['10MHz'][i]['PN']
    r1,r10=10**(p1/10)*1e12,10**(p10/10)*1e14; K=(r1-r10)/(9e-7); a_=r1-K*1e-6; f_=K/a_
    tab.append((g,a_,f_,np.sqrt(c8.fn2_free_running(a_,f_)),FNp(a_,f_,1e6,False)))
    print(f" {g:.2f} GHz alpha={a_:.3f} fc={f_/1e6:5.1f} MHz FNfree={tab[-1][3]/1e3:.2f} FN@1M={tab[-1][4]/1e3:.2f} kHz")
print("== K = alpha*fc vs temperature (Fig 22 fits; Fig 24 corner x fitted alpha) ==")
c24={round(t):c for t,c in [(r['T_K'],r['corner_Hz']) for r in D['fig24_corner_Hz_vs_T_K']]}
for T,Tk in (('4K',4),('300K',300)):
    a_,f_=fits[T]; print(f" {T}: fit K={a_*f_:.3g}, with Fig24 corner K={a_*[v for k,v in c24.items() if abs(k-Tk)<2][0]:.3g} Hz^2")
print("unit test (closed form vs quad):",c8.unit_test())
# figure
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,3,figsize=(15,4.3))
for T,c in (('4K','r'),('165K','g'),('300K','b')):
    f,y=curve(T,1.5e5,2.9e7); ax[0].semilogx(f,y,c+'.',ms=4); a_,f_=fits[T]; ff=np.logspace(5,7.5,100); ax[0].semilogx(ff,Ldb(ff,a_,f_),c+'-',lw=1,label=f"{T}: a={a_:.2f}, fc={f_/1e6:.2f} MHz")
ax[0].set_xlabel('offset [Hz]'); ax[0].set_ylabel('PN [dBc/Hz]'); ax[0].legend(fontsize=8); ax[0].set_title('(a) Fig 22 digitized + fit')
T24=np.array([[r['T_K'],r['corner_Hz']] for r in D['fig24_corner_Hz_vs_T_K']])
ax[1].loglog(T24[:,0],T24[:,1]/1e6,'ko-',label='Fig 24 (digitized)')
for T,Tk,c in (('4K',4,'r'),('165K',165,'g'),('300K',300,'b')): ax[1].loglog(Tk,fits[T][1]/1e6,c+'s',label=f'fit {T}')
ax[1].set_xlabel('T [K]'); ax[1].set_ylabel('1/f^3 corner [MHz]'); ax[1].legend(fontsize=8); ax[1].set_title('(b) corner: Fig 24 vs fit')
fp=np.logspace(5,6.9,100)
for cont,ls in ((False,'-'),(True,':')): ax[2].semilogx(fp/1e6,[FNp(al,fc,x,cont)/1e3 for x in fp],'r'+ls,label=f"4 K fit, {'continuous' if cont else 'fixed'} floor")
ax[2].axhline(1.9,color='k',ls='--',label='1.9 kHz'); ax[2].axvline(1,color='gray',ls=':'); ax[2].axhline(np.sqrt(c8.fn2_free_running(al,fc))/1e3,color='gray',ls='-.',label='free-running')
ax[2].set_xlabel('f_PLL [MHz]'); ax[2].set_ylabel('FN [kHz]'); ax[2].legend(fontsize=8); ax[2].set_title('(c) FN vs PLL bandwidth')
fig.tight_layout(); os.makedirs(os.path.join(here,'..','figures'),exist_ok=True); fig.savefig(os.path.join(here,'..','figures','patra_fit.png'),dpi=140); print('saved figure')
