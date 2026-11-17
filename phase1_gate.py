"""
Phase 1 verification gate: single-cell current clamp for every cell type.

Checks, per cell: resting Vm, input resistance, Ih sag, f-I curve, AP shape
(peak, half-width, AHP), spike-frequency adaptation; thalamic burst vs tonic
mode and rebound bursts; PT BAC firing (somatic AP + distal tuft input).

Usage:  python phase1_gate.py        (needs compiled mod/ in the cwd)
Writes: results/phase1_gate.json (all values + checks), results/phase1_gate.png (PT BAC)
"""
import json, os, sys
import numpy as np
from neuron import h
from netpyne import specs, sim
from cells import all_cell_params

h.load_file('stdrun.hoc')
os.makedirs('results', exist_ok=True)
DT = 0.025

# ------------------------------------------------------------ build cells
netParams = specs.NetParams()
cfg = specs.SimConfig()
cfg.duration, cfg.dt, cfg.hParams = 100, DT, {'celsius': 36, 'v_init': -70}
cfg.verbose, cfg.printRunTime = False, False
cfg.recordStep = DT
for name, cp in all_cell_params().items():
    netParams.cellParams[name] = cp
    netParams.popParams[name] = {'cellType': name, 'numCells': 1}
sim.initialize(netParams=netParams, simConfig=cfg)
sim.net.createPops(); sim.net.createCells()
h.celsius = 36
CELLS = {c.tags['pop']: c for c in sim.net.cells}


def sec(cell, name='soma'):
    return CELLS[cell].secs[name]['hObj']


_ACTIVE, _PENDING = [], []


def run(tstop, v_init=-70.0, settle=True):
    # stop vectors from earlier protocols so they don't get overwritten
    global _ACTIVE, _PENDING
    for vec in _ACTIVE:
        vec.play_remove()
    _ACTIVE, _PENDING = _PENDING, []
    h.dt = DT
    h.finitialize(v_init)
    if settle:  # let channels settle before the protocol starts
        h.t = -1000.0; h.dt = 1.0
        while h.t < -1.0: h.fadvance()
        h.dt = DT; h.t = 0.0
        h.fcurrent()
    h.frecord_init()
    h.continuerun(tstop)


def record(s, x=0.5):
    vt, tt = h.Vector(), h.Vector()
    vt.record(s(x)._ref_v); tt.record(h._ref_t)
    _PENDING.extend([vt, tt])
    return tt, vt


def iclamp(s, amp, delay, dur, x=0.5):
    ic = h.IClamp(s(x)); ic.amp, ic.delay, ic.dur = amp, delay, dur
    return ic


def spikes(t, v, thr=-10.0):
    t, v = np.array(t), np.array(v)
    idx = np.where((v[:-1] < thr) & (v[1:] >= thr))[0]
    return t[idx]


def ap_shape(t, v, tsp):
    """Peak, half-width, AHP depth of the first spike."""
    t, v = np.array(t), np.array(v)
    i0 = np.searchsorted(t, tsp - 3); i1 = np.searchsorted(t, tsp + 15)
    seg_t, seg_v = t[i0:i1], v[i0:i1]
    dv = np.gradient(seg_v, seg_t)
    ith = np.where(dv > 20)[0]
    vth = seg_v[ith[0]] if len(ith) else seg_v[0]
    ipk = np.argmax(seg_v); vpk = seg_v[ipk]
    half = (vpk + vth) / 2
    # contiguous region above half-amplitude around the peak only
    lo = ipk
    while lo > 0 and seg_v[lo - 1] > half: lo -= 1
    hi = ipk
    while hi < len(seg_v) - 1 and seg_v[hi + 1] > half: hi += 1
    hw = seg_t[hi] - seg_t[lo]
    ahp = vth - seg_v[ipk:].min()
    return dict(threshold=float(vth), peak=float(vpk), halfwidth=float(hw), ahp=float(ahp))


def passive(cell):
    s = sec(cell)
    tt, vt = record(s)
    ic = iclamp(s, -0.05, 200, 400)
    run(700)
    t, v = np.array(tt), np.array(vt)
    vrest = float(v[(t > 150) & (t < 200)].mean())
    vss = float(v[(t > 550) & (t < 600)].mean())
    vmin = float(v[(t > 200) & (t < 600)].min())
    rin = (vrest - vss) / 0.05  # MOhm
    sag = (vss - vmin) / (vrest - vmin) if vrest != vmin else 0.0
    return dict(vrest=vrest, rin_MOhm=float(rin), sag_ratio=float(sag)), (t, v)


def fi_curve(cell, amps):
    s = sec(cell); out, traces = [], {}
    for a in amps:
        tt, vt = record(s)
        ic = iclamp(s, a, 100, 500)
        run(700)
        t, v = np.array(tt), np.array(vt)
        sp = spikes(t, v); sp = sp[(sp > 100) & (sp < 600)]
        out.append(len(sp) / 0.5)
        traces[a] = (t, v, sp)
    return np.array(out), traces


def holding_current(cell, target, lo=-1.0, hi=1.0):
    """Bisection for the DC current that holds soma at `target` mV (no spikes)."""
    s = sec(cell)
    for _ in range(18):
        mid = (lo + hi) / 2
        tt, vt = record(s); ic = iclamp(s, mid, -2000, 1e9)  # on during settle
        run(800)
        v = np.array(vt)[-4000:].mean()
        lo, hi = (mid, hi) if v < target else (lo, mid)
    return (lo + hi) / 2


def thalamic_modes(cell, step):
    """Burst mode from -75 mV, tonic mode from -58 mV, rebound after hyperpolarization."""
    s = sec(cell); res = {}; traces = {}
    for mode, vh in (('burst', -75.0), ('tonic', -58.0)):
        ih = holding_current(cell, vh)
        tt, vt = record(s)
        a = iclamp(s, ih, -2000, 1e9); b = iclamp(s, step, 300, 200)  # hold on during settle
        run(600)
        t, v = np.array(tt), np.array(vt)
        sp = spikes(t, v); sp = sp[sp > 300]
        isi = np.diff(sp)
        res[mode] = dict(hold_nA=float(ih), n_spikes=int(len(sp)),
                         first_isi_ms=float(isi[0]) if len(isi) else None,
                         latency_ms=float(sp[0] - 300) if len(sp) else None)
        traces[mode] = (t, v)
    tt, vt = record(s)
    a = iclamp(s, -0.3 if cell == 'TC' else -0.2, 200, 300)
    run(800)
    t, v = np.array(tt), np.array(vt)
    sp = spikes(t, v); sp = sp[(sp > 500) & (sp < 700)]
    res['rebound'] = dict(n_spikes=int(len(sp)),
                          first_isi_ms=float(np.diff(sp)[0]) if len(sp) > 1 else None)
    traces['rebound'] = (t, v)
    return res, traces


def bac_test():
    """Larkum 1999 protocol: somatic 5 ms pulse; EPSC-shaped tuft current; both (5 ms lag)."""
    soma, tuft = sec('PT'), sec('PT', 'apic_tuft')
    res, traces = {}, {}
    for cond in ('soma', 'tuft', 'both'):
        ts, vs = record(soma); _, vd = record(tuft, 0.3)
        objs = []
        if cond in ('soma', 'both'):
            objs.append(iclamp(soma, 1.0, 100, 5))
        if cond in ('tuft', 'both'):
            # double-exponential EPSC (tau 0.5/5 ms) via a Vector-played IClamp
            ic = h.IClamp(tuft(0.3)); ic.delay, ic.dur = 0, 1e9
            tvec = np.arange(0, 400, DT)
            onset = 105.0
            x = np.clip(tvec - onset, 0, None)
            wave = np.where(tvec >= onset, np.exp(-x / 5.0) - np.exp(-x / 0.5), 0.0)
            wave = 0.8 * wave / wave.max()   # 0.8 nA peak
            iv, tv = h.Vector(wave), h.Vector(tvec)
            iv.play(ic._ref_amp, tv, True)
            objs += [ic, iv, tv]
        run(300)
        t, v, vdd = np.array(ts), np.array(vs), np.array(vd)
        sp = spikes(t, v); sp = sp[sp > 95]
        ca_ms = float(((vdd > -20) & (t > 95)).sum() * DT)
        res[cond] = dict(n_soma_spikes=int(len(sp)), tuft_peak_mV=float(vdd[t > 95].max()),
                         tuft_ms_above_m20=ca_ms)
        traces[cond] = (t, v, vdd)
    return res, traces


# --------------------------------------------------------------- run gate
AMPS = [0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6]
RANGES = {  # physiological target windows used by the gate
    'vrest':     {'IT': (-75, -62), 'CT': (-75, -62), 'PT': (-72, -60), 'PV': (-75, -62),
                  'SOM': (-75, -58), 'VIP': (-72, -58), 'TC': (-72, -58), 'TRN': (-85, -65)},
    'rin_MOhm':  {'IT': (60, 300), 'CT': (60, 300), 'PT': (25, 150), 'PV': (60, 250),
                  'SOM': (100, 450), 'VIP': (150, 700), 'TC': (50, 300), 'TRN': (50, 400)},
    'halfwidth': {'IT': (0.5, 2.0), 'CT': (0.5, 2.0), 'PT': (0.5, 2.0), 'PV': (0.15, 0.8),
                  'SOM': (0.4, 2.0), 'VIP': (0.4, 2.0), 'TC': (0.3, 2.0), 'TRN': (0.2, 1.5)},
    'adaptation':{'IT': (1.3, 20), 'CT': (1.2, 20), 'PT': (1.0, 20), 'PV': (0.7, 1.3),
                  'SOM': (1.5, 50), 'VIP': (1.2, 30)},
}

report, checks = {}, []


def check(cell, key, val):
    lo, hi = RANGES[key][cell]
    ok = val is not None and not np.isnan(val) and lo <= val <= hi
    checks.append((cell, key, val, (lo, hi), ok))
    return ok


for cell in CELLS:
    r = {}
    r['passive'], _ = passive(cell)
    rates, ftr = fi_curve(cell, AMPS)
    r['fI_Hz'] = dict(zip(map(str, AMPS), rates.tolist()))
    spiking = [a for a, f in zip(AMPS, rates) if f > 0]
    r['rheobase_nA'] = spiking[0] if spiking else None
    # AP shape at rheobase; adaptation at the first step with >= 4 spikes
    if spiking:
        t, v, sp = ftr[spiking[0]]
        r['ap'] = ap_shape(t, v, sp[0])
    adapt = None
    for a in AMPS:
        sp = ftr[a][2]
        if len(sp) >= 4:
            isi = np.diff(sp); adapt = float(isi[-1] / isi[0]); break
    r['adaptation_ratio'] = adapt
    check(cell, 'vrest', r['passive']['vrest'])
    check(cell, 'rin_MOhm', r['passive']['rin_MOhm'])
    if 'ap' in r:
        check(cell, 'halfwidth', r['ap']['halfwidth'])
        checks.append((cell, 'AP peak > 0 mV', r['ap']['peak'], (0, 60), r['ap']['peak'] > 0))
    if cell in RANGES['adaptation']:
        check(cell, 'adaptation', adapt)
    if cell in ('IT', 'CT', 'PT'):  # no depolarization block: still firing at end of a 0.4 nA step
        sp = ftr[0.4][2]; last = float(sp[-1]) if len(sp) else 0.0
        checks.append((cell, 'sustained firing @0.4 nA (last spike >550 ms)', last, (550, 600), last > 550))
    if cell in ('TC', 'TRN'):
        m, _ = thalamic_modes(cell, step=0.2 if cell == 'TC' else 0.08)
        r['modes'] = m
        b, tn, rb = m['burst'], m['tonic'], m['rebound']
        checks.append((cell, 'burst: >=2 spikes at ISI<8ms from -75',
                       b['first_isi_ms'], (0, 8),
                       b['n_spikes'] >= 2 and b['first_isi_ms'] is not None and b['first_isi_ms'] < 8))
        checks.append((cell, 'tonic: regular firing (1st ISI>10ms) from -58',
                       tn['first_isi_ms'], (10, 1e3),
                       tn['n_spikes'] >= 2 and tn['first_isi_ms'] is not None and tn['first_isi_ms'] > 10))
        checks.append((cell, 'rebound burst after hyperpolarization',
                       rb['n_spikes'], (2, 50), rb['n_spikes'] >= 2))
    if cell == 'SOM':
        s = sec('SOM'); tt, vt = record(s); ic = iclamp(s, -0.1, 200, 300); run(800)
        sp = spikes(np.array(tt), np.array(vt)); n = int(((sp > 500) & (sp < 650)).sum())
        r['rebound_spikes'] = n
        checks.append(('SOM', 'LTS rebound spike', n, (1, 50), n >= 1))
    report[cell] = r

bac, bactr = bac_test()
report['PT']['BAC'] = bac
checks.append(('PT', 'BAC: soma pulse alone -> 1 spike', bac['soma']['n_soma_spikes'], (1, 1),
               bac['soma']['n_soma_spikes'] == 1))
checks.append(('PT', 'BAC: tuft input alone -> no somatic spike', bac['tuft']['n_soma_spikes'], (0, 0),
               bac['tuft']['n_soma_spikes'] == 0))
checks.append(('PT', 'BAC: paired -> >=3 spikes + tuft Ca spike >10ms',
               bac['both']['n_soma_spikes'], (3, 20),
               bac['both']['n_soma_spikes'] >= 3 and bac['both']['tuft_ms_above_m20'] > 10))

# ------------------------------------------------------------------ print
print('\n=== PHASE 1 GATE ===')
for cell, key, val, rng, ok in checks:
    vs = 'None' if val is None else f'{val:8.2f}'
    print(f"{'PASS' if ok else 'FAIL'}  {cell:4s} {key:48s} {vs:>9s}   target {rng}")
npass = sum(c[-1] for c in checks)
print(f'\n{npass}/{len(checks)} checks passed')
report['_checks'] = [dict(cell=c, check=k, value=v, target=list(r), passed=bool(o))
                     for c, k, v, r, o in checks]
with open('results/phase1_gate.json', 'w') as f:
    json.dump(report, f, indent=1, default=float)

# ------------------------------------------------------------------- plot
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharex=True, sharey=True)
for ax, idx, where in ((axes[0], 1, 'soma'), (axes[1], 2, 'tuft')):
    for cond, c, lab in (('soma', 'C0', 'soma pulse only'), ('tuft', 'C1', 'tuft input only'),
                         ('both', 'k', 'soma pulse + tuft input')):
        ax.plot(bactr[cond][0], bactr[cond][idx], color=c, lw=0.8, label=lab)
    ax.set_xlim(90, 180); ax.set_xlabel('time (ms)'); ax.set_title(f'PT BAC firing ({where} Vm)')
    ax.legend(fontsize=8)
axes[0].set_ylabel('Vm (mV)')
plt.tight_layout(); plt.savefig('results/phase1_gate.png', dpi=150)
print('saved results/phase1_gate.png')
sys.exit(0 if npass == len(checks) else 1)
