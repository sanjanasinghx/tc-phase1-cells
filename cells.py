"""
Phase 1 cell library: reduced-morphology Hodgkin-Huxley cells for the
core/matrix thalamocortical model, written as NetPyNE cellParams.

Cortex : IT (RS pyramid, L2/3, L4, L5A, L6), PT (L5B, 4-compartment with an
         apical Ca hot zone for BAC firing), CT (L6 corticothalamic),
         PV (fast spiking), SOM (low-threshold spiking), VIP (adapting).
Thalamus: TC (relay; same biophysics for core and matrix, they differ only
          in projection pattern, set in Phase 2) and TRN (reticular).

All channels live in mod/ (hh2, km, itcal, hcur, cahva, cadyn, sk).
Run at celsius = 36.
"""

IONS = {'k': {'e': -90.0, 'i': 54.4, 'o': 2.5},
        'na': {'e': 50.0, 'i': 10.0, 'o': 140.0}}


def _sec(L, diam, mechs, nseg=1, Ra=100.0, cm=1.0, parent=None, parentX=1.0, ions=True):
    s = {'geom': {'L': L, 'diam': diam, 'nseg': nseg, 'Ra': Ra, 'cm': cm},
         'mechs': mechs}
    if ions:
        s['ions'] = {k: dict(v) for k, v in IONS.items()}
    if parent is not None:
        s['topol'] = {'parentSec': parent, 'parentX': parentX, 'childX': 0.0}
    return s


# ------------------------------------------------------------------ cortex
def rs_pyramid(gkm=4e-4, gl=1e-4, el=-70.0, vtraub=-55.0):
    """2-compartment regular-spiking pyramid (IT and CT)."""
    soma = _sec(25, 25, {'hh2': {'gnabar': 0.05, 'gkbar': 0.01, 'vtraub': vtraub},
                         'km': {'gbar': gkm},
                         'pas': {'g': gl, 'e': el}})
    dend = _sec(300, 3, {'pas': {'g': gl, 'e': el}}, nseg=3, parent='soma', parentX=0.0)
    return {'secs': {'soma': soma, 'dend': dend}}


def pt_cell(p=None):
    """4-compartment L5B PT neuron: soma, basal, apical trunk, tuft.
    Distal trunk/tuft carries the Ca hot zone (cahva + itcal + sk + cadyn) and
    a dense Ih, so coincident somatic AP + distal input gives BAC firing."""
    q = dict(gl=5e-5, el=-64.0, gkm=2e-4,
             g_hva_hot=0.004, g_t_hot=0.006, g_sk_hot=0.002,  # SK<0.002 -> runaway plateaus
             gh_soma=2e-5, gh_trunk=1e-4, gh_tuft=3e-4,
             gna_dend=0.004, gk_dend=0.002)
    if p:
        q.update(p)
    gl, el = q['gl'], q['el']
    soma = _sec(30, 30, {'hh2': {'gnabar': 0.05, 'gkbar': 0.01, 'vtraub': -60.0},
                         'km': {'gbar': q['gkm']},
                         'hcur': {'ghbar': q['gh_soma']},
                         'pas': {'g': gl, 'e': el}})
    basal = _sec(250, 4, {'hh2': {'gnabar': 0.005, 'gkbar': 0.0005, 'vtraub': -55.0},
                          'pas': {'g': gl, 'e': el}}, nseg=3, parent='soma', parentX=0.0)
    trunk = _sec(600, 3, {'hh2': {'gnabar': q['gna_dend'], 'gkbar': q['gk_dend'], 'vtraub': -55.0},
                          'hcur': {'ghbar': q['gh_trunk']},
                          'pas': {'g': gl, 'e': el}}, nseg=7, parent='soma', parentX=1.0)
    tuft = _sec(300, 2, {'hh2': {'gnabar': q['gna_dend'], 'gkbar': q['gk_dend'], 'vtraub': -55.0},
                         'hcur': {'ghbar': q['gh_tuft']},
                         'cahva': {'gbar': q['g_hva_hot']},
                         'itcal': {'gmax': q['g_t_hot'], 'vshift': 0.0},
                         'sk': {'gbar': q['g_sk_hot']},
                         'cadyn': {'taur': 80.0, 'gamma': 0.05},
                         'pas': {'g': gl, 'e': el}}, nseg=5, parent='apic_trunk', parentX=1.0)
    return {'secs': {'soma': soma, 'dend': basal, 'apic_trunk': trunk, 'apic_tuft': tuft}}


def pv_cell():
    """Fast-spiking basket cell: no adaptation, narrow spikes."""
    soma = _sec(35, 35, {'hh2': {'gnabar': 0.05, 'gkbar': 0.02, 'vtraub': -57.0},
                         'pas': {'g': 1.5e-4, 'e': -70.0}})
    return {'secs': {'soma': soma}}


def som_cell():
    """Low-threshold spiking (Martinotti-like): T current + strong adaptation."""
    soma = _sec(20, 20, {'hh2': {'gnabar': 0.05, 'gkbar': 0.01, 'vtraub': -55.0},
                         'km': {'gbar': 3e-4},
                         'itcal': {'gmax': 1.5e-3, 'vshift': 0.0},
                         'pas': {'g': 1e-4, 'e': -72.0}})
    dend = _sec(200, 2, {'pas': {'g': 1e-4, 'e': -72.0}}, nseg=3, parent='soma', parentX=0.0)
    return {'secs': {'soma': soma, 'dend': dend}}


def vip_cell():
    """Small high-input-resistance adapting interneuron."""
    soma = _sec(22, 22, {'hh2': {'gnabar': 0.05, 'gkbar': 0.01, 'vtraub': -55.0},
                         'km': {'gbar': 1.5e-4},
                         'pas': {'g': 1.5e-4, 'e': -66.0}})
    return {'secs': {'soma': soma}}


# ---------------------------------------------------------------- thalamus
def tc_cell(gT=0.002, gh=1e-5, gl=4e-5, el=-70.0):
    """Thalamocortical relay cell (Destexhe et al. 1998), single compartment."""
    soma = _sec(96, 96, {'hh2': {'gnabar': 0.09, 'gkbar': 0.01, 'vtraub': -50.0},
                         'itcal': {'gmax': gT, 'vshift': 0.0},
                         'hcur': {'ghbar': gh, 'eh': -40.0},
                         'pas': {'g': gl, 'e': el}})
    return {'secs': {'soma': soma}}


def trn_cell(gT=0.003, gl=5e-5, el=-75.0):
    """Thalamic reticular neuron: depolarized-shifted T current (ITs)."""
    soma = _sec(64.86, 70, {'hh2': {'gnabar': 0.2, 'gkbar': 0.02, 'vtraub': -55.0},
                            'itcal': {'gmax': gT, 'vshift': 5.0},
                            'pas': {'g': gl, 'e': el}})
    return {'secs': {'soma': soma}}


def all_cell_params():
    return {
        'IT':  rs_pyramid(),
        'CT':  rs_pyramid(gkm=2e-4, el=-72.0),
        'PT':  pt_cell(),
        'PV':  pv_cell(),
        'SOM': som_cell(),
        'VIP': vip_cell(),
        'TC':  tc_cell(),
        'TRN': trn_cell(),
    }
