#!/usr/bin/env python3
"""
Generate all paper figures for the whale thermoregulation manuscript.

Figures referenced by paper/main.tex:
  Fig 2  : Species thermoregulatory strategy portrait (radar chart)
  Fig 3  : TES vs water temperature profiles (all species)
  Fig 4  : Regional heat loss anatomy (stacked bars)
  Fig 5  : FA_min vulnerability heatmap (species x T_water)
  Fig 6  : Food-stress core temperature response
  Fig 7  : Scenario grid (all temperatures)

Optional diagnostic figure:
  Fig 1  : 3D FEM thermal field, not currently referenced in paper/main.tex

Usage:
  python scripts/postprocessing/generate_paper_figures.py [--species ...] [--fig 1 2 3 ...]
"""

import sys
import json
import argparse
import warnings
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize, BoundaryNorm, LogNorm
from matplotlib.cm import ScalarMappable
from pathlib import Path

PROJECT_DIR = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_DIR))

# ── species metadata ─────────────────────────────────────────────────────────

SPECIES = ['humpback', 'right_whale', 'blue_whale', 'bowhead', 'fin_whale', 'gray_whale', 'sei_whale']
LABELS  = ['Humpback', 'Right whale', 'Blue whale', 'Bowhead', 'Fin whale', 'Gray whale', 'Sei whale']
LABELS_SHORT = ['Hump.', 'Right', 'Blue', 'Bow.', 'Fin', 'Gray', 'Sei']
FAMILY  = ['Balaenopteridae', 'Balaenidae', 'Balaenopteridae',
           'Balaenidae', 'Balaenopteridae', 'Eschrichtiidae', 'Balaenopteridae']
SP_COLOR = {
    'humpback':    '#56b4e9',
    'right_whale': '#009e73',
    'blue_whale':  '#cc79a7',
    'bowhead':     '#d55e00',
    'fin_whale':   '#0072b2',
    'gray_whale':  '#e69f00',
    'sei_whale':   '#000000',
}
FAM_COLOR = {'Balaenopteridae': '#d73027', 'Balaenidae': '#4575b4',
             'Eschrichtiidae': '#66a61e'}
FAM_LS = {'Balaenopteridae': '-', 'Balaenidae': '--', 'Eschrichtiidae': ':'}

# Representative adult body mass (tonnes); sources in data/species/*.yaml
BODY_MASS_T = {
    'humpback':    30.0,
    'right_whale': 60.0,
    'blue_whale':  100.0,
    'bowhead':     80.0,
    'fin_whale':   70.0,
    'gray_whale':  25.0,
    'sei_whale':   20.0,
    'minke':        7.0,
}

# ── Minke validation species ─────────────────────────────────────────────────
# Minke is carried through the identical scenario grid but is the validation
# species (scaled sei-whale geometry, Folkow & Blix 1992 calibration).  It is
# overlaid on the comparative figures as a distinct grey dashed line where the
# metric is directly comparable (food, heat loss, core temperature).  For TES it
# is shown only as a control-law diagnostic, not as part of the comparative
# ranking, because the validation geometry is a scaled sei-whale proxy.
MINKE        = 'minke'
MINKE_LABEL  = 'Minke'
MINKE_COLOR  = '#555555'


def temp_label(value):
    """Format water temperatures so negative signs survive PDF text extraction."""
    return f'{value:g}'
MINKE_LS     = (0, (5, 2))
SP_COLOR['minke'] = MINKE_COLOR

# Thermal-maintenance metabolizable energy requirement.
# 1 dietary kcal = 4184 J; values are reported as Mcal/day for readability.
J_PER_MCAL = 4.184e6
def energy_mcal_per_day(Q_loss_W):
    return Q_loss_W * 86400.0 / J_PER_MCAL

OUTPUT_DIR = PROJECT_DIR / 'paper' / 'figures'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TEMPS = [-2, 0, 5, 10, 15, 20, 25, 28, 30]

# matplotlib style
plt.rcParams.update({
    'font.size': 9,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'font.family': 'serif',
    'axes.spines.top': False,
    'axes.spines.right': False,
})


# ── data helpers ─────────────────────────────────────────────────────────────

def load_all():
    d = {}
    for sp in SPECIES + [MINKE]:
        path = PROJECT_DIR / 'results' / sp / 'multicase_comparison_data.json'
        if not path.exists():
            print(f"  WARNING: no data for {sp} ({path}); skipping")
            continue
        d[sp] = json.loads(path.read_text())['cases']
    return d


def get_cases(cases, mode, fa=1.0, bci=1.0):
    return sorted(
        [c for c in cases
         if c['mode'] == mode
         and round(c.get('food_availability', 1.0), 2) == round(fa, 2)
         and round(c.get('body_condition_index', 1.0), 2) == round(bci, 2)],
        key=lambda x: x['T_water_C']
    )


def bm(case, key, default=np.nan):
    return case.get('biomarkers', {}).get(key, default)


def at_temp(cases, mode, T, fa=1.0, bci=1.0):
    hits = [c for c in get_cases(cases, mode, fa, bci)
            if abs(c['T_water_C'] - T) < 0.1]
    return hits[0] if hits else None


def minke_control_tes(data, T):
    """Return the current minke TES control-law value.

    The minke result export predates the epsilon-based TES post-processing and
    stores TES as zero.  Its Ttn and CCHE control ranges match the 10 C rorqual
    mapping, so the sei-whale control profile gives the correct diagnostic TES.
    """
    if 'sei_whale' not in data:
        return np.nan
    case = at_temp(data['sei_whale'], 'active', T)
    return bm(case, 'TES', np.nan) if case is not None else np.nan


# ── Figure 1: 3D FEM thermal field ───────────────────────────────────────────

def find_vtu(species, T_water, mode='active', fa=1.0, bci=1.0):
    results_dir = PROJECT_DIR / 'results' / species
    prefix = f"thermal_T{T_water:.0f}C_{mode}"
    if fa < 1.0 or bci < 1.0:
        prefix += f"_FA{fa:.2f}_BCI{bci:.2f}"
    candidates = sorted(results_dir.glob(f"{prefix}_*.vtu"))
    # Return the most recent FA=1.0/BCI=1.0 file (no suffix in name)
    plain = [f for f in candidates
             if '_FA' not in f.name and '_BCI' not in f.name]
    return plain[-1] if plain else (candidates[-1] if candidates else None)


def fig1_fem_thermal():
    """3D FEM thermal field: humpback at T=-2°C and T=30°C, side views."""
    try:
        import pyvista as pv
        pv.OFF_SCREEN = True
    except ImportError:
        print("  PyVista not available, skipping Fig 1")
        return

    cold_vtu = find_vtu('humpback', -2)
    warm_vtu = find_vtu('humpback', 30)
    if cold_vtu is None or warm_vtu is None:
        print(f"  VTU files not found for humpback cold/warm, skipping Fig 1")
        return

    print(f"  Cold VTU: {cold_vtu.name}")
    print(f"  Warm VTU: {warm_vtu.name}")

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.5))
    fig.subplots_adjust(left=0.02, right=0.88, top=0.93, bottom=0.05,
                        hspace=0.08, wspace=0.05)

    vmin, vmax = -2.0, 37.5
    cmap = 'coolwarm'
    norm = Normalize(vmin=vmin, vmax=vmax)

    configs = [
        (cold_vtu, f'Tw = {temp_label(-2)} °C (max vasoconstriction)', 'lateral'),
        (cold_vtu, f'Tw = {temp_label(-2)} °C (cross-section)',        'cross'),
        (warm_vtu, f'Tw = {temp_label(30)} °C (vasodilation)',         'lateral'),
        (warm_vtu, f'Tw = {temp_label(30)} °C (cross-section)',        'cross'),
    ]

    for idx, (vtu_path, title, view) in enumerate(configs):
        ax = axes[idx // 2][idx % 2]
        ax.set_title(title, pad=3)
        ax.set_xticks([])
        ax.set_yticks([])

        try:
            mesh = pv.read(str(vtu_path))
            T = mesh['Temperature_C']

            pl = pv.Plotter(off_screen=True, window_size=(600, 300))
            pl.set_background('white')
            pl.enable_anti_aliasing()

            if view == 'lateral':
                surface = mesh.extract_surface()
                pl.add_mesh(surface, scalars='Temperature_C', cmap=cmap,
                            clim=[vmin, vmax], show_scalar_bar=False,
                            smooth_shading=True, lighting=True)
                # lateral: look along Y axis
                pl.view_xz()
                pl.camera.zoom(1.1)
            else:
                # Cross-section through mid-body (y=0)
                bounds = mesh.bounds
                clipped = mesh.clip('y', origin=(0, 0, 0), invert=False)
                pl.add_mesh(clipped, scalars='Temperature_C', cmap=cmap,
                            clim=[vmin, vmax], show_scalar_bar=False,
                            lighting=False)
                pl.view_xz()
                pl.camera.zoom(1.1)

            screenshot = pl.screenshot(None, return_img=True)
            pl.close()
            ax.imshow(screenshot)
        except Exception as e:
            ax.text(0.5, 0.5, f'Render error:\n{e}',
                    ha='center', va='center', transform=ax.transAxes,
                    fontsize=7, color='red')

    # Shared colorbar
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.70])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label('Temperature (°C)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    fig.text(0.01, 0.97, 'a', fontsize=12, fontweight='bold', va='top')
    fig.text(0.01, 0.50, 'b', fontsize=12, fontweight='bold', va='top')

    out = OUTPUT_DIR / 'fig1_fem_thermal.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Figure 2: Strategy portrait (radar chart) ────────────────────────────────

def fig2_strategy_portrait(data):
    """Radar chart: 5 biomarkers at T=10°C for all species."""
    categories = ['|TES|\n(effort)', 'log η\n(efficiency)',
                  'TWF\n(windows)', 'log R_eff\n(insulation)',
                  'energy\n(per tonne)']
    N = len(categories)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(4.5, 4.5),
                           subplot_kw=dict(polar=True))
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)

    # Collect raw values
    raw = {}
    plot_species = list(SPECIES) + ([MINKE] if MINKE in data else [])
    plot_labels = list(LABELS) + ([MINKE_LABEL] if MINKE in data else [])

    for sp in plot_species:
        c = at_temp(data[sp], 'active', 10.0)
        if c is None:
            continue
        TES      = minke_control_tes(data, 10.0) if sp == MINKE else bm(c, 'TES', 0.0)
        eta      = bm(c, 'thermal_efficiency', 1.0)
        twf      = bm(c, 'thermal_window_fraction', 0.0)
        reff     = bm(c, 'R_eff_KperW', 1e-4)
        food_t   = energy_mcal_per_day(bm(c, 'Q_loss_W', 0.0)) / BODY_MASS_T[sp]
        raw[sp] = [TES, eta, twf, reff, food_t]

    # Normalise each dimension to [0,1] across species
    raw_arr = np.array([raw[sp] for sp in plot_species if sp in raw])
    # Radar portraits show effort magnitude.  Signed TES is retained in
    # the tables and TES profile figure; here |TES| makes high
    # vasoconstrictive effort visually extend outward.
    norm_TES   = np.abs(np.clip(raw_arr[:, 0], -1, 1))
    # log η: log10, then min-max
    log_eta    = np.log10(np.clip(raw_arr[:, 1], 0.01, None))
    norm_eta   = (log_eta - log_eta.min()) / (log_eta.max() - log_eta.min() + 1e-9)
    # TWF: already 0–1
    norm_twf   = raw_arr[:, 2]
    # log R_eff: log10, then min-max (higher = more insulation)
    log_reff   = np.log10(np.clip(raw_arr[:, 3], 1e-6, None))
    norm_reff  = (log_reff - log_reff.min()) / (log_reff.max() - log_reff.min() + 1e-9)
    # energy per tonne: higher = greater thermoregulatory energy cost
    food_col   = raw_arr[:, 4]
    norm_food  = (food_col - food_col.min()) / (food_col.max() - food_col.min() + 1e-9)

    normalised = np.column_stack([norm_TES, norm_eta, norm_twf, norm_reff, norm_food])

    for i, sp in enumerate(sp for sp in plot_species if sp in raw):
        vals = normalised[i].tolist() + [normalised[i, 0]]
        color = SP_COLOR[sp]
        label = plot_labels[plot_species.index(sp)]
        lw = 1.6 if sp == MINKE else 1.8
        ls = MINKE_LS if sp == MINKE else '-'
        alpha = 0.07 if sp == MINKE else 0.12
        ax.plot(angles, vals, color=color, linewidth=lw, linestyle=ls, label=label)
        ax.fill(angles, vals, color=color, alpha=alpha)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, size=8)
    ax.tick_params(axis='x', pad=14)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(['0.25', '0.5', '0.75', '1.0'], size=6, color='gray')
    ax.set_ylim(0, 1)
    ax.spines['polar'].set_visible(False)
    ax.grid(color='gray', linewidth=0.5, linestyle='--', alpha=0.5)

    # Outer ring labels for orientation
    ax.set_title('Thermoregulatory strategy portrait\n'
                 r'($T_w = 10\,°C$, active mode)',
                 pad=15, fontsize=9)

    handles = [mpatches.Patch(color=SP_COLOR[sp], label=plot_labels[i])
               for i, sp in enumerate(plot_species)]
    ax.legend(handles=handles, loc='upper right',
              bbox_to_anchor=(1.45, 1.10), frameon=False)

    out = OUTPUT_DIR / 'fig2_strategy_portrait.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Figure 3: TES profiles ───────────────────────────────────────────────────

def fig3_tes_profiles(data):
    """TES vs T_water for all species, active mode."""
    fig, ax = plt.subplots(figsize=(5.5, 3.35))

    for sp, lab in zip(SPECIES, LABELS):
        cases = get_cases(data[sp], 'active', fa=1.0)
        Tw  = [c['T_water_C'] for c in cases]
        tes = [bm(c, 'TES', np.nan) for c in cases]
        lw = 2.2 if sp in ('humpback', 'bowhead') else 1.5
        fam = FAMILY[SPECIES.index(sp)]
        ls = {'Balaenopteridae': '-', 'Balaenidae': '--',
              'Eschrichtiidae': ':'}.get(fam, '-')
        ax.plot(Tw, tes, color=SP_COLOR[sp], lw=lw, ls=ls,
                marker='o', ms=3.5, label=lab, zorder=3)

    if MINKE in data:
        cases = get_cases(data[MINKE], 'active', fa=1.0)
        Tw = [c['T_water_C'] for c in cases]
        tes = [minke_control_tes(data, t) for t in Tw]
        ax.plot(Tw, tes, color=MINKE_COLOR, lw=1.4, ls=MINKE_LS,
                marker='s', ms=3.2, label=MINKE_LABEL, zorder=2)

    ax.axhline(0, color='gray', lw=0.8, ls=':', zorder=1, label='Thermoneutral')
    ax.axhline(-1, color='#4575b4', lw=0.6, ls=':', alpha=0.5)
    ax.fill_between([-3, 31], [-1.05, -1.05], [-0.9, -0.9],
                    color='#4575b4', alpha=0.08, label='Vasoconstriction limit')

    ax.set_xlim(-4, 31)
    ax.set_ylim(-1.05, 0.05)
    ax.set_xlabel('Water temperature (°C)')
    ax.set_ylabel('Thermoregulatory Effort Score (TES)')
    ax.set_xticks(TEMPS)
    ax.set_xticklabels([temp_label(t) for t in TEMPS])
    ax.set_yticks([-1, -0.75, -0.5, -0.25, 0])
    ax.set_yticklabels(['−1\n(max vaso-\nconstriction)', '−0.75', '−0.5',
                        '−0.25', '0\n(thermo-\nneutral)'])

    # Family legend markers
    from matplotlib.lines import Line2D
    handles, hlabels = ax.get_legend_handles_labels()
    # Add family line-style note
    custom = [Line2D([0], [0], color='k', lw=1.5, ls='-',  label='Balaenopteridae'),
              Line2D([0], [0], color='k', lw=1.5, ls='--', label='Balaenidae'),
              Line2D([0], [0], color='k', lw=1.5, ls=':',  label='Eschrichtiidae')]
    n = len(SPECIES) + (1 if MINKE in data else 0)
    ax.legend(handles=handles[:n] + custom,
              labels=hlabels[:n] + ['Balaenopteridae', 'Balaenidae', 'Eschrichtiidae'],
              ncol=3, fontsize=6.7, frameon=False, loc='upper center',
              bbox_to_anchor=(0.5, -0.24), columnspacing=1.0, handlelength=1.8)

    ax.set_title('Thermoregulatory effort as a function of water temperature')
    out = OUTPUT_DIR / 'fig3_tes_profiles.pdf'
    fig.tight_layout(rect=[0, 0.14, 1, 1])
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Figure 4: Regional heat loss stacked bars ────────────────────────────────

REGION_MAP = {
    # raw region name  → display group, display colour
    'Body':                  ('Body trunk',     '#8ecae6'),
    'Peduncle':              ('Peduncle',        '#219ebc'),
    'Rostrum':               ('Rostrum',         '#023047'),
    'LeftPectoralFlipper':   ('Flippers',        '#fb8500'),
    'RightPectoralFlipper':  ('Flippers',        '#fb8500'),
    'Fluke':                 ('Flukes',          '#ffb703'),
    'DorsalFin':             ('Dorsal fin',      '#e76f51'),
    'VentralGrooves':        ('Ventral grooves', '#8338ec'),
}
GROUP_ORDER = ['Body trunk', 'Peduncle', 'Rostrum',
               'Flippers', 'Flukes', 'Dorsal fin', 'Ventral grooves']
GROUP_COLOR = {g: c for _, (g, c) in REGION_MAP.items()}


def fig4_regional_heat_loss(data):
    """Stacked bar chart: fractional heat loss by anatomical group, T=10°C."""
    # Comparative species plus the minke validation species (final bar).
    species = SPECIES + ([MINKE] if MINKE in data else [])
    labels  = LABELS + ([MINKE_LABEL] if MINKE in data else [])
    short_labels = LABELS_SHORT + (['Minke'] if MINKE in data else [])
    family  = FAMILY + (['Balaenopteridae'] if MINKE in data else [])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.6, 3.5),
                                    sharey=False)
    fig.subplots_adjust(wspace=0.35)

    x = np.arange(len(species))
    width = 0.6

    for ax, use_frac, ylabel in [
        (ax1, True,  'Fraction of total heat loss'),
        (ax2, False, 'Heat flux density (W m⁻²)'),
    ]:
        bottoms = np.zeros(len(species))
        patch_handles = {}

        for grp in GROUP_ORDER:
            vals = np.zeros(len(species))
            for si, sp in enumerate(species):
                c = at_temp(data[sp], 'active', 10.0)
                if c is None:
                    continue
                hf = c.get('heat_flux', {})
                Q_loss = c.get('total_heat_loss_W', 1.0)
                grp_W = 0.0
                grp_area = 0.0
                for raw_reg, (mapped_grp, _) in REGION_MAP.items():
                    if mapped_grp == grp and raw_reg in hf:
                        grp_W    += hf[raw_reg].get('total_W', 0.0)
                        grp_area += hf[raw_reg].get('area_m2', 0.0)
                if use_frac:
                    vals[si] = grp_W / Q_loss if Q_loss > 0 else 0
                else:
                    vals[si] = grp_W / grp_area if grp_area > 1e-6 else 0

            color = GROUP_COLOR.get(grp, '#aaa')
            bar = ax.bar(x, vals, width, bottom=bottoms, color=color,
                         label=grp, edgecolor='white', linewidth=0.4)
            patch_handles[grp] = bar
            bottoms = bottoms + vals

        ax.set_xticks(x)
        ax.set_xticklabels(short_labels, rotation=45, ha='right', fontsize=7)
        ax.set_ylabel(ylabel)
        ax.set_title('T$_w$ = 10°C, active mode')
        if use_frac:
            ax.set_ylim(0, 1.05)
            ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(xmax=1))
            ax.axhline(0.5, color='gray', lw=0.5, ls='--', zorder=0)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # Colour bars by family (minke shaded grey to mark the validation species)
        for si, sp in enumerate(species):
            fc = MINKE_COLOR if sp == MINKE else FAM_COLOR[family[si]]
            ax.axvspan(si - 0.4, si + 0.4, alpha=0.06 if sp == MINKE else 0.04,
                       color=fc, zorder=0)

    # Legend from ax1
    handles = [mpatches.Patch(color=GROUP_COLOR[g], label=g)
               for g in GROUP_ORDER]
    ax1.legend(handles=handles, fontsize=7, frameon=False,
               bbox_to_anchor=(0.0, -0.42), loc='lower left', ncol=4)

    # Family key: lifted above the per-axis subtitle to avoid overlapping it.
    for ax in [ax1, ax2]:
        ax.text(0.0, 1.16, '■ Balaenidae', color=FAM_COLOR['Balaenidae'],
                transform=ax.transAxes, ha='left', fontsize=7)
        ax.text(1.0, 1.16, '▲ Balaenopt.', color=FAM_COLOR['Balaenopteridae'],
                transform=ax.transAxes, ha='right', fontsize=7)

    fig.suptitle('Regional heat loss distribution', y=1.10, fontsize=10)
    out = OUTPUT_DIR / 'fig4_regional_heat_loss.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Figure 5: FA_min vulnerability heatmap ──────────────────────────────────

def fig5_famin_heatmap(data):
    """Heatmap: FA_min for each species × T_water."""
    # Build matrix
    mat = np.zeros((len(SPECIES), len(TEMPS)))
    for si, sp in enumerate(SPECIES):
        for ti, T in enumerate(TEMPS):
            c = at_temp(data[sp], 'active', T)
            if c:
                mat[si, ti] = bm(c, 'minimum_FA_to_balance', 1.0)

    fig, ax = plt.subplots(figsize=(6.5, 3.0))
    fig.subplots_adjust(left=0.18, right=0.87, top=0.90, bottom=0.22)

    # Use a custom colormap: green (safe) → yellow → red (critical)
    from matplotlib.colors import LinearSegmentedColormap
    cmap_custom = LinearSegmentedColormap.from_list(
        'vulnerability',
        [(0.0,  '#1a9850'),   # green: FA_min=0 (very safe)
         (0.35, '#fee08b'),   # yellow: FA_min=0.35
         (0.65, '#fc8d59'),   # orange: FA_min=0.65
         (0.99, '#d73027'),   # red: FA_min=0.99
         (1.0,  '#6e0000')],  # dark red: FA_min=1.0 (uncoupled)
    )

    im = ax.imshow(mat, cmap=cmap_custom, vmin=0, vmax=1,
                   aspect='auto', interpolation='nearest')

    # Contour lines at ecologically meaningful thresholds
    X = np.arange(len(TEMPS))
    Y = np.arange(len(SPECIES))
    Xg, Yg = np.meshgrid(X, Y)
    for level, ls, lw in [(0.5, '--', 1.0), (0.7, ':', 0.8)]:
        try:
            cs = ax.contour(Xg, Yg, mat, levels=[level],
                            colors='white', linewidths=lw, linestyles=ls)
            ax.clabel(cs, fmt=f'FA={level}', fontsize=6.5, colors='white',
                      inline=True)
        except Exception:
            pass

    # Annotate cells with value (only non-1.0 cells readable)
    for si in range(len(SPECIES)):
        for ti in range(len(TEMPS)):
            v = mat[si, ti]
            txt = '—' if v >= 0.999 else f'{v:.2f}'
            color = 'white' if v > 0.55 or v >= 0.999 else 'black'
            ax.text(ti, si, txt, ha='center', va='center',
                    fontsize=6.5, color=color, fontweight='normal')

    ax.set_xticks(range(len(TEMPS)))
    ax.set_xticklabels([temp_label(t) for t in TEMPS], fontsize=8)
    ax.set_yticks(range(len(SPECIES)))
    ax.set_yticklabels(LABELS, fontsize=8)
    ax.set_xlabel('Water temperature (°C)')
    ax.set_title(r'Minimum food availability ($\mathrm{FA}_\mathrm{min}$) for'
                 ' metabolic heat balance\n'
                 r'Active mode, FA = 1 baseline.  "—" = metabolic balance'
                 ' impossible (rorquals all temperatures)',
                 fontsize=8.5, pad=6)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label(r'$\mathrm{FA}_\mathrm{min}$', fontsize=8)
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    cbar.ax.tick_params(labelsize=7)

    # Family bracket on y-axis
    ax.axhline(1.5, color='gray', lw=0.5, ls=':')  # between right_whale and blue_whale
    ax.axhline(3.5, color='gray', lw=0.5, ls=':')  # between bowhead and fin_whale
    ax.text(-0.5, 0.5,  'Balaenopt.', ha='right', va='center', fontsize=6.5,
            color=FAM_COLOR['Balaenopteridae'], rotation=90)
    ax.text(-0.5, 2.5,  'Balaenidae', ha='right', va='center', fontsize=6.5,
            color=FAM_COLOR['Balaenidae'], rotation=90)
    ax.text(-0.5, 4.0,  'Balaenopt.', ha='right', va='center', fontsize=6.5,
            color=FAM_COLOR['Balaenopteridae'], rotation=90)

    out = OUTPUT_DIR / 'fig5_famin_heatmap.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Figure 6: Food-stress core temperature response ──────────────────────────

def fig6_food_stress(data):
    """T_max vs water temperature for 4 FA scenarios for each species."""
    # Scenarios restricted to the simulated (FA, BCI) grid: BCI in {1.0, 0.6} only.
    FA_CONFIGS = [
        (1.0, 1.0, 'Well-fed\n(FA=1.0, BCI=1.0)',     '#1a9850', '-',  2.0),
        (0.7, 1.0, 'Lean season\n(FA=0.7, BCI=1.0)',   '#fee08b', '--', 1.5),
        (0.5, 0.6, 'Food crisis\n(FA=0.5, BCI=0.6)',   '#fc8d59', '-.',  1.5),
        (0.3, 0.6, 'Starvation\n(FA=0.3, BCI=0.6)',    '#d73027', ':',  1.5),
    ]

    # Seven comparative species plus the minke validation species (8th panel).
    plot_species = list(SPECIES) + ([MINKE] if MINKE in data else [])
    plot_labels  = list(LABELS) + ([MINKE_LABEL] if MINKE in data else [])

    fig, axes = plt.subplots(4, 2, figsize=(7.2, 6.5), sharey=False)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.92, bottom=0.16,
                        hspace=0.72, wspace=0.30)
    axflat = axes.flatten()
    # Hide the unused 8th panel (7 species).
    for j in range(len(plot_species), len(axflat)):
        axflat[j].set_visible(False)

    for idx, (sp, lab) in enumerate(zip(plot_species, plot_labels)):
        ax = axflat[idx]
        ax.set_title(lab, pad=3)

        for fa, bci, label, color, ls, lw in FA_CONFIGS:
            cases = get_cases(data[sp], 'active', fa=fa, bci=bci)
            if not cases:
                continue
            Tw    = [c['T_water_C'] for c in cases]
            Tcore = [c.get('core_temp_C', c['T_max_C']) for c in cases]
            ax.plot(Tw, Tcore, color=color, lw=lw, ls=ls,
                    marker='o', ms=2.5, label=label if idx == 0 else '_nolegend_')

        # Target temperature line
        from scripts.utils.species_config import load_species_config
        try:
            cfg = load_species_config(sp)
            ax.axhline(cfg.core_temperature_C, color='gray', lw=0.7,
                       ls=':', label='T_target' if idx == 0 else '_nolegend_')
        except Exception:
            pass

        ax.set_xlabel('Water temp. (°C)', fontsize=8)
        ax.set_ylabel('Core temp. (°C)', fontsize=8)
        ax.set_xticks([-2, 5, 10, 15, 20, 25, 30])
        ax.tick_params(labelsize=7.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Shared legend
    handles = []
    for fa, bci, label, color, ls, lw in FA_CONFIGS:
        from matplotlib.lines import Line2D
        handles.append(Line2D([0], [0], color=color, lw=lw, ls=ls,
                               label=label.replace('\n', ' ')))
    from matplotlib.lines import Line2D
    handles.append(Line2D([0], [0], color='gray', lw=0.7, ls=':',
                          label='Physiological target'))
    fig.legend(handles=handles, loc='lower center', ncol=3,
               fontsize=7.2, frameon=False,
               bbox_to_anchor=(0.5, 0.02))

    fig.suptitle('Core temperature response to food-availability stress',
                 fontsize=10, y=0.985)

    out = OUTPUT_DIR / 'fig6_food_stress.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Figure 7 (scenario grid): regional heat loss grid, all temperatures ──────

def fig_scenario_grid(data):
    """
    3×3 grid of stacked bar charts (one per T_water).
    Each panel: all species, stacked regional kW.
    Core temp annotated above each bar; flipper-perfusion dot below.
    Saved as paper/figures/scenario_grid.pdf
    """
    # Comparative species plus the minke validation species (final bar).
    species = SPECIES + ([MINKE] if MINKE in data else [])
    labels  = LABELS + ([MINKE_LABEL] if MINKE in data else [])
    short_labels = LABELS_SHORT + (['Minke'] if MINKE in data else [])

    ncols, nrows = 3, 3
    fig, axes = plt.subplots(nrows, ncols, figsize=(10.5, 8.5),
                             sharey=False)
    fig.subplots_adjust(hspace=0.55, wspace=0.35)

    x = np.arange(len(species))
    width = 0.65

    for ti, T in enumerate(TEMPS):
        ax = axes[ti // ncols][ti % ncols]
        ax.set_title(f'Tw = {temp_label(T)}°C', fontsize=8, pad=2)

        bottoms = np.zeros(len(species))
        Q_totals = np.zeros(len(species))
        core_temps = [np.nan] * len(species)
        flipper_perf = [np.nan] * len(species)

        for grp in GROUP_ORDER:
            vals = np.zeros(len(species))
            for si, sp in enumerate(species):
                c = at_temp(data[sp], 'active', T)
                if c is None:
                    continue
                hf = c.get('heat_flux', {})
                Q_loss = c.get('total_heat_loss_W', 1.0)
                Q_totals[si] = Q_loss / 1000  # kW
                core_temps[si] = c.get('core_temp_C', np.nan)
                pm = c.get('perfusion_scale_map', {})
                flipper_perf[si] = pm.get('LeftPectoralFlipper',
                                   pm.get('RightPectoralFlipper', np.nan))

                grp_W = sum(
                    hf[r].get('total_W', 0.0)
                    for r, (g, _) in REGION_MAP.items()
                    if g == grp and r in hf
                ) / 1000  # kW
                vals[si] = grp_W

            ax.bar(x, vals, width, bottom=bottoms,
                   color=GROUP_COLOR.get(grp, '#aaa'),
                   edgecolor='white', linewidth=0.3, zorder=3)
            bottoms += vals

        # Headroom so the core-temp labels do not collide with the panel title
        _qmax = max(Q_totals) if np.isfinite(Q_totals).any() and max(Q_totals) > 0 else 1.0
        ax.set_ylim(top=_qmax * 1.22)

        # Annotations: core temperature above each bar
        for si in range(len(species)):
            if np.isfinite(core_temps[si]) and np.isfinite(Q_totals[si]):
                ax.text(x[si], Q_totals[si] + 0.02 * max(Q_totals) + 0.5,
                        f'{core_temps[si]:.1f}',
                        ha='center', va='bottom', fontsize=5.5,
                        rotation=0, color='#222')

        # Flipper-perfusion dots just below x-axis
        ax_ymin = ax.get_ylim()[0]
        for si in range(len(species)):
            if np.isfinite(flipper_perf[si]):
                # dot area ∝ perfusion scale (min→5pt, max→40pt)
                p = flipper_perf[si]
                sz = 5 + 120 * (p / 3.0)  # 3.0 = max perfusion scale
                ax.scatter(x[si], -0.05 * max(Q_totals) - 0.5, s=sz,
                           color=SP_COLOR[species[si]], zorder=4, clip_on=False)

        ax.set_xticks(x)
        ax.set_xticklabels(
            short_labels,
            fontsize=5.5
        )
        ax.set_ylabel('Heat loss (kW)', fontsize=6.5)
        ax.tick_params(axis='y', labelsize=6)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(bottom=-0.08 * max(Q_totals) - 1)

    # Shared region legend
    handles = [mpatches.Patch(color=GROUP_COLOR[g], label=g)
               for g in GROUP_ORDER]
    fig.legend(handles=handles, loc='lower center', ncol=4,
               fontsize=7, frameon=False, bbox_to_anchor=(0.5, -0.02))

    # Overall title
    fig.suptitle(
        'Regional heat loss by water temperature — active thermoregulatory control\n'
        '(numbers = core temp °C; dot size ∝ flipper perfusion scale)',
        fontsize=9, y=1.005
    )

    out = OUTPUT_DIR / 'scenario_grid.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Figure 8: mass-specific energy requirement relative to body size ──────────

def fig_food_per_mass(data):
    """
    Thermal-maintenance energy requirement normalised by body size, versus water
    temperature, for all seven species (active, FA=1, BCI=1).  Two size
    normalisations, both plotted against temperature:
      (a) per unit body mass  (Mcal / tonne / day);
      (b) per unit body surface area (Mcal / m^2 / day) -- the physically
          meaningful normalisation for a surface heat-loss process, equivalent
          to the mean surface heat flux expressed in energy units.
    (Per unit body *volume* is not shown: for near-neutrally-buoyant whales it
    is proportional to the per-mass panel, so it carries no additional
    information.)
    Saved as paper/figures/food_per_mass.pdf
    """
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.2, 3.4))

    for sp, lab, fam in zip(SPECIES, LABELS, FAMILY):
        cases = get_cases(data[sp], 'active', fa=1.0, bci=1.0)
        Ts = [c['T_water_C'] for c in cases]
        food = [energy_mcal_per_day(bm(c, 'Q_loss_W')) for c in cases]
        area = [bm(c, 'surface_area_m2') for c in cases]
        y_mass = [f / BODY_MASS_T[sp] for f in food]
        y_area = [f / a for f, a in zip(food, area)]
        axA.plot(Ts, y_mass, color=SP_COLOR[sp], ls=FAM_LS[fam], lw=1.8,
                 marker='o', ms=3, label=lab)
        axB.plot(Ts, y_area, color=SP_COLOR[sp], ls=FAM_LS[fam], lw=1.8,
                 marker='o', ms=3, label=lab)

    # Minke validation overlay (distinct grey dashed line)
    if MINKE in data:
        mc = get_cases(data[MINKE], 'active', fa=1.0, bci=1.0)
        Ts = [c['T_water_C'] for c in mc]
        food = [energy_mcal_per_day(bm(c, 'Q_loss_W')) for c in mc]
        area = [bm(c, 'surface_area_m2') for c in mc]
        axA.plot(Ts, [f / BODY_MASS_T[MINKE] for f in food], color=MINKE_COLOR,
                 ls=MINKE_LS, lw=1.8, marker='s', ms=3, label=MINKE_LABEL)
        axB.plot(Ts, [f / a for f, a in zip(food, area)], color=MINKE_COLOR,
                 ls=MINKE_LS, lw=1.8, marker='s', ms=3, label=MINKE_LABEL)

    axA.set_xlabel('Water temperature $T_\\infty$ (°C)')
    axA.set_ylabel('Energy requirement\n(Mcal tonne$^{-1}$ day$^{-1}$)')
    axA.set_title('(a) Per unit body mass', fontsize=9)
    axA.grid(True, alpha=0.25, lw=0.5)
    axA.set_ylim(bottom=0)
    axA.legend(fontsize=6.5, loc='upper right', frameon=False, ncol=1)

    axB.set_xlabel('Water temperature $T_\\infty$ (°C)')
    axB.set_ylabel('Energy requirement\n(Mcal m$^{-2}$ day$^{-1}$)')
    axB.set_title('(b) Per unit body surface area', fontsize=9)
    axB.grid(True, alpha=0.25, lw=0.5)
    axB.set_ylim(bottom=0)

    fig.tight_layout()
    out = OUTPUT_DIR / 'food_per_mass.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Figure: absolute thermal-maintenance energy requirement vs temperature ────

def fig_food_vs_temperature(data):
    """
    Absolute thermal-maintenance energy requirement (Mcal per day) versus water
    temperature for the seven comparative species, with the minke validation
    species overlaid as a distinct grey dashed line.  Active mode, FA=1, BCI=1.
    Saved as paper/figures/food_vs_temperature.pdf
    """
    fig, ax = plt.subplots(figsize=(5.5, 3.6))

    for sp, lab, fam in zip(SPECIES, LABELS, FAMILY):
        cases = get_cases(data[sp], 'active', fa=1.0, bci=1.0)
        Ts   = [c['T_water_C'] for c in cases]
        food = [energy_mcal_per_day(bm(c, 'Q_loss_W')) for c in cases]
        ax.plot(Ts, food, color=SP_COLOR[sp], ls=FAM_LS[fam], lw=1.8,
                marker='o', ms=3, label=lab)

    if MINKE in data:
        mc = get_cases(data[MINKE], 'active', fa=1.0, bci=1.0)
        Ts   = [c['T_water_C'] for c in mc]
        food = [energy_mcal_per_day(bm(c, 'Q_loss_W')) for c in mc]
        ax.plot(Ts, food, color=MINKE_COLOR, ls=MINKE_LS, lw=1.8,
                marker='s', ms=3, label=MINKE_LABEL)

    ax.set_xlabel('Water temperature $T_\\infty$ (°C)')
    ax.set_ylabel('Energy requirement\n(Mcal day$^{-1}$)')
    ax.set_xticks(TEMPS)
    ax.set_xticklabels([temp_label(t) for t in TEMPS], rotation=30, ha='right')
    ax.tick_params(axis='x', labelsize=7)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.25, lw=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=6.5, loc='upper right', frameon=False, ncol=2)
    ax.set_title('Thermal-maintenance energy requirement '
                 '(active mode, FA = 1)', fontsize=9)

    fig.tight_layout()
    out = OUTPUT_DIR / 'food_vs_temperature.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Figure 9: strategy portrait across all temperatures (radar small multiples) ─

def fig_strategy_portrait_grid(data):
    """
    3x3 grid of radar charts (one per water temperature) showing how the
    thermoregulatory strategy portrait shifts with environment.  Same four
    biomarkers as Fig. 2 (TES, log eta, TWF, log R_eff), all seven species.
    Each axis is normalised GLOBALLY across every temperature and species so the
    panels are directly comparable.  Saved as figures/strategy_portrait_grid.pdf
    """
    cats = [r'$|$TES$|$', r'log $\eta$', 'TWF', r'log $R_\mathrm{eff}$', 'energy/t']
    N = len(cats)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]

    # ---- collect raw values for every (temp, species) ----------------------
    plot_species = list(SPECIES) + ([MINKE] if MINKE in data else [])
    plot_labels = list(LABELS) + ([MINKE_LABEL] if MINKE in data else [])
    plot_family = list(FAMILY) + (['Balaenopteridae'] if MINKE in data else [])

    raw = {}  # (T, sp) -> [TES, eta, twf, reff, food_per_tonne]
    for T in TEMPS:
        for sp in plot_species:
            c = at_temp(data[sp], 'active', T)
            if c is None:
                continue
            TES = minke_control_tes(data, T) if sp == MINKE else bm(c, 'TES', 0.0)
            raw[(T, sp)] = [TES,
                            bm(c, 'thermal_efficiency', 1.0),
                            bm(c, 'thermal_window_fraction', 0.0),
                            bm(c, 'R_eff_KperW', 1e-4),
                            energy_mcal_per_day(bm(c, 'Q_loss_W', 0.0)) / BODY_MASS_T[sp]]

    all_vals = np.array(list(raw.values()))
    log_eta_all = np.log10(np.clip(all_vals[:, 1], 0.01, None))
    log_reff_all = np.log10(np.clip(all_vals[:, 3], 1e-6, None))
    eta_lo, eta_hi = log_eta_all.min(), log_eta_all.max()
    reff_lo, reff_hi = log_reff_all.min(), log_reff_all.max()
    food_lo, food_hi = all_vals[:, 4].min(), all_vals[:, 4].max()

    def normalise(vec):
        tes = abs(np.clip(vec[0], -1, 1))
        eta = (np.log10(max(vec[1], 0.01)) - eta_lo) / (eta_hi - eta_lo + 1e-9)
        twf = np.clip(vec[2], 0, 1)
        reff = (np.log10(max(vec[3], 1e-6)) - reff_lo) / (reff_hi - reff_lo + 1e-9)
        food = (vec[4] - food_lo) / (food_hi - food_lo + 1e-9)
        return [tes, eta, twf, reff, food]

    fig, axes = plt.subplots(3, 3, figsize=(8.5, 9.0),
                             subplot_kw=dict(polar=True))
    for ti, T in enumerate(TEMPS):
        ax = axes[ti // 3][ti % 3]
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        for i, sp in enumerate(plot_species):
            if (T, sp) not in raw:
                continue
            vals = normalise(raw[(T, sp)])
            vals = vals + vals[:1]
            ls = MINKE_LS if sp == MINKE else FAM_LS[plot_family[i]]
            lw = 1.2 if sp == MINKE else 1.3
            ax.plot(angles, vals, color=SP_COLOR[sp], lw=lw, ls=ls)
            ax.fill(angles, vals, color=SP_COLOR[sp],
                    alpha=0.04 if sp == MINKE else 0.06)
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(cats, size=6.5)
        ax.set_yticks([0.5, 1.0])
        ax.set_yticklabels([], size=5)
        ax.set_ylim(0, 1)
        ax.spines['polar'].set_visible(False)
        ax.grid(color='gray', linewidth=0.4, linestyle='--', alpha=0.5)
        ax.set_title(f'Tw = {temp_label(T)} °C', fontsize=9, pad=8)

    handles = [mpatches.Patch(color=SP_COLOR[sp], label=plot_labels[i])
               for i, sp in enumerate(plot_species)]
    fig.legend(handles=handles, loc='lower center', ncol=8 if MINKE in data else 7,
               frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle('Thermoregulatory strategy portrait across water temperatures '
                 '(active mode; axes normalised globally)', fontsize=10, y=0.995)
    fig.tight_layout(rect=[0, 0.03, 1, 0.98])
    out = OUTPUT_DIR / 'strategy_portrait_grid.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ── Figure: simulation design infographic ────────────────────────────────────

def fig_simulation_design(data=None):
    """
    Schematic of the simulation matrix: 99 cases per whale model (9 water
    temperatures x [4 thermoregulatory modes + 7 food/condition combinations]),
    across 8 whale models (7 comparative + minke validation), plus one minke
    Folkow-Blix calibration run.  793 finite-element simulations total.
    Saved as paper/figures/simulation_design.pdf
    """
    import matplotlib.patches as mp
    temps = [-2, 0, 5, 10, 15, 20, 25, 28, 30]
    nrows = len(temps)                 # 9 temperatures
    mode_cols = 4                      # 4 modes (baseline FA=1, BCI=1)
    grid_cols = 7                      # 7 additional FA x BCI points (active)
    ncols = mode_cols + grid_cols      # 11 -> 99 cells per model

    C_MODE = '#8ec9e6'   # thermoregulatory modes (baseline)
    C_GRID = '#f4a86a'   # food-availability x body-condition grid (active)

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.set_xlim(-2.6, ncols + 0.5)
    ax.set_ylim(-3.4, nrows + 2.0)
    ax.axis('off')

    # --- per-model 9x11 grid (99 cells) ---
    for r in range(nrows):
        for c in range(ncols):
            col = C_MODE if c < mode_cols else C_GRID
            ax.add_patch(mp.Rectangle((c, r), 0.9, 0.9, facecolor=col,
                                      edgecolor='white', lw=0.6))
        ax.text(-0.35, r + 0.45, f'{temps[r]}', ha='right', va='center',
                fontsize=7, color='0.25')
    ax.text(-2.4, nrows / 2, 'water temperature (°C)', rotation=90,
            ha='center', va='center', fontsize=8, color='0.25')

    # column-group brackets
    def bracket(x0, x1, y, label, color):
        ax.plot([x0, x0, x1, x1], [y, y + 0.18, y + 0.18, y], color=color, lw=1.0)
        ax.text((x0 + x1) / 2, y + 0.35, label, ha='center', va='bottom',
                fontsize=7.5, color=color)
    bracket(0, mode_cols - 0.1, nrows + 0.15,
            '4 thermoregulatory modes\n(FA=1, BCI=1)', '#2b6f9e')
    bracket(mode_cols, ncols - 0.1, nrows + 0.15,
            '7 food × condition combinations\n(active mode)', '#c56a25')

    ax.text(ncols / 2, nrows + 1.55, '99 finite-element cases per whale model',
            ha='center', va='center', fontsize=9.5, fontweight='bold')

    # --- x8 models strip ---
    y = -1.3
    ax.text(ncols / 2, y + 0.75, r'$\times$ 8 whale models',
            ha='center', fontsize=9.5, fontweight='bold')
    labels = ['Humpback', 'Right', 'Blue', 'Bowhead', 'Fin', 'Gray', 'Sei', 'Minke']
    cols = [SP_COLOR['humpback'], SP_COLOR['right_whale'], SP_COLOR['blue_whale'],
            SP_COLOR['bowhead'], SP_COLOR['fin_whale'], SP_COLOR['gray_whale'],
            SP_COLOR['sei_whale'], '#7f7f7f']
    x0 = (ncols - 8 * 1.35) / 2
    for i, (lab, cc) in enumerate(zip(labels, cols)):
        xx = x0 + i * 1.35
        hatch = '////' if lab == 'Minke' else None
        ax.add_patch(mp.Rectangle((xx, y), 1.15, 0.5, facecolor=cc, alpha=0.85,
                                  edgecolor='0.3', lw=0.5, hatch=hatch))
        ax.text(xx + 0.575, y + 0.25, lab, ha='center', va='center',
                fontsize=6.3, color='white' if lab != 'Sei' else 'white')
    ax.text(ncols / 2, y - 0.35, '7 comparative species  +  minke (validation, hatched)',
            ha='center', va='center', fontsize=7, color='0.35')

    # --- totals ---
    ax.text(ncols / 2, -2.7,
            r'$8 \times 99 = \mathbf{792}$ grid simulations $\;+\;$ '
            r'1 minke Folkow–Blix validation run $\;=\;$ '
            r'$\mathbf{793}$ finite-element simulations',
            ha='center', va='center', fontsize=9)

    fig.tight_layout()
    out = OUTPUT_DIR / 'simulation_design.pdf'
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {out}")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Generate paper figures')
    parser.add_argument('--fig', nargs='+', type=int, default=[2, 3, 4, 6, 7, 8, 9, 10, 11],
                        help='Which figures to generate (default: paper-referenced figures; 1 and 5 are optional/retired)')
    args = parser.parse_args()

    print("Loading data...")
    data = load_all()
    print(f"  Loaded {sum(len(v) for v in data.values())} cases across "
          f"{len(data)} species")

    figs = {
        1: ('3D FEM thermal field',               lambda: fig1_fem_thermal()),
        2: ('Strategy portrait (radar)',           lambda: fig2_strategy_portrait(data)),
        3: ('TES profiles',                        lambda: fig3_tes_profiles(data)),
        4: ('Regional heat loss stacked bars',     lambda: fig4_regional_heat_loss(data)),
        5: ('FA_min vulnerability heatmap',        lambda: fig5_famin_heatmap(data)),
        6: ('Food-stress core temp response',      lambda: fig6_food_stress(data)),
        7: ('Scenario grid (all temperatures)',    lambda: fig_scenario_grid(data)),
        8: ('Energy requirement per body mass',    lambda: fig_food_per_mass(data)),
        9: ('Strategy portrait grid (all temps)',  lambda: fig_strategy_portrait_grid(data)),
        10: ('Simulation design infographic',       lambda: fig_simulation_design()),
        11: ('Energy requirement vs temperature',   lambda: fig_food_vs_temperature(data)),
    }

    for n in sorted(args.fig):
        if n not in figs:
            print(f"  Unknown figure {n}")
            continue
        desc, fn = figs[n]
        print(f"\nFig {n}: {desc}")
        try:
            fn()
        except Exception as e:
            import traceback
            print(f"  ERROR: {e}")
            traceback.print_exc()

    print("\nDone. Figures in:", OUTPUT_DIR)


if __name__ == '__main__':
    main()
