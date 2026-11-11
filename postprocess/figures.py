"""Reproduce the axes, normalization and panel layout of free-film-rupture."""

import json
from pathlib import Path

from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from plot_style import setup

ROOT = Path(__file__).resolve().parent


def axis_reference():
    return json.loads((ROOT/'data/figure_reference.json').read_text(encoding='utf-8'))


def apply_reference(ax, reference, fit_y=False):
    ax.set(xscale=reference['xscale'], yscale=reference['yscale'],
           xlabel=reference['xlabel'], ylabel=reference['ylabel'])
    ax.set_xticks(reference['xticks'])
    ax.set_xlim(reference['xlim'])
    if fit_y:
        ax.margins(y=.08)
        ax.relim()
        ax.autoscale(axis='y')
    else:
        ax.set_yticks(reference['yticks'])
        ax.set_ylim(reference['ylim'])
    ax.grid(False, which='both')


def save(fig, name):
    fig.savefig(ROOT/'figures'/f'{name}.png', dpi=600)
    fig.savefig(ROOT/'figures'/f'{name}.pdf')


def normal_mode():
    setup()
    reference = axis_reference()
    history = np.genfromtxt(ROOT/'data/rate_history.csv', names=True, delimiter=',',
                            dtype=None, encoding='utf-8')
    rates = json.loads((ROOT/'data/rates.json').read_text(encoding='utf-8'))
    with np.load(ROOT/'data/rate_surfaces.npz') as data:
        field = 1e3*(data['attraction_basilisk'][-1].T-1)
    fig, (field_ax, rate_ax) = plt.subplots(1, 2, figsize=(32/3, 6), layout='constrained')
    lo, hi = reference['normal_color_range']
    image = field_ax.imshow(field, extent=(0, 128, 0, 128), origin='lower',
                            cmap='RdBu_r', vmin=lo, vmax=hi, interpolation='nearest')
    fig.colorbar(image, ax=field_ax, shrink=.78, label=r'$10^3(h/h_0-1)$')
    field_ax.set_title(r'$a=0.0084375,\quad t^*=640$')
    apply_reference(field_ax, reference['normal_mode'][0])
    for index, case in enumerate(('attraction', 'capillary')):
        rate = next(row['continuum'] for row in rates if row['case'] == case)
        for solver, color, marker, start in (('JAM', '#D62728', 'o', 0),
                                              ('Basilisk', '#1F77B4', 's', 1)):
            use = (history['grid'] == 64) & (history['case'] == case) & (history['solver'] == solver)
            time, amplitude = history['time'][use], history['amplitude'][use]
            assert len(time) == 18 and time[-1] == 640
            fitted = np.polyfit(time, np.log(amplitude/amplitude[0]), 1)[0]
            assert abs(fitted/rate-1) < .01
            if solver == 'JAM':
                rate_ax.plot(time, rate*time, 'k-', label='Continuum' if index == 0 else None)
            rate_ax.plot(time[start::2], np.log(amplitude/amplitude[0])[start::2], marker,
                         color=color, markerfacecolor='white', label=solver if index == 0 else None)
    rate_ax.axhline(0., color='.65', linewidth=.7, zorder=0)
    rate_ax.text(.96, .57, r'$a=0.0084375$', ha='right', transform=rate_ax.transAxes)
    rate_ax.text(.96, .22, r'$a=0$', ha='right', transform=rate_ax.transAxes)
    rate_ax.set_box_aspect(1)
    apply_reference(rate_ax, reference['normal_mode'][1], fit_y=True)
    rate_ax.legend(loc='lower center', bbox_to_anchor=(.5, 1.02), ncols=3,
                   handlelength=1.15, handletextpad=.35, columnspacing=.8, borderaxespad=0.)
    save(fig, 'normal_mode_2d')
    return fig


def tcw():
    setup()
    reference = axis_reference()['tcw'][0]
    data = np.genfromtxt(ROOT/'data/tcw.csv', names=True, delimiter=',')
    assert len(data) == 2*4*84
    fig, ax = plt.subplots(figsize=(32/3, 6), layout='constrained')
    colors = ('#000000', '#D62728', '#1F77B4')
    times = (80, 320, 1280)
    for protocol, line in ((0, '-'), (1, '--')):
        for time, color in zip(times, colors):
            rows = data[(data['protocol'] == protocol) & (data['time'] == time)]
            q = rows['q']
            ax.plot(np.r_[0., q], 1e3*np.r_[0., rows['rms_continuum']], line, color=color)
            for field, fill, offset in (('rms_jam', 'white', 0), ('rms_basilisk', color, 1)):
                ax.plot(q, 1e3*rows[field], line, color=color, marker='o',
                        markerfacecolor=fill, markevery=(offset, 2), linewidth=.8)
    handles = [
        Line2D([], [], color='black', label='Continuous theory'),
        Line2D([], [], color='black', marker='o', markerfacecolor='white', label='JAM'),
        Line2D([], [], color='black', marker='o', label='Basilisk'),
        Line2D([], [], color='black', linestyle='--', label=r'$\varphi=0$ after first step'),
        *(Line2D([], [], color=color, label=rf'$t^*={t}$') for color, t in zip(colors, times)),
    ]
    ax.legend(handles=handles, loc='lower right', bbox_to_anchor=(1., 1.01), ncols=2,
              columnspacing=1.1, labelspacing=.25, borderaxespad=0.)
    apply_reference(ax, reference, fit_y=True)
    ax.set_ylim(bottom=0.)
    ax.tick_params(direction='in', top=True, right=True)
    save(fig, 'tcw_sipf_512')
    return fig


if __name__ == '__main__':
    plt.close(normal_mode())
    plt.close(tcw())
