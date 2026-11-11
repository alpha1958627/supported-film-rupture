"""Render two thickness animations with a shared scale and camera."""

from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
from matplotlib import colormaps
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np
from PIL import Image


def style_surface(ax, limit):
    ax.set(xlim=(0, 256), ylim=(0, 256), zlim=(-limit, limit),
           xticks=(0, 128, 256), yticks=(128, 256), zticks=(-1, 0, 1),
           xlabel=r'$x/h_0$', ylabel=r'$y/h_0$', zlabel=r'$100(h/h_0-1)$')
    ax.zaxis.set_major_formatter(StrMethodFormatter('{x:g}'))
    ax.tick_params(axis='z', pad=7)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((.91, .93, .95, 1))
        axis.labelpad = 16
    ax.set_box_aspect((1, 1, .46))
    ax.view_init(elev=30, azim=-125)
    ax.grid(False)


def render(output, name, title, x, t, fields, limit):
    paired = len(fields) == 2
    fig = plt.figure(figsize=(15, 6) if paired else (8, 6), dpi=110, facecolor='#f5f7f9')
    layout = ((.06, .05, .40, .77), (.48, .05, .40, .77)) if paired else ((.12, .07, .70, .82),)
    axes = [fig.add_axes(bounds, projection='3d', facecolor=fig.get_facecolor()) for bounds in layout]
    for ax in axes:
        style_surface(ax, limit)
    if paired:
        for center, label in zip((.26, .68), ('JAM', 'Basilisk')):
            fig.text(center, .82, label, ha='center')

    palette = colormaps['viridis']
    norm = Normalize(-1., 1.)
    colorbar = fig.add_axes([.905, .29, .016, .43] if paired else [.84, .29, .023, .43])
    fig.colorbar(ScalarMappable(norm=norm, cmap=palette), cax=colorbar, ticks=(-1, -.5, 0, .5, 1),
                 format=StrMethodFormatter('{x:g}'), extend='both')
    heading = fig.suptitle('', y=.94)
    X, Y = np.meshgrid(x, x, indexing='ij')
    frames = []
    for time, values in zip(t, zip(*fields)):
        heading.set_text(rf'{title}, $t^*={time:g}$')
        surfaces = [ax.plot_surface(X, Y, field, rcount=len(x), ccount=len(x),
                                    cmap=palette, norm=norm, linewidth=0,
                                    antialiased=False, rasterized=True)
                    for ax, field in zip(axes, values)]
        fig.canvas.draw()
        frame = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert('RGB')
        frames.append(frame.quantize(colors=256, dither=Image.Dither.NONE))
        if not paired and time == t[-1]:
            fig.savefig(output / (name + '.png'), dpi=600, facecolor=fig.get_facecolor())
            fig.savefig(output / (name + '.pdf'), dpi=600, facecolor=fig.get_facecolor())
        for surface in surfaces:
            surface.remove()
    duration = [450] + [100]*(len(frames)-2) + [900]
    path = output / (name + '.gif')
    frames[0].save(path, save_all=True, append_images=frames[1:], duration=duration, loop=0, optimize=False)
    with Image.open(path) as gif:
        assert gif.n_frames == len(t) and gif.info['loop'] == 0
    plt.close(fig)
    print(f'{path.name}: {len(t)} frames, t*={t[0]:g} to {t[-1]:g}', flush=True)


def main():
    from plot_style import setup
    setup()
    root = Path(__file__).resolve().parent
    with np.load(root/'data/tcw_surface_frames.npz') as archive:
        x, t = archive['x'], archive['t']
        cases = {'continuous': 'Continuous noise', 'first_step': 'First-step noise only'}
        fields = {name: (100*(archive['jam_'+name]-1.), 100*(archive[name]-1.))
                  for name in cases}
    for name, title in cases.items():
        render(root/'figures', f'tcw_{name}_3d_paired', title, x, t, fields[name], 1.75)


if __name__ == '__main__':
    main()
