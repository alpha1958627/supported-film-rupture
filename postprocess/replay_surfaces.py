"""Replay seed 61100 at the same 82 times and 256 h0 window as free-film-rupture."""

import argparse
from pathlib import Path
import tempfile

import numpy as np

from run_tcw import WallFilm2D, run_pair, shell_rows


def replay(executable, output, data):
    model = WallFilm2D(512, 1024., 8., theta=1e-5)
    steps = np.array([0, 1, *range(2, 161, 2)])
    selected = np.searchsorted(steps, [1, 10, 40, 160])
    windows, power, errors = {}, [], []
    with tempfile.TemporaryDirectory(dir=output, prefix='replay-') as scratch:
        for name, noisy in (('continuous', 160), ('first_step', 1)):
            jam, basilisk, metrics = run_pair(executable, model, 160, noisy, 61100,
                                              scratch, 2)
            assert len(jam) == len(steps) == 82
            assert np.max(abs(jam.mean(axis=(1, 2))-1.)) < 1e-12
            power.append(abs(np.fft.fft2(np.stack((jam[selected], basilisk[selected]))-1.)
                             /model.n**2)**2)
            windows[name] = basilisk[:, :128, :128].copy()
            windows['jam_'+name] = jam[:, :128, :128].copy()
            errors.append(metrics['max_paired_error'])
            print(name, 'replayed, max paired field difference', errors[-1], flush=True)
    np.testing.assert_array_equal(windows['continuous'][:2], windows['first_step'][:2])
    rows = np.asarray(shell_rows(model, np.asarray(power), 8.*np.array([1, 10, 40, 160])))
    with np.load(data/'tcw_seed_power.npz') as archive:
        assert archive['seeds'][0] == 61100
        np.testing.assert_array_equal(rows[:, :5], archive['shell_rows'])
        expected = np.sqrt(archive['power'][0])
    np.testing.assert_allclose(rows[:, 5:], expected, rtol=2e-9, atol=1e-17)
    np.savez_compressed(data/'tcw_surface_frames.npz', **windows, seed=61100,
                        x=(np.arange(128)+.5)*2, t=8.*steps,
                        max_pair_error=max(errors))
    print('Replay matched all 672 archived shell rows.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    replay(args.executable.resolve(), root.parent/'build', root/'data')
