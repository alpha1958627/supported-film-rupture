"""Validate the supported-film midpoint solver, then run paired TCW ensembles."""

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from time import perf_counter

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'jam'))
from wall_film_2d import WallFilm2D


def run_pair(exe, model, steps, noisy, seed, scratch, stride, *, initial=None):
    n = model.n
    h = np.ones((n, n)) if initial is None else np.asarray(initial, dtype=float)
    if h.shape != (n, n):
        raise ValueError(f'expected initial height of shape {(n, n)}')
    saved_steps = sorted(set(range(0, steps+1, stride)) | {1, steps})
    save = set(saved_steps)
    frames = [h]
    jam_seconds = 0.
    preparation_seconds = 0.
    rng = np.random.default_rng(seed)
    with tempfile.TemporaryFile(dir=scratch) as file:
        h.astype('<f8', copy=False).tofile(file)
        for j in range(steps):
            start = perf_counter()
            white = rng.standard_normal(h.shape) if j < noisy else None
            if white is not None:
                white.astype('<f8', copy=False).tofile(file)
            preparation_seconds += perf_counter()-start
            start = perf_counter()
            h = model.advance(h, white)
            jam_seconds += perf_counter()-start
            if j+1 in save:
                frames.append(h)
        file.seek(0)
        start = perf_counter()
        args = (n, model.length, model.attraction, model.theta, model.dt, steps, stride, noisy)
        process = subprocess.run([str(exe), *map(str, args)], stdin=file,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if process.returncode:
            raise RuntimeError(f'Basilisk exited with {process.returncode}: '
                               f'{process.stderr.decode(errors="replace")}')
        basilisk_seconds = perf_counter()-start
    jam = np.asarray(frames)
    basilisk = np.frombuffer(process.stdout, dtype='<f8').reshape(jam.shape)
    error = float(np.max(abs(jam-basilisk)))
    assert error < 1e-10
    return jam, basilisk, dict(jam_seconds=jam_seconds, basilisk_seconds=basilisk_seconds,
                              preparation_seconds=preparation_seconds, max_paired_error=error)


def continuous_rate(model):
    return (model.q**2/3.)*(3*model.attraction-model.q**2)


def shell_rows(model, power, times):
    """The main theory uses the continuum spatial and temporal symbol."""
    q = model.q
    rate = continuous_rate(model)
    decay = -np.expm1(2*times[:, None, None]*rate)
    equilibrium = np.divide(model.theta, model.length**2*q**2,
                            out=np.zeros_like(q), where=q > 0)
    theory_on = equilibrium*decay
    theory_off = power[1, 0, 0]*np.exp(2*(times[:, None, None]-model.dt)*rate)
    shell = np.rint(q*model.length/(2*np.pi)).astype(int)
    max_shell = int(np.floor((np.pi/6)*model.length/(2*np.pi)-.5))
    counts = np.bincount(shell.ravel())
    rows = []
    for p, theory in enumerate((theory_on, theory_off)):
        for i, t in enumerate(times):
            spectra = [np.bincount(shell.ravel(), weights=field.ravel())
                       for field in (theory[i], *power[p, :, i])]
            for m in range(1, max_shell+1):
                rows.append((p, t, m, counts[m], 2*np.pi*m/model.length,
                             *(np.sqrt(s[m]/counts[m]) for s in spectra)))
    return rows



def validate_rates(exe, out, scratch):
    mode = (2, 2)
    length, end = 128., 640.
    records, history, surfaces = [], [], {}
    for n, dt in ((64, 8.), (128, 4.)):
        phase = 2*np.pi*mode[0]*(np.arange(n)+.5)/n
        initial = 1.+.001*np.cos(phase[:, None]+phase)
        steps, stride = int(end/dt), int(end/dt)//16
        time = dt*np.array([0, 1, *range(stride, steps+1, stride)])
        for case, attraction in (('attraction', .0084375), ('capillary', 0.)):
            model = WallFilm2D(n, length, dt, attraction)
            jam, basilisk, metrics = run_pair(exe, model, steps, 0, 0, scratch, stride,
                                              initial=initial)
            continuum = float(continuous_rate(model)[mode])
            midpoint = float(np.log(model.gain[mode])/dt)
            for solver, fields in (('JAM', jam), ('Basilisk', basilisk)):
                amplitude = 2*abs(np.fft.fft2(fields-1.)[:, mode[0], mode[1]])/n**2
                slope, intercept = np.polyfit(time, np.log(amplitude/amplitude[0]), 1)
                error = float(abs(slope/continuum-1))
                midpoint_error = float(abs(slope/midpoint-1))
                assert slope*continuum > 0 and error < .01 and midpoint_error < 1e-4
                assert np.max(abs(fields.mean(axis=(1, 2))-1.)) < 1e-12
                records.append(dict(case=case, grid=n, dt=dt, solver=solver, fitted=slope,
                                    continuum=continuum, midpoint=midpoint,
                                    continuum_error=error, midpoint_error=midpoint_error,
                                    fit_max_log_residual=float(np.max(abs(
                                        np.log(amplitude/amplitude[0])-slope*time-intercept))),
                                    max_paired_field_error=metrics['max_paired_error']))
                history.extend((case, n, solver, t, a) for t, a in zip(time, amplitude))
                if n == 64:
                    surfaces[f'{case}_{solver.lower()}'] = fields[[0, -1]]
    with (out/'rate_history.csv').open('w', encoding='utf-8') as file:
        file.write('case,grid,solver,time,amplitude\n')
        for case, n, solver, t, a in history:
            file.write(f'{case},{n},{solver},{t:.16g},{a:.16g}\n')
    np.savez_compressed(out/'rate_surfaces.npz', x=(np.arange(64)+.5)*length/64,
                        time=[0., end], **surfaces)
    (out/'rates.json').write_text(json.dumps(records, indent=2), encoding='utf-8')
    print('Both signed rates passed for both solvers at 64 and 128 cells per side.',
          flush=True)


def trajectory(exe, n, seed, scratch, save_fields):
    film = WallFilm2D(n, 2.*n, 8., theta=1e-5)
    samples = [1, 2, 5, 17]
    times = film.dt*np.array([1, 10, 40, 160])
    power, surfaces, metrics = [], {}, {}
    first = None
    for label, noisy in (('continuous', 160), ('first_step', 1)):
        jam, basilisk, metrics[label] = run_pair(exe, film, 160, noisy, seed, scratch, 10)
        fields = np.stack((jam[samples], basilisk[samples]))
        mean_error = float(np.max(abs(fields.mean(axis=(-1, -2))-1.)))
        assert mean_error < 1e-12 and np.all(np.isfinite(fields))
        metrics[label]['max_mean_error'] = mean_error
        if first is None:
            first = fields[:, 0].copy()
        else:
            np.testing.assert_array_equal(first, fields[:, 0])
        power.append(abs(np.fft.fft2(fields-1.)/n**2)**2)
        if save_fields:
            surfaces[label] = fields
    rows = np.asarray(shell_rows(film, np.asarray(power), times))
    return seed, rows, surfaces, metrics


def discrete_power(film, rows):
    # The midpoint recurrence preserves the capillary equilibrium variance exactly.
    qd2 = 3*film.mobility
    equilibrium = np.divide(film.theta, film.length**2*qd2,
                            out=np.zeros_like(qd2), where=qd2 > 0)
    one_step = film.noise**2/film.n**2
    shell = np.rint(film.q*film.length/(2*np.pi)).astype(int).ravel()
    exact = np.empty(len(rows))
    for protocol, t in np.unique(rows[:, :2], axis=0):
        step = round(t/film.dt)
        if protocol == 0:
            field = equilibrium*(1-film.gain**(2*step))
        else:
            field = one_step*film.gain**(2*(step-1))
        select = (rows[:, 0] == protocol) & (rows[:, 1] == t)
        shells = rows[select, 2].astype(int)
        exact[select] = np.bincount(shell, weights=field.ravel())[shells]/rows[select, 3]
    return exact


def spectra(exe, out, n, replicas, workers, scratch):
    start = perf_counter()
    records = {}
    first_seed = 61100
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for seed, rows, fields, metrics in pool.map(
            trajectory, [exe]*replicas, [n]*replicas, range(first_seed, first_seed+replicas),
            [scratch]*replicas, [True]+[False]*(replicas-1)
        ):
            records[seed] = (rows, metrics)
            if fields:
                np.savez_compressed(out/'tcw_surfaces_seed_61100.npz',
                                    x=(np.arange(n)+.5)*2., time=[8., 80., 320., 1280.],
                                    solvers=['JAM', 'Basilisk'], **fields)
            print(f'{len(records)}/{replicas} seeds, both noise protocols finished', flush=True)

    seeds = list(records)
    all_rows = np.stack([records[seed][0] for seed in seeds])
    rows = all_rows[0, :, :5]
    power = all_rows[:, :, 5:]**2
    mean = power.mean(axis=0)
    sem = power.std(axis=0, ddof=1)/np.sqrt(replicas) if replicas > 1 else np.full_like(mean, np.nan)
    np.savez_compressed(out/'tcw_seed_power.npz', seeds=seeds, shell_rows=rows,
                        power=power, columns=['continuum', 'jam', 'basilisk'])
    film = WallFilm2D(n, 2.*n, 8., theta=1e-5)
    midpoint = discrete_power(film, rows)
    shell_count = int(np.floor(n/6-.5))
    rms_floor = 1e-17
    rms = np.sqrt(mean)
    assert len(rows) == 8*shell_count and np.all(midpoint > 0)
    # The binary64 height stores 1 + disturbance, so decayed modes need an absolute tolerance.
    np.testing.assert_allclose(rms[:, 1], rms[:, 2], rtol=1e-9, atol=rms_floor)
    np.savetxt(out/'tcw.csv', np.column_stack((rows, rms[:, 0], np.sqrt(midpoint),
                                               rms[:, 1:], sem[:, 1:])),
               delimiter=',', comments='',
               header='protocol,time,shell,mode_count,q,rms_continuum,rms_midpoint,'
                      'rms_jam,rms_basilisk,power_sem_jam,power_sem_basilisk')
    plotted = rows[:, 1] > film.dt
    report = dict(grid=n, length=film.length, dt=film.dt, theta=film.theta, steps=160,
                  pairs_per_protocol=replicas, complete_shells=shell_count, workers=workers,
                  elapsed_seconds=perf_counter()-start, first_seed=first_seed,
                  numerical_rms_floor=rms_floor)
    for protocol, label in enumerate(('continuous', 'first_step')):
        select = plotted & (rows[:, 0] == protocol)
        use = select & (midpoint >= rms_floor**2)
        discrete_error = np.sqrt(mean[use, 1]/midpoint[use])-1
        continuum_error = np.sqrt(mean[use, 1]/mean[use, 0])-1
        report[label] = dict(
            max_rms_error_vs_midpoint=float(np.max(abs(discrete_error))),
            max_rms_error_vs_continuum=float(np.max(abs(continuum_error))),
            median_rms_error_vs_midpoint=float(np.median(abs(discrete_error))),
            excluded_tail_rows=int(np.sum(select & ~use)),
            max_paired_rms_absolute_difference=float(np.max(abs(rms[select, 1]-rms[select, 2]))),
            max_paired_field_error=max(m[label]['max_paired_error'] for _, m in records.values()),
            max_mean_error=max(m[label]['max_mean_error'] for _, m in records.values()))
    (out/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--scratch', type=Path)
    parser.add_argument('--grid', type=int, default=512)
    parser.add_argument('--replicas', type=int, default=30)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--rates-only', action='store_true')
    args = parser.parse_args()
    n = args.grid
    if n < 16 or n > 512 or n & (n-1) or min(args.replicas, args.workers) < 1:
        parser.error('grid must be a power of two between 16 and 512, counts must be positive')
    args.output.mkdir(parents=True)
    scratch = args.scratch or args.output
    validate_rates(args.executable.resolve(), args.output, scratch)
    if not args.rates_only:
        spectra(args.executable.resolve(), args.output, n, args.replicas, args.workers, scratch)
