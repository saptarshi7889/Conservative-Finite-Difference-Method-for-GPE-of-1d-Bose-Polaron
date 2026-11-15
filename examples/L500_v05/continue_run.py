"""Continue t=10 to t=20 in the same folder; plot at t=15."""
from concurrent.futures import ThreadPoolExecutor
import fcntl
import lzma
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from run import compress, POLL_SECONDS, COMPRESSION_WORKERS

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT/'results_L500_v05'
EXAMPLE = Path(__file__).resolve().parent
START = 10.0
END = 20.0
PLOT_STEP = 1500
START_STEP = 1000
EXPECTED_ROWS = 201
EXPECTED_SNAPSHOTS = 41


def main():
    with (ROOT/'build/L500_v05.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        history = (RESULTS/'velocity.dat').read_bytes()
        rows = [list(map(float, line.split())) for line in history.decode().splitlines() if line and not line.startswith('#')]
        if rows[-1][0] != START:
            raise SystemExit('Existing velocity history must end at t=10.')
        if (RESULTS/'parameters_continuation_1000.txt').exists():
            raise SystemExit('This continuation has already been started; refusing to duplicate output.')
        for name, expected in {'length':500, 'dx':0.001, 'dt':0.01, 'mass_ratio':0.333, 'coupling':1, 'gamma':0.1, 'velocity_0':0.5, 'iterations':10, 'tolerance':1e-8, 'save_every':50, 'velocity_every':10}.items():
            parameters = dict(line.split() for line in (RESULTS/'parameters.txt').read_text().splitlines())
            if float(parameters[name]) != expected:
                raise SystemExit(f'Original parameter mismatch: {name}')
        temporary = tempfile.NamedTemporaryFile(dir=ROOT/'build', prefix='restart_t10_', suffix='.dat', delete=False)
        restart = Path(temporary.name)
        try:
            with temporary, lzma.open(RESULTS/'phi_1000.dat.xz', 'rb') as stream:
                shutil.copyfileobj(stream, temporary)
            command = [str(ROOT/'build/gpe_research_continuation'), '--length','500','--dx','0.001',
                       '--dt','0.01','--mass-ratio','0.333','--coupling','1','--gamma','0.1',
                       '--velocity','0.5','--iterations','10','--tolerance','1e-8',
                       '--save-every','50','--velocity-every','10','--start',str(START),'--end',str(END),
                       '--restart',str(restart),'--restart-velocity',repr(rows[-1][1]),'--output',str(RESULTS)]
            pending = {}
            plotted = False
            plot_attempted = False
            errors = []
            env = dict(os.environ, MPLCONFIGDIR='/tmp/gpe-L500-mpl')
            with (EXAMPLE/'continuation.log').open('w') as log, ThreadPoolExecutor(max_workers=COMPRESSION_WORKERS) as worker:
                process = subprocess.Popen(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
                try:
                    while True:
                        if (RESULTS/'parameters_continuation_1000.txt').exists():
                            # The solver has loaded and closed the temporary checkpoint.
                            restart.unlink(missing_ok=True)
                        for path in RESULTS.glob('phi_*.dat'):
                            if int(path.stem.split('_')[1]) > START_STEP and path not in pending:
                                pending[path] = worker.submit(compress, path)
                        checkpoint = RESULTS/f'phi_{PLOT_STEP}.dat'
                        if not plot_attempted and (checkpoint.exists() or Path(str(checkpoint)+'.xz').exists()):
                            plot_attempted = True
                            result = subprocess.run([sys.executable, str(EXAMPLE/'plot_t15.py')], env=env,
                                                    capture_output=True, text=True)
                            (EXAMPLE/'plot_t15.log').write_text(result.stdout+result.stderr)
                            plotted = result.returncode == 0
                            print('t=15 plots generated.' if plotted else 'Plotting failed; solver continues. See plot_t15.log.', flush=True)
                        try:
                            code = process.wait(timeout=POLL_SECONDS)
                            break
                        except subprocess.TimeoutExpired:
                            pass
                    for path in RESULTS.glob('phi_*.dat'):
                        if int(path.stem.split('_')[1]) > START_STEP and path not in pending:
                            pending[path] = worker.submit(compress, path)
                finally:
                    if process.poll() is None:
                        process.terminate()
                        process.wait()
                for path, future in pending.items():
                    try:
                        future.result()
                    except Exception as error:
                        errors.append(f'{path.name}: {error}')
                        print('Compression failed; raw file retained:', path, error, flush=True)
            if code:
                raise SystemExit(f'Solver failed (exit {code}); see continuation.log. Partial data preserved.')
            if not plotted:
                subprocess.run([sys.executable, str(EXAMPLE/'plot_t15.py')], env=env, check=True)
            data = (RESULTS/'velocity.dat').read_bytes()
            assert data.startswith(history), 'Previous velocity rows changed!'
            final_rows = [list(map(float,line.split())) for line in data.decode().splitlines() if line and not line.startswith('#')]
            assert len(final_rows) == EXPECTED_ROWS and final_rows[-1][0] == END
            assert all(a[0] < b[0] for a,b in zip(final_rows,final_rows[1:]))
            assert len(list(RESULTS.glob('phi_*.dat.xz'))) == EXPECTED_SNAPSHOTS
            if errors:
                raise SystemExit('Compression errors: '+'; '.join(errors))
            assert not list(RESULTS.glob('phi_*.dat'))
            print(f'Completed t=20; original velocity history preserved; {len(final_rows)} velocity rows; {EXPECTED_SNAPSHOTS} compressed configurations.', flush=True)
            (EXAMPLE/'continuation_summary.txt').write_text(f'Completed t=20. Final velocity: {final_rows[-1][1]:.17g}\nOriginal 0<=t<=10 velocity history preserved byte-for-byte.\n201 velocity rows; 41 compressed configurations; no raw configurations remain.\nPlots at t=15: velocity_t15.png, density_t15.png, density_t15_full_system.png\n')
        finally:
            restart.unlink(missing_ok=True)


if __name__ == '__main__':
    main()
