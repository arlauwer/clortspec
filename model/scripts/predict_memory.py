"""
predict_memory.py
=================
Estimate the size of the XRayIonicGasMix tables of combined-grid SKIRT runs
from their shell.txt, to find the runs that need the most memory.

  python3 scripts/predict_memory.py RESOURCE_DIR run_00134 run_00135 ...
  python3 scripts/predict_memory.py RESOURCE_DIR --min-fraction 1e-6 run_*   # preview a tolerance
  python3 scripts/predict_memory.py RESOURCE_DIR --all-ions                  # worst case

RESOURCE_DIR holds Ionic_PA.txt, Ionic_FL.txt and Ionic_LN.txt
(SKIRT9_Resources_AtomsMolecules/XRayIonicGas).

SKIRT builds one XRayIonicGasMix per distinct (ion fractions, T) cell. Only
ions with a non-zero fraction are kept, and each mix stores a table of
numLambda x (ions + fluorescence + 2 x lines) doubles (_cumsigmascavv; every
line is a resonant and a recombination channel), where numLambda = 2500
points/dex over the energy range + 7 points around every photo-absorption
threshold + every fluorescence and recombination energy. Both factors grow
with the ions present, so memory grows faster than the number of ions.

--min-fraction applies steps/skirt_generator.prune_trace_ions to the existing
shell.txt (default: none, i.e. the file as written), to preview the effect of
config.SHELL_ION_MIN_FRACTION without regenerating. The printed GB is that
table summed over the mixes. Measured SKIRT peaks are ~0.8x this estimate
(four radial-test runs at logNHs 24: 17.1 GB vs 21.2-21.7 predicted, 33.5-33.8
vs 41.9-42.6), so the summary line reports 0.8x the largest estimate as the
expected peak, and with --ram GB how many simultaneous runs fit in 90% of it.
"""

import argparse
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config as C

E_MIN_EV, E_MAX_EV = C.E_MIN_KEV * 1e3, C.E_MAX_KEV * 1e3
PER_DEX = 2500
PEAK_FACTOR = 0.8  # measured SKIRT peak / predicted mix tables (radial test)
HC_EV_M = 1.239841984e-6  # h c / e: photon energy (eV) <-> wavelength (m)


def ion_key(label):
	"""'Fe+16' -> (26, 10): (Z, number of electrons)."""
	symbol, charge = label.split("+")
	z = ELEMENTS.index(symbol) + 1
	return z, z - int(charge)


ELEMENTS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn").split()


def read_resource(path, energy_col):
	"""Per (Z, N): number of records and the energies (eV) in the simulated range."""
	count, energies = Counter(), {}
	for line in open(path):
		if line.startswith("#") or not line.strip():
			continue
		cols = line.split()
		key = (int(float(cols[0])), int(float(cols[1])))
		count[key] += 1
		if energy_col is not None:
			energies.setdefault(key, []).append(float(cols[energy_col]))
	return count, energies


def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("resources", type=Path, help="directory with Ionic_PA/FL/LN.txt")
	parser.add_argument("runs", nargs="*", type=Path, help="run directories holding a shell.txt")
	parser.add_argument("--min-fraction", type=float, default=0.0,
						help="drop ions below this fraction of their element first (default 0: as written)")
	parser.add_argument("--all-ions", action="store_true", help="worst case: every ion in every cell")
	parser.add_argument("--ram", type=float, help="server RAM in GB: report how many runs fit at once")
	args = parser.parse_args()
	res = args.resources
	_, pa_e = read_resource(res / "Ionic_PA.txt", 4)
	fl_n, fl_e = read_resource(res / "Ionic_FL.txt", 5)
	# lines: resonant channel + a recombination "fluorescence" channel (if its yield is
	# non-zero at T, assumed here) whose energy is also added to the wavelength grid
	ln_n, ln_lam = read_resource(res / "Ionic_LN.txt", 5)
	fine = int(np.ceil(PER_DEX * np.log10(E_MAX_EV / E_MIN_EV))) + 1

	def mix_size(present):
		n_lambda = fine
		n_int = len(present)  # electron scattering: one channel per ion
		for key in present:
			# threshold points: 7 per threshold, roughly all within +-2 thermal widths of E_th
			n_lambda += 7 * sum(E_MIN_EV < e < E_MAX_EV for e in pa_e.get(key, []))
			n_lambda += sum(E_MIN_EV < e < E_MAX_EV for e in fl_e.get(key, []))
			n_lambda += sum(E_MIN_EV < HC_EV_M / lam < E_MAX_EV for lam in ln_lam.get(key, []))
			n_int += fl_n[key] + 2 * ln_n[key] + (key[1] == 3)  # + Li-like -> He-like z channel
		return n_lambda, n_int

	from steps.skirt_generator import prune_trace_ions, template_ions
	ions = template_ions()
	keys = [ion_key(i) for i in ions]

	if args.all_ions:
		n_lambda, n_int = mix_size(keys)
		gb = C.SHELL_N_RADIAL * 8.0 * n_lambda * n_int / 1e9
		print(f"all {len(keys)} ions: numLambda {n_lambda}, {n_int} channels per mix; "
			  f"x {C.SHELL_N_RADIAL} cells = {gb:.1f} GB of mix tables")
		return

	largest = (0.0, None)
	for run in args.runs:
		data = np.loadtxt(run / "shell.txt", comments="#", ndmin=2)
		fractions = prune_trace_ions(dict(zip(ions, data[:, 8:].T)), args.min_fraction)
		present = np.array([fractions[ion] for ion in ions]).T > 0
		cells = {(t, *row) for t, row in zip(data[:, 7], present)}  # the two walls share a mix
		total = sum(8.0 * np.prod(mix_size([k for k, p in zip(keys, cell[1:]) if p])) for cell in cells)
		mean_ions = np.mean([sum(cell[1:]) for cell in cells])
		print(f"{run}: {len(cells)} mixes, {mean_ions:.0f} ions/mix, {total / 1e9:.1f} GB of mix tables")
		largest = max(largest, (total / 1e9, str(run)))

	if len(args.runs) > 1:
		peak = PEAK_FACTOR * largest[0]
		print(f"largest: {largest[1]}, {largest[0]:.1f} GB of mix tables -> expected SKIRT peak ~{peak:.0f} GB")
	if args.ram and args.runs:
		peak = PEAK_FACTOR * largest[0]
		print(f"with {args.ram:.0f} GB RAM: up to {max(1, int(0.9 * args.ram // peak))} simultaneous runs (-j), "
			  f"assuming they all need the largest peak")


if __name__ == "__main__":
	main()
