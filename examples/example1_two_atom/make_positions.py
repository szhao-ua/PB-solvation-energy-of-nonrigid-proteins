#!/usr/bin/env python3
"""
Regenerate the structure pairs of Example 1 in data/solvent/ and data/vacuum/.

The solvent structure holds two unit charges of opposite sign at x = -2.2 and
+2.2 A, radius 2 A. For every seed and Cp, it is written to
data/solvent/2ATO_s<seed>_cp<Cp>.xyzqr (and .pqr), and its perturbed copy to the
file of the same name in data/vacuum/. Each atom is displaced by

    r_hat_{j,i} = r_{j,i} + Cp * (gamma_{j,i} - 1/2),    gamma_{j,i} ~ U(0,1),

with the gamma drawn, as in the DIPB/REG code, by gfortran's RANDOM_SEED(PUT=[s,...,s])
followed by RANDOM_NUMBER, three draws per atom. The coordinates are rounded to three
decimals, as the DIPB/REG code writes them. Because the stream is gfortran's own generator,
this script calls libgfortran through ctypes; the structures it produces are already
committed, so running it is only needed to verify or extend them.

Usage:
    python make_positions.py --libgfortran /path/to/libgfortran.5.dylib   (or .so.5)
    python make_positions.py --libgfortran ... --check     compare with data/, write nothing
"""
import argparse
import ctypes
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
# solvent structure: x, y, z, charge, radius
ATOMS = np.array([[-2.2, 0.0, 0.0, -1.0, 2.0],
                  [2.2, 0.0, 0.0, 1.0, 2.0]])
SEEDS = [12345, 2026, 99999, 123456789, 987654321]
CPS = [0.10, 0.20, 0.30, 0.40, 0.50]


# gfortran array descriptor for a rank-1 INTEGER(4) array (GCC >= 8 layout)
class _Dtype(ctypes.Structure):
    _fields_ = [("elem_len", ctypes.c_size_t), ("version", ctypes.c_int), ("rank", ctypes.c_byte),
                ("type", ctypes.c_byte), ("attribute", ctypes.c_short)]


class _Dim(ctypes.Structure):
    _fields_ = [("stride", ctypes.c_ssize_t), ("lbound", ctypes.c_ssize_t), ("ubound", ctypes.c_ssize_t)]


class _Desc(ctypes.Structure):
    _fields_ = [("base_addr", ctypes.c_void_p), ("offset", ctypes.c_ssize_t), ("dtype", _Dtype),
                ("span", ctypes.c_ssize_t), ("dim", _Dim * 1)]


class GfortranRandom:
    def __init__(self, libpath):
        self.lib = ctypes.CDLL(libpath)

    def put_seed(self, s):
        """CALL RANDOM_SEED(SIZE=n); CALL RANDOM_SEED(PUT=[s, s, ..., s])"""
        n = ctypes.c_int(0)
        self.lib._gfortran_random_seed_i4(ctypes.byref(n), None, None)
        self._seed = (ctypes.c_int * n.value)(*([s] * n.value))    # kept alive for the call
        desc = _Desc(ctypes.cast(self._seed, ctypes.c_void_p), -1, _Dtype(4, 0, 1, 1, 0), 4,
                     (_Dim * 1)(_Dim(1, 1, n.value)))
        self.lib._gfortran_random_seed_i4(None, ctypes.byref(desc), None)

    def draw(self):
        """CALL RANDOM_NUMBER(x) for a REAL(8) scalar"""
        x = ctypes.c_double()
        self.lib._gfortran_random_r8(ctypes.byref(x))
        return x.value


def perturb(rng, atoms, seed, cp):
    rng.put_seed(seed)
    out = atoms.copy()
    for j in range(len(atoms)):
        for i in range(3):
            out[j, i] = round(atoms[j, i] - 0.5 * cp + rng.draw() * cp, 3)
    return out


def write_xyzqr(path, a):
    with open(path, "w") as fh:
        for x, y, z, q, r in a:
            fh.write(f"{x:10.3f}{y:10.3f}{z:10.3f}{q:10.3f}{r:10.3f}\n")


def write_pqr(path, a):
    with open(path, "w") as fh:
        for i, (x, y, z, q, r) in enumerate(a, 1):
            fh.write(f"ATOM  {i:5d}  X   ION X{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f} {q:7.4f} {r:7.4f}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--libgfortran", required=True, help="path to libgfortran (the one the DIPB/REG code was built with)")
    ap.add_argument("--check", action="store_true", help="compare with the committed structures instead of writing")
    args = ap.parse_args()

    rng = GfortranRandom(args.libgfortran)
    worst = 0.0
    for seed in SEEDS:
        for cp in CPS:
            name = f"2ATO_s{seed}_cp{cp:.2f}"
            pair = {"solvent": ATOMS, "vacuum": perturb(rng, ATOMS, seed, cp)}
            for state, a in pair.items():
                stem = os.path.join(DATA, state, name)
                if args.check:
                    worst = max(worst, np.abs(np.loadtxt(stem + ".xyzqr", ndmin=2) - a).max())
                else:
                    os.makedirs(os.path.dirname(stem), exist_ok=True)
                    write_xyzqr(stem + ".xyzqr", a)
                    write_pqr(stem + ".pqr", a)
            print(("checked " if args.check else "wrote ") + f"solvent/{name} and vacuum/{name}")
    if args.check:
        print(f"largest difference from the committed structures: {worst:.1e} A")


if __name__ == "__main__":
    main()
