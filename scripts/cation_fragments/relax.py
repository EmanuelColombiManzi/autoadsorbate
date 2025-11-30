import time
import os
from glob import glob
from joblib import Parallel, delayed

import pandas as pd
import numpy as np

import ase
from ase.optimize import BFGS
from ase import units
from ase.constraints import FixAtoms

from mace.calculators import mace_mp
from pathlib import Path

import sys
from autoadsorbate.utils import freeze_atoms, save_relax_config, relaxatoms

rcut=40
steps = 500
fmax = 0.05
disp=True
prefix = 'MACE_prerelax'
freeze_bottom=False
traj = ase.io.read('./generated_traj.xyz', index = ':')


def main():

    macemp = mace_mp(model='medium',dispersion=disp,device="cpu",dispersion_cutoff=rcut* units.Bohr, default_dtype="float64")

    for atoms in traj:
        c = FixAtoms(mask=[atoms.arrays['fragments'][atom.index] == 0  for atom in atoms])
        atoms.set_constraint(c)

   
    for atoms in traj:
        relaxatoms(atoms, macemp, prefix, steps=steps, freeze_bottom=freeze_bottom, fmax=fmax)
        

if __name__ == '__main__':
    main()
