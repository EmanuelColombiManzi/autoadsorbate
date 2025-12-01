from ase.build import fcc111
from autoadsorbate import Fragment, Surface
from ase.io import write
import uuid

###########################
# ----- USER INPUT -------
###########################

slab_metal = 'Cu'
slab_size = [4, 4, 4]
fragment_conformers = 10
cations = ['Na', 'K', 'Li', 'Cs']
out_filename = 'XXgenerated_traj.xyz'

smiles = [
    'ClC(=O)[O-]',
    'Cl[O+]=C([O-])[O-]',
    'ClC(=C=O)[O-]',
    'S1SOC([O-])=C1[O-]',
]

smiles_plus = [
    'ClC(=O)OCl',
    'Cl[O+]=C(O2)OS2',
    'ClC(=C=O)OCl',
    'S1SOC(O2)=C1OS2',
]

###########################
# ----- FUNCTIONS --------
###########################

def generate_uid():
    """Generate a UUID4 string for unique identification."""
    return str(uuid.uuid4())

def make_smiles(surface, smiles_list):
    """
    Place fragments on all surface sites without cations.
    Returns list of ASE Atoms with updated info.
    """
    traj = []
    for smile in smiles_list:
        fragment = Fragment(smile, to_initialize=fragment_conformers)
        populated = surface.get_populated_sites(
            fragment,
            site_index='all',
            sample_rotation=True,
            mode='heuristic',
            conformers_per_site_cap=5,
            overlap_thr=1.5,
            verbose=False,
        )
        for atoms in populated:
            # ensure adsorbate_info exists
            if 'adsorbate_info' not in atoms.info:
                atoms.info['adsorbate_info'] = {}
            atoms.info['adsorbate_info']['cation'] = 'None'
            traj.append(atoms)
    return traj

def make_smiles_plus(surface, smiles_list, cations_list):
    """
    Place fragments on surface and replace dummy atom with each cation.
    Returns list of ASE Atoms with updated info.
    """
    traj = []
    fragment_pairs = []

    # First, generate all populated structures
    for smile in smiles_list:
        fragment = Fragment(smile, to_initialize=fragment_conformers)
        populated = surface.get_populated_sites(
            fragment,
            site_index='all',
            sample_rotation=True,
            mode='heuristic',
            conformers_per_site_cap=5,
            overlap_thr=1.5,
            verbose=False,
        )
        for atoms in populated:
            fragment_pairs.append((atoms, fragment))

    # Now assign cations and compute adsorbate formula
    out_traj = []
    for atoms, fragment in fragment_pairs:
        for cation in cations_list:
            a = atoms.copy()
            cf = fragment.get_conformer(0)
            # compute formula using CHO atoms only
            cho_indices = [atom.index for atom in cf if atom.symbol in ['C', 'H', 'O']]
            adsorbate_formula = cf[cho_indices].get_chemical_formula(empirical=True)

            if 'adsorbate_info' not in a.info:
                a.info['adsorbate_info'] = {}

            a.info['adsorbate_info']['adsorbate_formula'] = adsorbate_formula
            a.info['adsorbate_info']['cation'] = cation

            # assume last atom is dummy cation placeholder
            a[-1].symbol = cation

            out_traj.append(a)

    return out_traj

def main():
    # build slab and surface
    slab = fcc111(slab_metal, size=slab_size, vacuum=10)
    surface = Surface(slab)
    surface.sym_reduce()

    # generate trajectories
    out_traj = []
    out_traj += make_smiles(surface, smiles)
    out_traj += make_smiles_plus(surface, smiles_plus, cations)

    # assign unique IDs
    for a in out_traj:
        a.info['uid'] = generate_uid()

    print(f'Writing... {len(out_traj)} structures to {out_filename}')
    write(out_filename, out_traj)

###########################
# ----- ENTRY POINT ------
###########################

if __name__ == '__main__':
    main()
