from ase.build import fcc111
from autoadsorbate import Fragment, Surface
from ase.io import write
import uuid
import base64

def generate_uid(length=5):
    """Generate a short UID from UUID4 using base64 encoding."""
    uid = str(uuid.uuid4())
    # Convert to base64 and take first N characters
    encoded = base64.urlsafe_b64encode(uid.encode()).decode()
    return encoded[:length]

def make_smiles_plus(surface, smiles_plus):
    
    traj = []
    for smile in smiles_plus:

        fragment = Fragment(smile, to_initialize=fragment_conformers)

        traj += surface.get_populated_sites(
            fragment,
            site_index='all',
            sample_rotation=True, mode='heuristic',
            conformers_per_site_cap=5, overlap_thr=1.5, verbose=False
            )

    out_traj = []
    for cation in cations:
        for atoms in traj:
            a = atoms.copy()
            cf = fragment.get_conformer(0)
            adsorbate_formula = cf[
                [atom.index for atom in cf if atom.symbol in ['C', 'H', 'O']]
                ].get_chemical_formula(empirical=True)
            a.info['adsorbate_info']['adsorbate_formula'] = adsorbate_formula
            a.info['adsorbate_info']['cation'] = cation
            a[-1].symbol = cation
            out_traj+=[a]
    return out_traj
    

def make_smiles(surface, smiles):
    
    traj = []
    for smile in smiles:

        fragment = Fragment(smile, to_initialize=fragment_conformers)

        traj += surface.get_populated_sites(
            fragment,
            site_index='all',
            sample_rotation=True, mode='heuristic',
            conformers_per_site_cap=5, overlap_thr=1.5, verbose=False
            )

    for atoms in traj:
        atoms.info['adsorbate_info']['cation'] = 'None'
        
    return traj


def main():
    slab = fcc111(slab_metal, size=slab_size, vacuum=10)
    surface = Surface(slab)
    surface.sym_reduce()

    out_traj = []
    out_traj += make_smiles(surface, smiles)
    out_traj += make_smiles_plus(surface, smiles_plus)

    for a in out_traj:
        a.info['uid'] = generate_uid()

    print(f'Writing. . . {len(out_traj) = }, {out_filename = }')
    write(out_filename, out_traj)


###################### user input

slab_metal = 'Cu'
slab_size = [4,4,4]
fragment_conformers = 10
cations = ['Na','K','Li','Cs']
out_filename = 'generated_traj.xyz'

smiles = [
    'ClC(=O)[O-]',
    'Cl[O+]=C([O-])[O-]',
    'ClC(=C=O)[O-]',
    'S1SOC([O-])=C1[O-]',
]

smiles_plus = [ #using a similar hack as for the surrogate atom.
    'ClC(=O)OCl',
    'Cl[O+]=C(O2)OS2',
    'ClC(=C=O)OCl',
    'S1SOC(O2)=C1OS2',
]

if __name__ == '__main__':
    main()