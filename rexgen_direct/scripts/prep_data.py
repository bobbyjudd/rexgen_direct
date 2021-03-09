import rdkit.Chem as Chem
import numpy as np
from tqdm import tqdm

'''
This script prepares the data used in Wengong Jin's NIPS paper on predicting reaction outcomes for the modified
forward prediction script. Rather than just training to predict which bonds change, we make a direct prediction
on HOW those bonds change
'''
bond_fdim = 6
binary_fdim = 4 + bond_fdim

def get_changed_bonds(rxn_smi):
    reactants = Chem.MolFromSmiles(rxn_smi.split('>')[0])
    products  = Chem.MolFromSmiles(rxn_smi.split('>')[2])

    conserved_maps = [a.GetProp('molAtomMapNumber') for a in products.GetAtoms() if a.HasProp('molAtomMapNumber')]
    bond_changes = set() # keep track of bond changes

    # Look at changed bonds
    bonds_prev = {}
    for bond in reactants.GetBonds():
        nums = sorted(
            [bond.GetBeginAtom().GetProp('molAtomMapNumber'),
             bond.GetEndAtom().GetProp('molAtomMapNumber')])
        if (nums[0] not in conserved_maps) and (nums[1] not in conserved_maps): continue
        bonds_prev['{}~{}'.format(nums[0], nums[1])] = bond.GetBondTypeAsDouble()
    bonds_new = {}
    for bond in products.GetBonds():
        nums = sorted(
            [bond.GetBeginAtom().GetProp('molAtomMapNumber'),
             bond.GetEndAtom().GetProp('molAtomMapNumber')])
        bonds_new['{}~{}'.format(nums[0], nums[1])] = bond.GetBondTypeAsDouble()


    for bond in bonds_prev:
        if bond not in bonds_new:
            bond_changes.add((bond.split('~')[0], bond.split('~')[1], 0.0)) # lost bond
        else:
            if bonds_prev[bond] != bonds_new[bond]:
                bond_changes.add((bond.split('~')[0], bond.split('~')[1], bonds_new[bond])) # changed bond
    for bond in bonds_new:
        if bond not in bonds_prev:
            bond_changes.add((bond.split('~')[0], bond.split('~')[1], bonds_new[bond]))  # new bond

    return bond_changes
    
def process_patent_list(fpath):
    with open(fpath, 'r') as fid_in, open(fpath + '.proc', 'w') as fid_out, open(fpath + '.error', 'w') as err_out:
        for line in tqdm(fid_in):
            if line.startswith('ReactionSmiles'):
                continue
            rxn_smi = line.strip().split()[0]
            try:
                bond_changes = get_changed_bonds(rxn_smi)
            except:
                err_out.write(rxn_smi)
                continue
            if bond_changes:
                fid_out.write('{} {}\n'.format(rxn_smi, ';'.join(['{}-{}-{}'.format(x[0], x[1], x[2]) for x in bond_changes])))
            else:
                err_out.write(rxn_smi)

def get_bin_feature(r, max_natoms):
    comp = {}
    for i, s in enumerate(r.split('.')):
        mol = Chem.MolFromSmiles(s)
        for atom in mol.GetAtoms():
            comp[atom.GetIntProp('molAtomMapNumber') - 1] = i
    n_comp = len(r.split('.'))
    rmol = Chem.MolFromSmiles(r)
    n_atoms = rmol.GetNumAtoms()
    bond_map = {}
    for bond in rmol.GetBonds():
        a1 = bond.GetBeginAtom().GetIntProp('molAtomMapNumber') - 1
        a2 = bond.GetEndAtom().GetIntProp('molAtomMapNumber') - 1
        bond_map[(a1,a2)] = bond_map[(a2,a1)] = bond
        
    features = []
    for i in range(max_natoms):
        for j in range(max_natoms):
            f = np.zeros((binary_fdim,))
            if i >= n_atoms or j >= n_atoms or i == j:
                features.append(f)
                continue
            if (i,j) in bond_map:
                bond = bond_map[(i,j)]
                f[1:1+bond_fdim] = bond_features(bond)
            else:
                f[0] = 1.0
            f[-4] = 1.0 if comp[i] != comp[j] else 0.0
            f[-3] = 1.0 if comp[i] == comp[j] else 0.0
            f[-2] = 1.0 if n_comp == 1 else 0.0
            f[-1] = 1.0 if n_comp > 1 else 0.0
            features.append(f)

def bond_features(bond):
    bt = bond.GetBondType()
    return np.array([bt == Chem.rdchem.BondType.SINGLE, bt == Chem.rdchem.BondType.DOUBLE, bt == Chem.rdchem.BondType.TRIPLE, bt == Chem.rdchem.BondType.AROMATIC, bond.GetIsConjugated(), bond.IsInRing()], dtype=np.float32)

def graph_encoding_capable(path, idxfunc=lambda x:x.GetIdx()):
    data = []
    max_natoms = 0
    with open(path, 'r') as f:
        for line in f.readlines():
            r,e = line.split()
            
            max_natoms = max(max_natoms, Chem.MolFromSmiles(r.split('>')[0]).GetNumAtoms())
            data.append((r,e))
    with open('../data/custom_filtered.rsmi.proc', 'w') as f, open('../data/custom_filtered.rsmi.error', 'w') as err:
        for r,e in data:
            try:
                react,_,p = r.split('>')
                mol = Chem.MolFromSmiles(react)
                if not mol:
                    raise ValueError("Could not parse smiles string:", s)
                
                # Test mol_graph
                n_atoms = mol.GetNumAtoms()
                n_bonds = max(mol.GetNumBonds(), 1)

                max_natoms = max(n_atoms, max_natoms)
                for atom in mol.GetAtoms():
                    idx = idxfunc(atom)
                    if idx >= n_atoms:
                        raise Exception(smiles)

                for bond in mol.GetBonds():
                    a1 = idxfunc(bond.GetBeginAtom())
                    a2 = idxfunc(bond.GetEndAtom())
                    idx = bond.GetIdx()
                    bond_features(bond)
                
                # Test get_bin_feature
                get_bin_feature(react,max_natoms)

                f.write('{} {}\n'.format(r, e))
            except:
                err.write('{} {}\n'.format(r, e))

def remove_duplicates(fpath):
    with open(fpath) as f1, open('../data/train.txt.proc') as f2, open('../data/test.txt.proc') as f3, open('../data/valid.txt.proc') as f4, open('../data/custom_filtered_unique.rsmi.proc', 'w') as fout:
        custom_set = set(f1.readlines())
        train_set = set(f2.readlines())
        test_set = set(f3.readlines())
        valid_set = set(f4.readlines())

        full_orig = train_set | test_set | valid_set
        unique_custom = custom_set - full_orig

        fout.write(''.join(unique_custom))

def filter_custom_data(fpath):
    # Format custom dataset
    process_patent_list(fpath)

    # Filter datapoints that do not meet the graph nn input format
    graph_encoding_capable(fpath+'.proc', idxfunc=lambda x:x.GetIntProp('molAtomMapNumber') - 1)

    # Remove any datapoints that existed in the included dataset
    remove_duplicates('../data/custom_filtered.rsmi.proc')

def process_file(fpath):
    with open(fpath, 'r') as fid_in, open(fpath + '.proc', 'w') as fid_out:
        for line in tqdm(fid_in):
            rxn_smi = line.strip().split(' ')[0]
            bond_changes = get_changed_bonds(rxn_smi)
            fid_out.write('{} {}\n'.format(rxn_smi, ';'.join(['{}-{}-{}'.format(x[0], x[1], x[2]) for x in bond_changes])))
    print('Finished processing {}'.format(fpath))

if __name__ == '__main__':
    # Test summarization
    """
    for rxn_smi in [
            '[CH2:15]([CH:16]([CH3:17])[CH3:18])[Mg+:19].[CH2:20]1[O:21][CH2:22][CH2:23][CH2:24]1.[Cl-:14].[OH:1][c:2]1[n:3][cH:4][c:5]([C:6](=[O:7])[N:8]([O:9][CH3:10])[CH3:11])[cH:12][cH:13]1>>[OH:1][c:2]1[n:3][cH:4][c:5]([C:6](=[O:7])[CH2:15][CH:16]([CH3:17])[CH3:18])[cH:12][cH:13]1',
            '[CH3:14][NH2:15].[N+:1](=[O:2])([O-:3])[c:4]1[cH:5][c:6]([C:7](=[O:8])[OH:9])[cH:10][cH:11][c:12]1[Cl:13].[OH2:16]>>[N+:1](=[O:2])([O-:3])[c:4]1[cH:5][c:6]([C:7](=[O:8])[OH:9])[cH:10][cH:11][c:12]1[NH:15][CH3:14]',
            '[CH2:1]([CH3:2])[n:3]1[cH:4][c:5]([C:22](=[O:23])[OH:24])[c:6](=[O:21])[c:7]2[cH:8][c:9]([F:20])[c:10](-[c:13]3[cH:14][cH:15][c:16]([NH2:19])[cH:17][cH:18]3)[cH:11][c:12]12.[CH:25](=[O:26])[OH:27]>>[CH2:1]([CH3:2])[n:3]1[cH:4][c:5]([C:22](=[O:23])[OH:24])[c:6](=[O:21])[c:7]2[cH:8][c:9]([F:20])[c:10](-[c:13]3[cH:14][cH:15][c:16]([NH:19][CH:25]=[O:26])[cH:17][cH:18]3)[cH:11][c:12]12',
            ]:
        print(rxn_smi)
        print(get_changed_bonds(rxn_smi))
    """
    # Process files
    
    #process_file('../data/train.txt')
    #process_file('../data/valid.txt')
    #process_file('../data/test.txt')
    #process_file('../data/test_human.txt')
    filter_custom_data('../data/1976_Sep2016_USPTOgrants_smiles.rsmi')