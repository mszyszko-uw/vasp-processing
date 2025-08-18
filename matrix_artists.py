import numpy as np

def MatrixSculptor(mat: np.ndarray, 
                   projection = {'direction':[0], 'ion':..., 'orbital':...}
                   ):
    try: mat = np.sum(mat[..., projection['orbital']], axis=-1)
    except IndexError: mat = np.sum(mat, axis=-1)
    try: mat = np.sum(mat[..., projection['ion']], axis=-1)
    except IndexError: mat = np.sum(mat, axis=-1)
    try: mat = np.sum(mat[..., projection['direction']], axis=-1)
    except IndexError: mat = np.sum(mat, axis=-1)
    return mat


def ProjectionPainter(ions:list, num_ions:list, description:str):
    orbital_names = {'s':[0], 'p':[1,2,3], 'd':[4,5,6,7,8], 'f':[i for i in range(9,16)],
                     'sp':[0,1,2,3,4], 'pd':[i for i in range(1,9)], 'spd':[i for i in range(0,9)],
                     'py':[1], 'pz':[2], 'px':[3],
                     'dxy':[4], 'dyz':[5], 'dz2':[6], 'dxz':[7], 'x2-y2':[8],
                     'fy3x2':[9], 'fxyz':[10], 'fyz2':[11], 'fz3':[12], 'fxz2':[13], 'fzx2':[14], 'fx3':[15]}
        
    direction_names = {'up':[0], 'down':[1], 
                    'tot':[0], 'x':[1], 'y':[2], 'z':[3]}

    # ions = vaspout['input']['poscar']['ion_types'].asstr()[:]
    # num_ions = vaspout['input']['poscar']['number_ion_types'][:]

    now = 0
    list_ions = []
    for n in num_ions:
        list_ions.append([i for i in range(now, now+n)])
        now += n

    atom_indexes = {}
    for name,pos in zip(ions,list_ions): atom_indexes[name] = pos

    projection = {'direction':[], 'ion':[], 'orbital':[]}
    if description == 'clean': return projection

    description = description.split()

    iatoms = []
    for i in description:
        if i.isnumeric(): iatoms.append(i)

    for atom in atom_indexes.keys():
        if atom in description:
            if len(iatoms):
                for iatom in iatoms:
                    projection['ion'].append(atom_indexes[atom][int(iatom)-1])
            else: projection['ion'] = atom_indexes[atom]
    if projection['ion']==[]: projection['ion'] = ...

    for dir in direction_names.keys():
        if dir in description: projection['direction'] += direction_names[dir]
    if projection['direction']==[]: projection['direction'] = [0]

    for orb in orbital_names.keys():
        if orb in description: projection['orbital'] += orbital_names[orb]
    if projection['orbital']==[]: projection['orbital'] = ...
    return projection


def SegmentArranger(nSeg:int, mat: np.ndarray, arrangement: str = ''):
    if arrangement == '': return mat
    # nSeg = vaspout['input']['kpoints']['number_line_segments'][()]/2
    mat = np.split(mat, nSeg)
    arrangement = arrangement.split()

    newmat = []
    for i in arrangement:
        i = int(i)
        if i > 0: 
            try: newmat = np.concatenate([newmat,mat[i-1]])
            except ValueError: newmat = mat[i-1]
        if i < 0: 
            try: newmat = np.concatenate([newmat,mat[np.abs(i)-1][::-1,...]])
            except ValueError: newmat = mat[np.abs(i)-1][::-1,...]

    return newmat


def KlineEngineer(Lvec, Scal, Kpts, Klen, Klab, 
                  nSeg, 
                  arrangement: str = ''
                  ):
    # Lvec = vaspout['input']['poscar']['lattice_vectors'][:]
    # Scal = vaspout['input']['poscar']['scale'][()]
    # Kpts = vaspout['results']['electron_eigenvalues']['kpoint_coords'][:]
    # Klen = vaspout['input']['kpoints']['number_kpoints'][()]
    # Klab = vaspout['input']['kpoints']['labels_kpoints'].asstr()[:]

    Kpts = SegmentArranger(nSeg, Kpts, arrangement)
    Dvec = Kpts.copy()
    a1,a2,a3 = Lvec*Scal
    bM = np.zeros((3,3))
    V = np.dot(a1, np.cross(a2, a3))
    bM[0,...] = 2*np.pi*np.cross(a2, a3) / V
    bM[1,...] = 2*np.pi*np.cross(a3, a1) / V
    bM[2,...] = 2*np.pi*np.cross(a1, a2) / V
    Rval = np.sum(bM**2, axis = 1)**.5
    Dvec *= Rval
    Dvec[1:] -= Dvec[:-1]
    Dvec = np.cumsum(np.sum((Dvec)**2, axis = 1)**.5)

    Nkpt = len(Dvec)
    Ktik = []
    for i in range(0,Nkpt//Klen):
        Ktik.append(i*Klen)
        Ktik.append((i+1)*Klen-1)
    Ktik = Dvec[Ktik]
    Klab = SegmentArranger(nSeg, Klab, arrangement)

    return Kpts, Dvec, Ktik, Klab