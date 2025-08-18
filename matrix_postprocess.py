import numpy as np
import matplotlib.pyplot as plt
from matrix_artists import MatrixSculptor, SegmentArranger, KlineEngineer, ProjectionPainter

def BandToGnuplot(Pmat, Emat, Kpts, Dvec, 
                  projection, bands, description, arrangement, 
                  folder=''
                  ):
    orbit_names=['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'x2-y2']
    chosen_orbits = ''
    orbitals = projection['orbital'] if projection['orbital'] != ... else np.arange(Pmat.shape[-1])
    directions = projection['direction']
    ions = projection['ion']
    bands = bands if bands != ... else np.arange(Pmat.shape[1])

    for i in orbitals:
        chosen_orbits += ' '+orbit_names[i]
    Dvec = Dvec.copy()
    Dvec -= Dvec[0]
    with open(folder+description.replace(' ', '_'), 'w') as file:
        file.write(description+' '+arrangement+' '+(str(bands) if len(bands)<Pmat.shape[1] else '')+'\n')
        header = 'k point coordinates x y z             | Distance   | Energy     |'
        # file.write(' Dir 1 Ion 1 px py pz                 : Ion 2 px py pz\n')
        for band in bands:
            mat = np.concatenate([Kpts, Dvec[:,None], Emat[:,band][:,None]], axis = 1)
            for in_dir, direction in enumerate(directions):
                stub = f' Dir {in_dir+1} '
                for in_ion, ion in enumerate(ions):
                    stub += f': Ion {in_ion+1}' + chosen_orbits
                    stub += ' '*(len(orbitals)*13*(in_ion+1) - len(stub)-1)
                    mat = np.concatenate([mat, Pmat[:,band,direction,ion,orbitals]], axis = 1)
                if header: header += stub +'|'
            if header: file.write(header+'\n')
            header = ''
            np.savetxt(file, mat, '%12.8f')
            file.write('\n')


def BandStructurePlot(Pmat, Emat, Kpts, Dvec, 
                      Ktik, Klab, 
                      Evalmax, Egap, 
                      ions, num_ions,
                      description: str = '', arrangement: str = '', bands: list = ..., 
                      min_diff=.1, ax=None, E0=None, gnuplot=False, folder='', save=False
                      ):
    # Pmat, Emat, Kpts, Dvec, Ktik, Klab, Evalmax, Egap = read_BS(vaspout, arrangement)
    if E0 == None: Emat -= Evalmax
    else: Emat -= E0

    projection = ProjectionPainter(ions, num_ions, description)

    try: Emat = np.sum(Emat[..., projection['direction']], axis=-1)
    except IndexError: Emat = Emat[...,0]

    if gnuplot: BandToGnuplot(Pmat, Emat, Kpts, Dvec, projection, bands, description, arrangement,folder)

    Pmat = MatrixSculptor(Pmat, projection)
    if ax==None: ax = plt.gca()

    ### Band labelling
    diff = Emat[-1, 0]
    band_prev = 0
    for band in range(Pmat.shape[1]):
        diff = abs(Emat[-1, band] - diff)
        if diff > min_diff:
            if band-1 != band_prev:
                ax.annotate(str(band_prev+1)+'-'+str(band-1+1), xy = (Dvec[-1], Emat[-1, band-1]))
            else: ax.annotate(str(band_prev+1), xy = (Dvec[-1], Emat[-1, band-1]))
            band_prev = band
        diff = Emat[-1, band]
    band = Pmat.shape[1]
    if band-1 != band_prev:
        ax.annotate(str(band_prev+1)+'-'+str(band-1+1), xy = (Dvec[-1], Emat[-1, band-1]))
    else: ax.annotate(str(band_prev+1), xy = (Dvec[-1], Emat[-1, band-1]))

    Emat = Emat[:,bands]
    Pmat = Pmat[:,bands]

    ax.plot(Dvec, Emat, color="black", linewidth=0.5, linestyle=':')
    no_proj_kwds = ['clean', 'clear', 'pure']
    if description not in no_proj_kwds:
        # print(np.stack([Dvec]*Emat.shape[1], axis=1).shape, Emat.shape, Pmat.shape)
        ax.scatter(np.stack([Dvec]*Emat.shape[1], axis=1), Emat, np.maximum(0,Pmat)*100, 
                   alpha=0.10, label=description)
    ax.vlines(Ktik[1:-1], [Emat.min()]*(len(Ktik)-2), [Emat.max()]*(len(Ktik)-2),
            color="black", linestyle=":")
    ax.axhline(0, color="black", linestyle=":")
    if E0==None:
        ax.axhline(Egap, color="black", linestyle=":")
        ax.annotate(f'{Egap:.3f}', xy = (Dvec[0], Egap))
    ax.set_xticks(Ktik, Klab)
    ax.set_xlim(Dvec[0], Dvec[-1])
    if save: plt.savefig(folder+description.replace(' ', '_'))
    return Evalmax, Egap


def DensityOfStatesPlot(Vdos, Pdos, Edos, Efer, 
                        ions, num_ions,
                        description, 
                        ax=None, E0=None
                        ):
    # Vdos, Pdos, Edos, Efer = read_DOS(DOS_vaspout)
    if E0 == None: Edos -= Efer
    else: Edos -= E0

    projection = ProjectionPainter(ions, num_ions, description)
    Pdos = MatrixSculptor(Pdos, projection)

    Vdos = np.sum(Vdos[..., projection['direction']], axis=-1)

    if ax==None: ax = plt.gca()
    if description == 'total':
        ax.plot(Vdos, Edos, label=description)
        ax.set_xlim([np.percentile(Vdos, .1),np.percentile(Vdos, 99.9)])
    else:
        ax.plot(Pdos, Edos, label=description)
        ax.set_xlim([np.percentile(Pdos, .1),np.percentile(Pdos, 99.9)])
    
    Emin = np.min(Edos)
    Emax = np.max(Edos)

    
    ax.set_ylim([Emin, Emax])

    return Emax - Emin, Efer