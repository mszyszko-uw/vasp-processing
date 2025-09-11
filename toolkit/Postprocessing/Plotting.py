import numpy as np
import matplotlib.pyplot as plt
from InformationCollector import InformationCollector
from Artists import MatrixSculptor, ProjectionPainter, SegmentArranger, KlineEngineer, BandWriter

def BandToGnuplot(Pmat, Emat, Kpts, Dvec, 
                  projection, bands, description, arrangement, 
                  folder=''
                  ):
    """_summary_

    Args:
        Pmat (_type_): _description_
        Emat (_type_): _description_
        Kpts (_type_): _description_
        Dvec (_type_): _description_
        projection (_type_): _description_
        bands (_type_): _description_
        description (_type_): _description_
        arrangement (_type_): _description_
        folder (str, optional): _description_. Defaults to ''.
    """    
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


def DOSToGnuplot(Pdos, Edos, 
                  projection, description, 
                  folder=''
                  ):
    """_summary_

    Args:
        Pmat (_type_): _description_
        Emat (_type_): _description_
        Kpts (_type_): _description_
        Dvec (_type_): _description_
        projection (_type_): _description_
        bands (_type_): _description_
        description (_type_): _description_
        arrangement (_type_): _description_
        folder (str, optional): _description_. Defaults to ''.
    """    
    orbit_names=['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'x2-y2']
    chosen_orbits = ''
    orbitals = projection['orbital'] if projection['orbital'] != ... else np.arange(Pdos.shape[-1])
    directions = projection['direction']
    ions = projection['ion']

    for i in orbitals:
        chosen_orbits += ' '+orbit_names[i]
    with open(folder+description.replace(' ', '_'), 'w') as file:
        file.write(description+'\n')
        header = ' Energy     |'
        # file.write(' Dir 1 Ion 1 px py pz                 : Ion 2 px py pz\n')
        mat = Edos[:,None]
        for in_dir, direction in enumerate(directions):
            stub = f' Dir {in_dir+1} '
            for in_ion, ion in enumerate(ions):
                stub += f': Ion {in_ion+1}' + chosen_orbits
                stub += ' '*(len(orbitals)*13*(in_ion+1) - len(stub)-1)
                mat = np.concatenate([mat, Pdos[:,direction,ion,orbitals]], axis = 1)
            if header: header += stub +'|'
        if header: file.write(header+'\n')
        header = ''
        np.savetxt(file, mat, '%12.8f')
        file.write('\n')


def BandStructurePlot(ic: InformationCollector,
                      min_diff=.1, 
                      ax=None, 
                      E0=None, 
                      gnuplot=False, 
                      folder='', 
                      save=False,
                      bandnums=False, 
                      color = None, 
                      mult=1, 
                      alpha=.1,
                      **kwargs
                      ):
    """_summary_

    Args:
        ic (InformationCollector): _description_
        min_diff (float, optional): _description_. Defaults to .1.
        ax (_type_, optional): _description_. Defaults to None.
        E0 (_type_, optional): _description_. Defaults to None.
        gnuplot (bool, optional): _description_. Defaults to False.
        folder (str, optional): _description_. Defaults to ''.
        save (bool, optional): _description_. Defaults to False.
        bandnums (bool, optional): _description_. Defaults to False.
        color (_type_, optional): _description_. Defaults to None.
        mult (int, optional): _description_. Defaults to 1.
        alpha (float, optional): _description_. Defaults to .1.

    Returns:
        _type_: _description_
    """    
    Pmat = SegmentArranger(ic.nSeg, ic.Pmat, ic.arrangement)
    Emat = SegmentArranger(ic.nSeg, ic.Emat, ic.arrangement)
    kline_args = (ic.Lvec, ic.Scal, ic.Kpts, ic.Klen, ic.Klab, ic.nSeg)
    Kpts, Dvec, Ktik, Klab = KlineEngineer(*kline_args, ic.arrangement)

    Evalmax  = ic.Evalmax
    Egap, ions, num_ions = ic.Egap, ic.ions, ic.num_ions
    description, arrangement, bands = ic.description, ic.arrangement, BandWriter(ic.bands)

    if not folder: folder = ic.folder

    if E0 == None: Emat -= Evalmax
    elif E0 == 'fermi': Emat -= ic.Efer
    else: Emat -= E0

    projection = ProjectionPainter(ions, num_ions, description)

    try: Emat = np.sum(Emat[..., projection['direction']], axis=-1)
    except IndexError: Emat = Emat[...,0]

    if gnuplot: BandToGnuplot(Pmat, Emat, Kpts, Dvec, projection, bands, description, arrangement,folder)

    Pmat = MatrixSculptor(Pmat, projection)
    if ax==None: ax = plt.gca()

    ### Band labelling
    if bandnums:
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

    if 'clean' not in description:
        ax.plot(Dvec, Emat, color="black", linewidth=0.5, linestyle=':')
        ax.scatter(np.stack([Dvec]*Emat.shape[1], axis=1), Emat, np.abs(Pmat)*100*mult, 
                   alpha=alpha, label=description, c=color)
    else: ax.plot(Dvec, Emat, color=color, linewidth=1, label=description)

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


def DensityOfStatesPlot(ic: InformationCollector,
                        ax=None, E0=None, gnuplot=False,
                        **kwargs
                        ):
    """_summary_

    Args:
        ic (InformationCollector): _description_
        ax (_type_, optional): _description_. Defaults to None.
        E0 (_type_, optional): _description_. Defaults to None.
    """    
    Vdos, Pdos, Edos, Efer = ic.Vdos, ic.Pdos, ic.Edos, ic.Efer
    ions, num_ions, description = ic.ions, ic.num_ions, ic.description
    Emin, Emax = ic.Emin, ic.Emax

    if E0 == None: 
        Edos -= Efer
        Emin -= Efer
        Emax -= Efer
    else: 
        Edos -= E0
        Emin -= E0
        Emax -= E0

    projection = ProjectionPainter(ions, num_ions, description)

    if gnuplot: DOSToGnuplot(Pdos, Edos,
                             projection, description, ic.folder)

    Pdos = MatrixSculptor(Pdos, projection)

    Vdos = np.sum(Vdos[..., projection['direction']], axis=-1)

    if ax==None: ax = plt.gca()
    if description in ['total', 'clean']:
        ax.plot(Vdos, Edos, label='total')
        ax.set_xlim([np.percentile(Vdos, .1),np.percentile(Vdos, 99.9)])
    else:
        ax.plot(Pdos, Edos, label=description)
        ax.set_xlim([np.percentile(Pdos, .1),np.percentile(Pdos, 99.9)])

    
    ax.set_ylim([Emin, Emax])