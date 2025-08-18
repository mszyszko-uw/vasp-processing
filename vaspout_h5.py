import h5py
import numpy as np
import matplotlib.pyplot as plt
from matrix_postprocess import BandStructurePlot, DensityOfStatesPlot
from matrix_artists import MatrixSculptor, SegmentArranger, KlineEngineer, ProjectionPainter

class vaspout_h5:
    def __init__(self, vaspout: str):
        self.vaspout = h5py.File(vaspout, 'r')
        self.folder = vaspout[:vaspout.rindex('/')+1]

    def read_BS(self, arrangement:str = ''):
        vaspout = self.vaspout
        Pmat = vaspout['results']['projectors']['par'][:]
        Emat = vaspout['results']['electron_eigenvalues']['eigenvalues'][:]
        Occp = vaspout['results']['electron_eigenvalues']['fermiweights'][:]
        nSeg = vaspout['input']['kpoints']['number_line_segments'][()]/2
        Lvec = vaspout['input']['poscar']['lattice_vectors'][:]
        Scal = vaspout['input']['poscar']['scale'][()]
        Kpts = vaspout['results']['electron_eigenvalues']['kpoint_coords'][:]
        Klen = vaspout['input']['kpoints']['number_kpoints'][()]
        Klab = vaspout['input']['kpoints']['labels_kpoints'].asstr()[:]

        Pmat = np.moveaxis(Pmat, [-2,-1], [0,1])
        Pmat = SegmentArranger(nSeg, Pmat, arrangement)

        Evalmax = np.max(Emat*Occp)
        Econmin = np.min(Emat[Occp==0])
        Egap = Econmin-Evalmax

        Emat = np.moveaxis(Emat, [-2,-1], [0,1])
        Emat = SegmentArranger(nSeg, Emat, arrangement)

        Kpts, Dvec, Ktik, Klab = KlineEngineer(Lvec, Scal, Kpts, Klen, Klab, nSeg, arrangement)

        return Pmat, Emat, Kpts, Dvec, Ktik, Klab, Evalmax, Egap


    def read_DOS(self):
        vaspout = self.vaspout

        Vdos = vaspout['results']['electron_dos']['dos'][:]
        Vdos = np.moveaxis(Vdos, [-1], [0])

        Pdos = vaspout['results']['electron_dos']['dospar'][:]
        Pdos = np.moveaxis(Pdos, [-1], [0])

        Edos = vaspout['results']['electron_dos']['energies'][:]

        Efer = vaspout['results']['electron_dos']['efermi'][()]
        return Vdos, Pdos, Edos, Efer        
    

    def plot_BS(self,  
                description: str = '', 
                arrangement: str = '', 
                bands: list = ..., 
                min_diff=.1, 
                ax=None, 
                E0=None, 
                gnuplot=False, 
                folder=''
                ):
        if not folder: folder = self.folder
        ions = self.vaspout['input']['poscar']['ion_types'].asstr()[:]
        num_ions = self.vaspout['input']['poscar']['number_ion_types'][:]
        Pmat, Emat, Kpts, Dvec, Ktik, Klab, Evalmax, Egap = self.read_BS(arrangement)
        return BandStructurePlot(Pmat, Emat, Kpts, Dvec, Ktik, Klab, Evalmax, Egap, ions, num_ions, 
                                 description, arrangement, bands, min_diff, ax, E0, gnuplot, folder)
    

    def plot_DOS(self, 
                 description='', 
                 ax=None, 
                 E0=None
                 ):
        ions = self.vaspout['input']['poscar']['ion_types'].asstr()[:]
        num_ions = self.vaspout['input']['poscar']['number_ion_types'][:]
        Vdos, Pdos, Edos, Efer = self.read_DOS()
        return DensityOfStatesPlot(Vdos, Pdos, Edos, Efer, ions, num_ions, 
                                   description, ax=ax, E0=E0)
    

def plot_BS_DOS(BS:vaspout_h5, 
                DOS:vaspout_h5, 
                description='', 
                arrangement='', 
                bands: list = ...
                ):
    fig, axes = plt.subplots(ncols=2, sharey=True, gridspec_kw={'width_ratios': [3, 1]})
    Ediff, Efer = DOS.plot_DOS(description, axes[1])
    Evalmax, Egap = BS.plot_BS(description, arrangement, bands, Ediff*1e-2, axes[0])
    axes[1].clear()
    DOS.plot_DOS(description, ax=axes[1], E0=Evalmax)
    plt.tight_layout()