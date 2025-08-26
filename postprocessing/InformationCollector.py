class InformationCollector:
    """A class to collect all relevant DFT information
    """    
    def __init__(self):
        ## General
        self.description = ''
        self.folder = ''

        ## Band Structure
        self.arrangement = ''
        self.bands = ''
        self.Pmat = [[[[[]]]]]  # shape ( nkpts, nbands, ndirs, nions, norbits )
        self.Emat = [[[]]]      # shape ( nkpts, nbands, 2 if SpinPolarized else 1  )
        self.Occp = [[]]        # shape ( nkpts, nbands )
        self.Lvec = []          # shape ( nkpts )
        self.Kpts = []          # shape ( nkpts )
        self.Dvec = []          # shape ( nkpts )
        self.Klab = []
        self.Ktik = []
        self.nSeg = 1
        self.Scal = 1
        self.Klen = 1
        self.Evalmax = 0
        self.Econmin = 0
        self.Egap = 0
        
        # Density of States
        self.Vdos = []          # shape ( nedos )
        self.Pdos = [[[[]]]]    # shape ( nedos, ndirs, nions, norbits )
        self.Edos = []          # shape ( nedos )
        self.ions = []          # list
        self.num_ions = []      # list
        self.Efer = 0
        self.Emin = 0
        self.Emax = 0
        self.Ediff = 0