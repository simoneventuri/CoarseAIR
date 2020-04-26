import numpy

class Input(object):

    PathToOutput  = 
    System        =
    T0_Vec        =
    NBins         =
    KinMthd       =
    NTtra         =
    NTint         =
    TInt_Vec      =
    NMolecules    =
    StartBin      =
    FinalBin      =
    NBinnedMol    =
    BinnedMolName =

    SystemPath    = PathToOutput + System
    RatesPath     = PathToOutput + System + '/N2/Rates/'         
    RunKONIGFldr  = 'RunKonig_Inel_w_Exch'
    OutputPath    = PathToOutput + RunKONIGFldr + '/output_0595/'

    iFigure       = 2001
    SaveFigs      = 2
    FigsPath      = './CGQCT-' + System + '-Figures/'

    def __init__(self, 
                 PathToOutput,
                 System,
                 T0_Vec,
                 NBins,
                 KinMthd,
                 NTtra,
                 NTint,
                 TInt_Vec,
                 NMolecules,
                 StartBin,
                 FinalBin,
                 NBinnedMol,
                 BinnedMolName,
                 SystemPath,
                 RatesPath,
                 RunKONIGFldr,
                 OutputPath,
                 iFigure,
                 SaveFig,
                 FigsPath):
            
        self.PathToOutput  = PathToOutput
        self.System        = System
        self.T0_Vec        = T0_Vec
        self.NBins         = NBins
        self.KinMthd       = KinMthd
        self.NTtra         = NTtra
        self.NTint         = NTint
        self.TInt_Vec      = TInt_Vec
        self.NMolecules    = NMolecules
        self.StartBin      = StartBin
        self.FinalBin      = FinalBin
        self.NBinnedMol    = NBinnedMol
        self.BinnedMolName = BinnedMolName 

        self.SystemPath    = SystemPath
        self.RatesPath     = RatesPath
        self.RunKONIGFldr  = RunKONIGFldr
        self.OutputPath    = OutputPath

        self.iFigure       = iFigure
        self.SaveFigs      = SaveFigs
        self.FigsPath      = FigsPath