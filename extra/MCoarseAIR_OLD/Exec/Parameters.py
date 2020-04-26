class PhysParams(object):

    def __init__():
            
        self.Plnck    = 6.62607004e-34
        self.UKb      = 1.380658e-23
        self.Ue       = 1.602191e-19
        self.KeV      = 8.617330e-05
        self.AvN      = 6.0221409e+23
        self.AMUToKg  = 1.e0 / AvN * 1.e-3
        self.EhToeV   = 27.2114e0
        self.DSWtoKg  = 1.e-3 / 1.8208e+03
        self.ATMToPa  = 101325.e0



class Paths(object):

    def __init__():

        self.System   = Input.PathToOutput + '/' + Input.System                                             # For qnsEnBin.dat (Levels Info)
        self.Rates    = Input.PathToOutput + '/' + Input.System   + '/' + Input.MoleculesName + '/Rates/'   # For Levels / Bins Rates
        self.Database = Input.PathToOutput + '/' + Input.RunKONIG + '/database/'                            # For thermo (Bins Info)
        self.Output   = Input.PathToOutput + '/' + Input.RunKONIG + '/output/'                              # For box.dat, Tint.dat, 

        Paths.RunKONIG
        Paths.Figs


class PlotParams(object):

    def __init__():

        self.linS          = {'--','-.',':','-'};
        self.linST         = {'-','-.',':','--'};

        self.AxisFontSz    = 36;
        self.AxisFontNm    = 'Arial';
        self.AxisLabelSz   = 44;
        self.AxisLabelNm   = 'Arial';
        self.LegendFontSz  = 40;
        self.LegendFontNm  = 'Arial';

        self.XLimPlot
        self.YLimPlot 

        self.RCVec = [255  50  20] ./ 255;
        self.BCVec = [  0  70 200] ./ 255;
        self.GCVec = [  0 140  50] ./ 255;
        self.KCVec = [  0   0   0] ./ 255;
        self.OCVec = [255 105  45] ./ 255;
        self.PCVec = [155  45 175] ./ 255;
        self.WCVec = [  1   1   1] ./ 255;
        self.JCVec = [100 100 100] ./ 255;
        self.YCVec = [255 255   0] ./ 255;
        self.CCVec = [205 205 205] ./ 255;
        self.MCVec = [100  25  15] ./ 255;


        