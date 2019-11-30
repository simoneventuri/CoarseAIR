import numpy
import pandas

PathToInputFldr       = '/Users/sventuri/WORKSPACE/CoarseAIR/coarseair/dtb/N3/PESs/BNN/5000Points/CalibratedParams/'
PathToOutputFldrPlus  = '/Users/sventuri/WORKSPACE/CoarseAIR/coarseair/dtb/N3/PESs/NN/5000PointsPlus1Sigma/CalibratedParams/'
PathToOutputFldrMinus = '/Users/sventuri/WORKSPACE/CoarseAIR/coarseair/dtb/N3/PESs/NN/5000PointsMinus1Sigma/CalibratedParams/'

PathToFinalW = PathToInputFldr + '/W1_mu.csv'
WData = pandas.read_csv(PathToFinalW, header=None)
WData = WData.apply(pandas.to_numeric, errors='coerce')
#yData  = numpy.transpose(yyData.values)
W1_MU  = WData.values

PathToFinalb = PathToInputFldr + '/b1_mu.csv'
bData = pandas.read_csv(PathToFinalb, header=None)
bData = bData.apply(pandas.to_numeric, errors='coerce')
#yData  = numpy.transpose(yyData.values)
b1_MU  = bData.values[:,0]


PathToFinalW = PathToInputFldr + '/W2_mu.csv'
WData = pandas.read_csv(PathToFinalW, header=None)
WData = WData.apply(pandas.to_numeric, errors='coerce')
#yData  = numpy.transpose(yyData.values)
W2_MU  = WData.values

PathToFinalb = PathToInputFldr + '/b2_mu.csv'
bData = pandas.read_csv(PathToFinalb, header=None)
bData = bData.apply(pandas.to_numeric, errors='coerce')
#yData  = numpy.transpose(yyData.values)
b2_MU  = bData.values[:,0]


PathToFinalW = PathToInputFldr + '/W3_mu.csv'
WData = pandas.read_csv(PathToFinalW, header=None)
WData = WData.apply(pandas.to_numeric, errors='coerce')
#yData  = numpy.transpose(yyData.values)
W3_MU  = WData.values

PathToFinalb = PathToInputFldr + '/b3_mu.csv'
bData = pandas.read_csv(PathToFinalb, header=None)
bData = bData.apply(pandas.to_numeric, errors='coerce')
#yData  = numpy.transpose(yyData.values)
b3_MU  = bData.values[:,0]


PathToFinalb = PathToInputFldr + '/Sigma_mu.csv'
bData = pandas.read_csv(PathToFinalb, header=None)
bData = bData.apply(pandas.to_numeric, errors='coerce')
#yData  = numpy.transpose(yyData.values)
Sigma_MU  = bData.values[:,0]





PathToFinalW = PathToInputFldr + '/W1_sd.csv'
WData = pandas.read_csv(PathToFinalW, header=None)
WData = WData.apply(pandas.to_numeric, errors='coerce')
#yData  = numpy.transpose(yyData.values)
W1_SD  = WData.values

PathToFinalb = PathToInputFldr + '/b1_sd.csv'
bData = pandas.read_csv(PathToFinalb, header=None)
bData = bData.apply(pandas.to_numeric, errors='coerce')
#yData  = numpy.transpose(yyData.values)
b1_SD  = bData.values[:,0]


PathToFinalW = PathToInputFldr + '/W2_sd.csv'
WData = pandas.read_csv(PathToFinalW, header=None)
WData = WData.apply(pandas.to_numeric, errors='coerce')
#yData  = numpy.transpose(yyData.values)
W2_SD  = WData.values

PathToFinalb = PathToInputFldr + '/b2_sd.csv'
bData = pandas.read_csv(PathToFinalb, header=None)
bData = bData.apply(pandas.to_numeric, errors='coerce')
#yData  = numpy.transpose(yyData.values)
b2_SD  = bData.values[:,0]


PathToFinalW = PathToInputFldr + '/W3_sd.csv'
WData = pandas.read_csv(PathToFinalW, header=None)
WData = WData.apply(pandas.to_numeric, errors='coerce')
#yData  = numpy.transpose(yyData.values)
W3_SD  = WData.values

PathToFinalb = PathToInputFldr + '/b3_sd.csv'
bData = pandas.read_csv(PathToFinalb, header=None)
bData = bData.apply(pandas.to_numeric, errors='coerce')
#yData  = numpy.transpose(yyData.values)
b3_SD  = bData.values[:,0]


PathToFinalb = PathToInputFldr + '/Sigma_sd.csv'
bData = pandas.read_csv(PathToFinalb, header=None)
bData = bData.apply(pandas.to_numeric, errors='coerce')
#yData  = numpy.transpose(yyData.values)
Sigma_SD  = bData.values[:,0]



W1 = W1_MU + 1.0*W1_SD
W2 = W2_MU + 1.0*W2_SD
W3 = W3_MU + 1.0*W3_SD

b1 = b1_MU + 1.0*b1_SD
b2 = b2_MU + 1.0*b2_SD
b3 = b3_MU + 1.0*b3_SD

Sigma = Sigma_MU + 1.0*Sigma_SD


PathToFinalW = PathToOutputFldrPlus + '/W1.csv'
numpy.savetxt(PathToFinalW, W1, delimiter=",")

PathToFinalb = PathToOutputFldrPlus + '/b1.csv'
numpy.savetxt(PathToFinalb, b1, delimiter=",")

PathToFinalW = PathToOutputFldrPlus + '/W2.csv'
numpy.savetxt(PathToFinalW, W2, delimiter=",")

PathToFinalb = PathToOutputFldrPlus + '/b2.csv'
numpy.savetxt(PathToFinalb, b2, delimiter=",")

PathToFinalW = PathToOutputFldrPlus + '/W3.csv'
numpy.savetxt(PathToFinalW, W3, delimiter=",")

PathToFinalb = PathToOutputFldrPlus + '/b3.csv'
numpy.savetxt(PathToFinalb, b3, delimiter=",")

PathToFinalSigma = PathToOutputFldrPlus + '/Sigma.csv'
numpy.savetxt(PathToFinalSigma, Sigma, delimiter=",")




W1 = W1_MU - 1.0*W1_SD
W2 = W2_MU - 1.0*W2_SD
W3 = W3_MU - 1.0*W3_SD

b1 = b1_MU - 1.0*b1_SD
b2 = b2_MU - 1.0*b2_SD
b3 = b3_MU - 1.0*b3_SD

Sigma = Sigma_MU - 1.0*Sigma_SD


PathToFinalW = PathToOutputFldrMinus + '/W1.csv'
numpy.savetxt(PathToFinalW, W1, delimiter=",")

PathToFinalb = PathToOutputFldrMinus + '/b1.csv'
numpy.savetxt(PathToFinalb, b1, delimiter=",")

PathToFinalW = PathToOutputFldrMinus + '/W2.csv'
numpy.savetxt(PathToFinalW, W2, delimiter=",")

PathToFinalb = PathToOutputFldrMinus + '/b2.csv'
numpy.savetxt(PathToFinalb, b2, delimiter=",")

PathToFinalW = PathToOutputFldrMinus + '/W3.csv'
numpy.savetxt(PathToFinalW, W3, delimiter=",")

PathToFinalb = PathToOutputFldrMinus + '/b3.csv'
numpy.savetxt(PathToFinalb, b3, delimiter=",")

PathToFinalSigma = PathToOutputFldrMinus + '/Sigma.csv'
numpy.savetxt(PathToFinalSigma, Sigma, delimiter=",")