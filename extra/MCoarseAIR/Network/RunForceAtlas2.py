import networkx as nx
#import matplotlib.pyplot as plt
import numpy as np
import csv
from fa2 import ForceAtlas2
import time
#from infomap import infomap

System       = 'O3'
Molecule     = 'O2'
RatesType    = 'Inel'
PathToResults= '/home/venturi/WORKSPACE/CG-QCT/RESULTS/' + System + '/Postprocessing/' + Molecule + '_Network/'

t0 = time.time()
G = nx.DiGraph(Rates=RatesType)

#G      = nx.read_pajek("/Users/sventuri/Dropbox/TempRes/InelRates.net")
#G      = nx.read_pajek("/home/venturi/WORKSPACE/CG-QCT/RESULTS/CO2/4Networks/InelExchange/InelRates1.net")

Id        = []
vqn       = []
Longitude = []
Latitude  = []
rIn       = []  
EeVVib    = []
EeVRot    = []
#with open('/Users/sventuri/Dropbox/TempRes/InelLevels.csv') as csvDataFile:
FileName = PathToResults + RatesType + 'Levels.csv'
with open(FileName) as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    i = -1
    for row in csvReader:
        if i>-1:
            Id        = np.append(Id,int(row[0]))
            vqn       = np.append(vqn,int(row[1]))
            Longitude = np.append(Longitude,float(row[2]))
            Latitude  = np.append(Latitude,float(row[3]))
            rIn       = np.append(rIn,float(row[4]))
            EeVVib    = np.append(EeVVib,float(row[5]))
            EeVRot    = np.append(EeVRot,float(row[6]))
            G.add_node(Id[i], v=int(vqn[i]), Longitude=float(Longitude[i]), Latitude=float(Latitude[i]), rIn=float(rIn[i]), EeVVib=float(EeVVib[i]), EeVRot=float(EeVRot[i]))
        i=i+1


Source = []
End    = []
Weight = []
#with open('/Users/sventuri/Dropbox/TempRes/Edges.net') as csvDataFile:
FileName = PathToResults + RatesType + 'Edges.net'
with open(FileName) as csvDataFile:
    csvReader = csv.reader(csvDataFile, delimiter=' ')
    #i=0
    for row in csvReader:
        Source = int(row[0])
        End    = int(row[1])
        Weight = float(row[2])
        G.add_edge(Source, End, weight=float(Weight))
        #Source     = np.append(Source,int(row[0]))
        #End        = np.append(End,int(row[1]))
        #Weight     = np.append(Weight,float(row[2]))
        #G.add_edge(Source[i], End[i], weight=Weight[i])
        #if Source%500 == 0:
        #    print(i)
        #i=i+1
LoadingTime = time.time() - t0
print('Loading done. Elapsed time: %s' % LoadingTime)



t0 = time.time()
forceatlas2 = ForceAtlas2(
                          # Behavior alternatives
                          outboundAttractionDistribution=False,  # Dissuade hubs
                          linLogMode=False,  # NOT IMPLEMENTED
                          adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
                          edgeWeightInfluence=1.0,

                          # Performance
                          jitterTolerance=1.e-7,  # Tolerance
                          barnesHutOptimize=False,
                          barnesHutTheta=1.2,
                          multiThreaded=False,  # NOT IMPLEMENTED

                          # Tuning
                          scalingRatio=1.e7,
                          strongGravityMode=False,
                          gravity=1.0,

                          # Log
                          verbose=True)
AtlasLoadingTime = time.time() - t0
print('Atlas is done Loading. Elapsed time: %s' % AtlasLoadingTime)


t0 = time.time()
positions = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=20000)

i = 0
X = []
Y = []
for node,(x,y) in positions.items():
    G.node[node]['x']         = float(x)
    G.node[node]['y']         = float(y)
    X = np.append(X,float(x))
    Y = np.append(Y,float(y))
    i    = i+1

#nx.write_graphml(G, "G_2000.graphml")
FileName = PathToResults + RatesType + '_G20000.csv'
with open(FileName, mode='w') as csvFile:
    csvWriter = csv.writer(csvFile, delimiter=',')
    csvWriter.writerow(['id','v','Longitude','Latitude','rIn','EeVVib','EeVRot','x','y'])
    for i in range(0, len(vqn)):
        csvWriter.writerow([str(i+1), str(vqn[i]), str(Longitude[i]), str(Latitude[i]), str(rIn[i]), str(EeVVib[i]), str(EeVRot[i]), str(X[i]), str(Y[i])])
Atlas3Time = time.time() - t0
print('Atlas 20000 is done. Elapsed time: %s' % Atlas3Time)



#nx.draw_networkx(G, positions, cmap=plt.get_cmap('jet'), node_size=50, with_labels=False)
#plt.show()

#G.get_edge_data('Name 1','Name 2')

