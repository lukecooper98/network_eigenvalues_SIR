import functions
path = 'Animal_Networks/*.edges' #change this directory to wherever you're storing the networks

files = glob.glob(path)
graphList = []
graphNames = []
for name in files:
    try:
        with open(name) as f:
            firstLine=f.readline()
            if firstLine.count(" ") > 0 and firstLine.count(" ") <= 2 and (firstLine[0] =='1' or firstLine[0] =='0'):
              G=max(connected_component_subgraphs(nx.read_weighted_edgelist(f,nodetype=int)),key=len)
              
              graphList.append(nx.convert_node_labels_to_integers(G)) #all networks will be stored in graphList after this code is run
              graphNames.append(name) 
              print(name)
    except IOError as exc:
        if exc.errno != errno.EISDIR:
            raise

tangoList=[]
unifList=[]
for i in range(0,len(graphList)):
  
  n=len(list(graphList[i].nodes))
  tangoList.append((n-1)*(secLargestEig(normTangoMatrix(graphList[i]))-(n-2)/(n-1)))
  unifList.append((n-1)*(secLargestEig(normUniformMatrix(graphList[i]))-(n-2)/(n-1)))

newGraphList=[graphList[17],graphList[np.argmax(unifList)],graphList[np.argmin(tangoList)]]
newGraphNames=['Ant','Tortoise','Hyena']
