
n=100 #fixed number of vertices for all randomly generated graphs

beta = 1 #infection rate along edge ij is ijth entry of the matrix U or T times beta (default 1)
gamma= 0.001 #constant recovery rate in SIR simulations
iterations=100 #number of simulations run for each randomly generated network
maxInfectedVecTango=[] #records the epidemic size for degree-based and uniform weights for each simulation
maxInfectedVecUniform=[]
secLargestEigTango=[] #records second largest eigenvalue of degree-based and uniform for each randomly generated graph
secLargestEigUniform=[]
maxInfectedVecTangoAverage=[] #records the average epidemic size over all simulations on a specific random graph
maxInfectedVecUniformAverage=[]
maxTimeVecTango=[] #records the simulation time for degree-based and uniform weights for each simulation
maxTimeVecUniform=[]
maxTimeVecTangoAverage=[]#records the average simulation time over all simulations on a specific random graph
maxTimeVecUniformAverage=[]
maxInfectTimeVecTango=[]#records the epidemic time for degree-based and uniform weights for each simulation
maxInfectTimeVecUniform=[]
maxInfectTimeVecUniformAverage=[]#records the average epidemic time over all simulations on a specific random graph
maxInfectTimeVecTangoAverage=[]
mVec=[m for m in range(1,n)] #range of m values to consider (BA parameter) (here the full range is used)
for m in mVec:
    for i in range(0,10):#generate 10 random graphs for each value of p
        G=nx.barabasi_albert_graph(n, m)
        while not nx.is_connected(G):
            G=nx.barabasi_albert_graph(n, m)#keep generating graphs until a connected graph is created
    
        print(m)#print each p to record code progress
        T=normTangoMatrix(G) #define degree-based and uniform matrix for each graph
        U=normUniformMatrix(G)
        secLargestEigTango.append(secLargestEig(T)) #get second largest eigs of degree-based and uniform matrices
        secLargestEigUniform.append(secLargestEig(U))
        for e in list(G.edges):
                G.adj[e[0]][e[1]]['weightU']=U[e[0],e[1]] #set the edge weights of the graph to the entries of uniform and degree-based matrices
    
        for e in list(G.edges):
                G.adj[e[0]][e[1]]['weightT']=T[e[0],e[1]]
        
        maxInfectedVecTango=[] #reset all of the following lists for each individual randomly generated graph
        maxTimeVecTango=[]
        maxInfectTimeVecTango=[]
        maxInfectedVecUniform=[]
        maxTimeVecUniform=[]
        maxInfectTimeVecUniform=[]     
        for i in range(0,iterations): #do 100 simulations on each graph
            #t is a list that records running time of simulation, S the number of susceptible at each t, I the number of infected at each t, and R the number of recovered at each t
            #a new time is added to the list t each time a transition between S, I and R occurs
            t, S, I, R = EoN.fast_SIR(G, beta, gamma,rho=None ,transmission_weight='weightT')
            
            maxInfectedVecTango.append(max(I)) #epidemic size is max of vector I
            maxTimeVecTango.append(max(t)) #simulation time is the largest (and last) entry in t
            maxInfectTimeVecTango.append(t[np.argmax(I)]) #epidemic time is the entry of t in the position of the largest entry of I

            t, S, I, R = EoN.fast_SIR(G, beta, gamma, rho=None,transmission_weight='weightU')
            
            maxInfectedVecUniform.append(max(I))
            maxTimeVecUniform.append(max(t))
            maxInfectTimeVecUniform.append(t[np.argmax(I)])
        
        #take averages for each quantity over all simulations on the graph
        maxInfectedVecTangoAverage.append(np.mean(maxInfectedVecTango))
        maxInfectedVecUniformAverage.append(np.mean(maxInfectedVecUniform))
        maxTimeVecTangoAverage.append(np.mean(maxTimeVecTango))
        maxTimeVecUniformAverage.append(np.mean(maxTimeVecUniform))
        maxInfectTimeVecTangoAverage.append(np.mean(maxInfectTimeVecTango))
        maxInfectTimeVecUniformAverage.append(np.mean(maxInfectTimeVecUniform))

#this just fixes the fact that 10 graphs are produced for each m (otherwise mVec and secLargestEigTango have different lengths)
mmVec=[]
for m in mVec:
    for i in range(0,10):
        mmVec.append(m)

#plotting
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

plt.rc('font', **font)
plt.figure(figsize=(14, 10))
plt.plot(mmVec,secLargestEigTango,'r.',label='Degree Dependent')
plt.plot(mmVec,secLargestEigUniform,'b.',label='Uniform')
plt.ylabel('$\lambda_2$')
plt.xlabel('$m_{BA}$')
plt.title("Number of Edge Attachments ($m_{BA}$) against $\lambda_2$")
plt.legend()
plt.savefig('mAgainstLambda2100verticesBA.pdf')

plt.figure(figsize=(14, 10))
plt.plot(secLargestEigTango,maxInfectedVecTangoAverage,'r.',label='Degree Dependent')
plt.plot(secLargestEigUniform,maxInfectedVecUniformAverage,'b.',label='Uniform')
plt.xlabel('$\lambda_2$')
plt.ylabel('Number infected')
plt.title("$\lambda_2$ against Epidemic Size with $\gamma=0.001$")
plt.legend()
plt.savefig('Epidemic_Against_Eig100verticesBAgamma0_001iterations100mincent.pdf')

plt.figure(figsize=(14, 10))
plt.plot(secLargestEigTango,maxInfectTimeVecTangoAverage,'r.',label='Degree Dependent')
plt.plot(secLargestEigUniform,maxInfectTimeVecUniformAverage,'b.',label='Uniform')
plt.xlabel('$\lambda_2$')
plt.ylabel('Time Until Max Infected')
plt.title("$\lambda_2$ against Time Until Max Infected with $\gamma=0.001$")
plt.legend()
plt.savefig('Time_Against_Eig100verticesBAgamma0_001iterations100mincent.pdf')

plt.figure(figsize=(14, 10))
plt.plot(secLargestEigTango,maxTimeVecTangoAverage,'r.',label='Degree Dependent')
plt.plot(secLargestEigUniform,maxTimeVecUniformAverage,'b.',label='Uniform')
plt.xlabel('$\lambda_2$')
plt.ylabel('Convergence Time')
plt.title("$\lambda_2$ against Convergence Time with $\gamma=0.001$")
plt.legend()
plt.savefig('ConvTime_Against_Eig100verticesBAgamma0_001iterations100mincent.pdf')
