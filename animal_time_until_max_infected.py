import functions
beta = 1 #infection rate along edge ij is ijth entry of the matrix U or T times beta (default 1)

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

maxInfectedVecTangoStddev=[]#records the standard deviation of all the quantities for each specific random graph
maxInfectedVecUniformStddev=[]
maxTimeVecTangoStddev=[]
maxTimeVecUniformStddev=[]
maxInfectTimeVecTangoStddev=[]
maxInfectTimeVecUniformStddev=[]
plt.figure(figsize=(14, 10))
colVecT=['b','g','y','k']
colVecU=['r','c','m','0.7']
for i in range(0,len(newGraphList)):
        G=newGraphList[i]
        T=normTangoMatrix(G)
        U=normUniformMatrix(G)
        n=len(list(G.nodes))
        gamma= 1/(10*n)
        print(n)
        secLargestEigTango.append((n-1)*(secLargestEig(normTangoMatrix(G))-(n-2)/(n-1)))
        secLargestEigUniform.append((n-1)*(secLargestEig(normUniformMatrix(G))-(n-2)/(n-1)))
        for e in list(G.edges):
                G.adj[e[0]][e[1]]['weightU']=U[e[0],e[1]]
    
        for e in list(G.edges):
                G.adj[e[0]][e[1]]['weightT']=T[e[0],e[1]]
        
        maxInfectedVecTango=[]
        maxTimeVecTango=[]
        maxInfectTimeVecTango=[]
        maxInfectedVecUniform=[]
        maxTimeVecUniform=[]
        maxInfectTimeVecUniform=[]   
        for j in range(0,iterations):
            
            t, S, I, R = EoN.fast_SIR(G, beta, gamma,rho=None ,transmission_weight='weightT')
            
            while max(I)/n < 0.1:
                t, S, I, R = EoN.fast_SIR(G, beta, gamma,rho=None,transmission_weight='weightT')
            maxInfectedVecTango.append(max(I)/n)
            maxTimeVecTango.append(max(t)/n)
            maxInfectTimeVecTango.append(t[np.argmax(I)]/max(t))
            t, S, I, R = EoN.fast_SIR(G, beta, gamma, rho=None,transmission_weight='weightU')
            
            while max(I)/n < 0.1:
                t, S, I, R = EoN.fast_SIR(G, beta, gamma, rho=None,transmission_weight='weightU')
            maxInfectedVecUniform.append(max(I)/n)
            maxTimeVecUniform.append(max(t)/n)
            maxInfectTimeVecUniform.append(t[np.argmax(I)]/max(t))
            
        vecT=[i for k in range(0,100)]  
        vecU=[i+0.5 for k in range(0,100)]  
        weightsT = [2*k for k in Counter(maxInfectTimeVecTango).values() for m in range(k)]
        weightsU = [2*k for k in Counter(maxInfectTimeVecUniform).values() for m in range(k)]
        if i ==2:
          plt.scatter(vecT,maxInfectTimeVecTango,c='r',label='Degree-based',s=weightsT)
          plt.scatter(vecU,maxInfectTimeVecUniform,c='b',label='Uniform',s=weightsU)
        else:
          plt.scatter(vecT,maxInfectTimeVecTango,c='r',s=weightsT)
          plt.scatter(vecU,maxInfectTimeVecUniform,c='b',s=weightsU)
        
x = np.array([0.25,1.25,2.25])

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

plt.rc('font', **font)

plt.xticks(x, newGraphNames)
plt.xlabel('')
plt.ylabel('Proportion of Time Until Maximum Infected')
plt.title("Proportion of time until maximum infected for wildlife contact networks")
plt.legend()
plt.savefig('TimeAntHyenaTortoise.pdf')            
plt.show()             
            
