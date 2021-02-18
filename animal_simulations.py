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

maxInfectedVecTangoStddev=[] #records the standard deviation of all the quantities for each specific random graph
maxInfectedVecUniformStddev=[]
maxTimeVecTangoStddev=[]
maxTimeVecUniformStddev=[]
maxInfectTimeVecTangoStddev=[]
maxInfectTimeVecUniformStddev=[]

for G in graphList:
    
        T=normTangoMatrix(G) #create the degree-based and uniform matrices for each graph G
        U=normUniformMatrix(G)
        n=len(list(G.nodes)) #get the number of vertices of G
        gamma= 1/(10*n) #define the recovery rate gamma according to number of vertices
        print(n) #print number of vertices to track code progress
        secLargestEigTango.append((n-1)*(secLargestEig(normTangoMatrix(G))-(n-2)/(n-1))) #record the normalised eigenvalues of degree-based and uniform
        secLargestEigUniform.append((n-1)*(secLargestEig(normUniformMatrix(G))-(n-2)/(n-1)))
        for e in list(G.edges): #apply the matrix entries to the edges of the graph
                G.adj[e[0]][e[1]]['weightU']=U[e[0],e[1]]
    
        for e in list(G.edges):
                G.adj[e[0]][e[1]]['weightT']=T[e[0],e[1]]
        
        #reset these vectors for each graph
        maxInfectedVecTango=[]
        maxTimeVecTango=[]
        maxInfectTimeVecTango=[]
        maxInfectedVecUniform=[]
        maxTimeVecUniform=[]
        maxInfectTimeVecUniform=[]   
        for i in range(0,iterations): #do 100 simulations on each graph
            #t is a list that records running time of simulation, S the number of susceptible at each t, I the number of infected at each t, and R the number of recovered at each t
            #a new time is added to the list t each time a transition between S, I and R occurs
            
            t, S, I, R = EoN.fast_SIR(G, beta, gamma,rho=None ,transmission_weight='weightT')
            
            while max(I)/n < 0.1: #continue running new simulations until an epidemic occurs (more than 10% in this case)
                t, S, I, R = EoN.fast_SIR(G, beta, gamma,rho=None,transmission_weight='weightT')
            maxInfectedVecTango.append(max(I)/n) #update the max proportion infected
            maxTimeVecTango.append(max(t)/n) #update simulation time relative to number of vertices
            maxInfectTimeVecTango.append(t[np.argmax(I)]/max(t)) #update proportion of time until max infected is reach (as proportion of total simulation time)

            t, S, I, R = EoN.fast_SIR(G, beta, gamma, rho=None,transmission_weight='weightU')
            
            while max(I)/n < 0.1:
                t, S, I, R = EoN.fast_SIR(G, beta, gamma, rho=None,transmission_weight='weightU')
            maxInfectedVecUniform.append(max(I)/n)
            maxTimeVecUniform.append(max(t)/n)
            maxInfectTimeVecUniform.append(t[np.argmax(I)]/max(t))
            
        #take averages over all simulations
        maxInfectedVecTangoAverage.append(np.mean(maxInfectedVecTango))
        maxInfectedVecUniformAverage.append(np.mean(maxInfectedVecUniform))
        maxTimeVecTangoAverage.append(np.mean(maxTimeVecTango))
        maxTimeVecUniformAverage.append(np.mean(maxTimeVecUniform))
        maxInfectTimeVecTangoAverage.append(np.mean(maxInfectTimeVecTango))
        maxInfectTimeVecUniformAverage.append(np.mean(maxInfectTimeVecUniform))
        
        #take approx 95% confidence intervals of these averages
        maxInfectedVecTangoStddev.append(1.96*np.std(maxInfectedVecTango)/np.sqrt(n))
        maxInfectedVecUniformStddev.append(1.96*np.std(maxInfectedVecUniform)/np.sqrt(n))
        maxTimeVecTangoStddev.append(1.96*np.std(maxTimeVecTango)/np.sqrt(n))
        maxTimeVecUniformStddev.append(1.96*np.std(maxTimeVecUniform)/np.sqrt(n))
        maxInfectTimeVecTangoStddev.append(1.96*np.std(maxInfectTimeVecTango)/np.sqrt(n))
        maxInfectTimeVecUniformStddev.append(1.96*np.std(maxInfectTimeVecUniform)/np.sqrt(n))
        
#plotting

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

plt.rc('font', **font)

plt.figure(figsize=(14, 10))
plt.errorbar(secLargestEigTango,maxInfectedVecTangoAverage,maxInfectedVecTangoStddev,None,'r.',label='Degree-based')
plt.errorbar(secLargestEigUniform,maxInfectedVecUniformAverage,maxInfectedVecUniformStddev,None,'b.',label='Uniform')
plt.xlabel('$\Lambda$')
plt.ylabel('Average Maximum Proportion Infected')
plt.title("Maximum proportion infected against $\Lambda$ for wildlife contact networks")
plt.legend()
plt.savefig('Epidemic_Against_EigAnimalNetworksgamma1over10niterations100error.pdf')

plt.figure(figsize=(14, 10))
plt.errorbar(secLargestEigTango,maxInfectTimeVecTangoAverage,maxInfectTimeVecTangoStddev,None,'r.',label='Degree-based')
plt.errorbar(secLargestEigUniform,maxInfectTimeVecUniformAverage,maxInfectTimeVecUniformStddev,None,'b.',label='Uniform')
plt.xlabel('$\Lambda$')
plt.ylabel('Average Proportion of Time Until Maximum Infected')
plt.title("Proportion of time until maximum infected against $\Lambda$ for wildlife contact networks")
plt.legend()
plt.savefig('Time_Against_EigAnimalNetworksgamma1over10niterations100error.pdf')

plt.figure(figsize=(14, 10))
plt.plot(secLargestEigTango,maxInfectedVecTangoAverage,'r.',label='Degree-based')
plt.plot(secLargestEigUniform,maxInfectedVecUniformAverage,'b.',label='Uniform')
plt.xlabel('$\Lambda$')
plt.ylabel('Average Maximum Proportion Infected')
plt.title("Maximum proportion infected against $\Lambda$ for wildlife contact networks")
plt.legend()
plt.savefig('Epidemic_Against_EigAnimalNetworksgamma1over10niterations100.pdf')

plt.figure(figsize=(14, 10))
plt.plot(secLargestEigTango,maxInfectTimeVecTangoAverage,'r.',label='Degree-based')
plt.plot(secLargestEigUniform,maxInfectTimeVecUniformAverage,'b.',label='Uniform')
plt.xlabel('$\Lambda$')
plt.ylabel('Average Proportion of Time Until Maximum Infected')
plt.title("Proportion of time until maximum infected against $\Lambda$ for wildlife contact networks")
plt.legend()
plt.savefig('Time_Against_EigAnimalNetworksgamma1over10niterations100.pdf')

plt.figure(figsize=(14, 10))
plt.errorbar(secLargestEigTango,maxTimeVecTangoAverage,maxTimeVecTangoStddev,None,'r.',label='Degree-based')
plt.errorbar(secLargestEigUniform,maxTimeVecUniformAverage,maxTimeVecUniformStddev,None,'b.',label='Uniform')
plt.xlabel('$\Lambda$')
plt.ylabel('Convergence Time')
plt.title("$\lambda_2$ against Convergence Time with $\gamma=0.001$")
plt.legend()
plt.savefig('ConvTime_Against_EigAnimalNetworksgamma1over10niterations100error.pdf')
