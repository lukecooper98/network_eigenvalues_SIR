import numpy as np
import sympy as sp
import scipy as sc
import networkx as nx
import matplotlib.pyplot as plt
import random
import glob
import errno
import EoN as EoN
from collections import Counter

#get the un-normalised degree-based matrix of a graph G 
def randWalkMatrix(G):
  
  G=nx.convert_node_labels_to_integers(G, first_label=0)
  
  A=np.zeros((len(list(G.nodes())),len(list(G.nodes()))))
 
  rowsum=[0 for i in list(G.nodes())]
  
  for e in list(G.edges()):
      A[e[0],e[1]]= 1/(G.degree(e[0])*G.degree(e[1]))
      A[e[1],e[0]]= A[e[0],e[1]]
      
      rowsum[e[0]] = rowsum[e[0]] + A[e[0],e[1]]
      rowsum[e[1]] = rowsum[e[1]] + A[e[1],e[0]]
  for i in range(0,len(list(G.nodes()))):
       A[i,i]= 1-rowsum[i]

  return A

#get the normalised degree-based matrix of a graph G
def normTangoMatrix(G):
  T=randWalkMatrix(G)
  offDiagSum=0
  rowsum=0
  for i in range(0,T.shape[0]):
    for j in range(0,T.shape[0]):
      if i != j:
        offDiagSum=offDiagSum+T[i,j]
  
  for i in range(0,T.shape[0]):
    for j in range(0,T.shape[0]):
      if i != j:
        T[i,j]= T[i,j]/offDiagSum
        rowsum = rowsum + T[i,j]
    T[i,i]= 1-rowsum
    rowsum=0
  return T

#get the uniform matrix of a graph G
def normUniformMatrix(G):
  G=nx.convert_node_labels_to_integers(G, first_label=0)
  A=np.zeros((len(list(G.nodes())),len(list(G.nodes()))))
  p=1/(2*len(list(G.edges)))
  rowsum=[0 for i in list(G.nodes())]
  
  for e in list(G.edges()):     
      A[e[0],e[1]]= p
      A[e[1],e[0]]= p
      
      rowsum[e[0]] = rowsum[e[0]] + p
      rowsum[e[1]] = rowsum[e[1]] + p
  for i in range(0,len(list(G.nodes()))):
       A[i,i]= 1-rowsum[i]
  
  
  return A

#get the laplace matrix of a graph G
def laplaceMatrix(G):
  A=nx.to_numpy_matrix(G,weight=None)
  D=np.identity(G.order())
  for i in range(0,G.order()):
    D[i,i]=G.degree[i]
  L=D-A
  return L

#get the unsigned laplace matrix of a graph G
def unsignedLaplaceMatrix(G):
  A=nx.to_numpy_matrix(G,weight=None)
  D=np.identity(G.order())
  for i in range(0,G.order()):
    D[i,i]=G.degree[i]
  L=D+A
  return L

#get the second largest eigenvalue of a matrix A
def secLargestEig(A):
  eigWalk=np.linalg.eig(A)
  lambda1=0
  for i in range(0,len(eigWalk[0])):
      if abs(eigWalk[0][i]) >= abs(lambda1):
        lambda1 = eigWalk[0][i]
        maxIndex1 = i
  
  lambda2=0
  for i in range(0,len(eigWalk[0])):
    if eigWalk[0][i] != lambda1:
      if abs(eigWalk[0][i]) >= abs(lambda2):
        lambda2 = eigWalk[0][i]
        maxIndex2 = i
  
  return lambda2

#get the second most positive eigenvalue of a matrix A (returns same result as secLargestEig if all eigenvalues of A are positive)
def secLargestEigPos(A):
  eigWalk=np.linalg.eig(A)
  lambda1=0
  for i in range(0,len(eigWalk[0])):
      if eigWalk[0][i] >= lambda1:
        lambda1 = eigWalk[0][i]
        maxIndex1 = i
  
  lambda2=0
  for i in range(0,len(eigWalk[0])):
    if eigWalk[0][i] != lambda1:
      if eigWalk[0][i] >= lambda2:
        lambda2 = eigWalk[0][i]
        maxIndex2 = i
  
  return lambda2

#return the largest eigenvalue of a matrix A
def largestEig(A):
  eigWalk=np.linalg.eig(A)
  lambda1=0
  for i in range(0,len(eigWalk[0])):
      if eigWalk[0][i] >= lambda1:
        lambda1 = eigWalk[0][i]
        maxIndex1 = i
  return lambda1

#return the second smallest eigenvalue of a matrix L (i.e. the fiedler value of L)
def fiedler(L):
  eigWalk=np.linalg.eig(L)
  lambda2=eigWalk[0][0]
  j=0
  while lambda2==0:
    j=j+1
    lambda2=eigWalk[0][j]
  for i in range(0,len(eigWalk[0])):
    if eigWalk[0][i] <= lambda2 and eigWalk[0][i] != 0 :
      lambda2 = eigWalk[0][i]
  
  return lambda2

#returns the connected components of a graph G as a list of graphs
def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c)
        
#reads in a list from a text file
def listReader(fileName):
  arr_intValues = []
  myFile = open(fileName, "r")
  for myLine in myFile:
    arr_intValues.append(float(myLine))
  return arr_intValues


