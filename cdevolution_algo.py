#Title:The cooperation-defection evolution on social networks [Algorithm version]
#Author: Bijan Sarkar
#Reference: https://sourceforge.net/projects/pycx/files/   
import matplotlib
matplotlib.use('TkAgg')
from pylab import *
import networkx as nx
import random as rd 
import sys
import timeit

start = timeit.default_timer()
 
N = 100 # a graph of N nodes; another value N=1000
#m = 5  # number of edges per new node; m=2,3,4,5
Pi_CC=1.5   
Pi_CD=-0.3
Pi_DC=1.8   
Pi_DD=0     
Gamma_max=0.77
p=0.5   # probability of attachment of an individual
q=0.4   # probability of leaving of an individual

sig_x=0.3*10**(-1) # coefficient of random drift
sig_y=0.3*10**(-1) # coefficient of random drift
#sig_x=0
#sig_y=0

def initialize():
    global g,nextg,z,v,countc,countd,countv,cdata,ddata
    cdata=[]
    ddata=[]
    g = nx.random_regular_graph(5,100)  #g = random_regular_graph(d, n, seed=None) 
    #g = nx.barabasi_albert_graph(n,m)
    g.pos = nx.spring_layout(g)
    g.count = 0  # number of rounds in each generation
    
    sum_of_edges = sum(g.degree(l) for l in g.nodes())  
    z= sum_of_edges//100    # mean=z
    print(z)
    
    v=sum((g.degree(l)-z)**2 for l in g.nodes())   # variance=v
    print(v)
    
    for i in g.nodes_iter():
      g.node[i]['state'] = 1 if random()<0.3 \
                                 else 0 if 0.3<random()<0.7 else 0.5
             
    countc= 0   # the counter counts number of cooperators (1s)
    countd= 0   # the counter counts number of defectors (0s)
    countv= 0   # the counter counts number of vacant places (0.5s)
    
    nod=[g.node[i]['state'] for i in g.nodes_iter()]
    countc=nod.count(1)
    countd=nod.count(0)
    
    countv=N-countc-countd 
    print(countc,countd,countv)
    
    nextg=g
    
def observe():
    global g,nextg,z,v,countc,countd,countv,cdata,ddata
    
    subplot(2, 1, 1)
    node_color = []
    for i in g.nodes_iter():
       if g.node[i]['state'] == 1:
        node_color.append('blue')
       elif g.node[i]['state'] == 0:
        node_color.append('red')
       elif g.node[i]['state'] == 0.5:
        node_color.append('white')   

    cla()
    nx.draw(g, with_labels=False, node_size=70, node_color=node_color,
            pos = g.pos)
    title('Cooperator-defector distribution',fontsize=9)
    
    subplot(2, 1, 2)
    cla()
    cdata.append(countc)
    ddata.append(countd)
    grid()
    plot(cdata, label = 'Cooperator')
    plot(ddata, label = 'Defector')
    legend(loc='best', prop={'size': 9})
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)

def find_out(i):
    global g,pre_neigh
    for j in g.neighbors(i):
              if g.node[j]['state'] ==0.5:
                 if pre_neigh.count(j) == 0:
                    return j 

def find_out2(i):
    global g,pre_neigh
    for c in range(20):    # c is an arbitrary counter
       j = rd.choice(g.neighbors(i))
       if pre_neigh.count(j) == 0:
          return j                  
      
def update(): 

    global g,nextg,z,v,countc,countd,countv,cdata,ddata,pre_neigh
    
    g.count += 1

#  birth-death effect in the cooperation-defection evolution    
    if g.count % 15 == 0: # birth-death occurs atmost once in each generation of 15 time steps
       pre_neigh = []
       for i in g.nodes_iter():
          if pre_neigh.count(i) == 0:  # 1st checking:if a node does not participate previously, then the node is allowable
             if g.degree(i)>0:   
                icountc=0   # the counter counts cooperators (1s) in a neighbourhood
                icountd=0   # the counter counts defectors (0s) in a neighbourhood
                
                nei=[g.node[j]['state'] for j in g.neighbors(i)]
                icountc=nei.count(1)
                icountd=nei.count(0) 
            
                j = find_out(i)
                if j:
                  nbd=j
                else:
                  nbd= find_out2(i)
               
                if nbd:  # 2nd checking: if nbd exists, then proceed
                   if g.node[i]['state']==1:
                      if random()<(icountc*Pi_CC+icountd*Pi_CD)/(z*Gamma_max): 
                        nextg.node[nbd]['state']=1
                        pre_neigh.append(nbd)
                      elif random()<icountc/z:
                        nextg.node[i]['state']=0.5   
                   elif g.node[i]['state']==0:
                      if random()<(icountc*Pi_DC+icountd*Pi_DD)/(z*Gamma_max):
                         nextg.node[nbd]['state']=0
                         pre_neigh.append(nbd)
                      elif random()<icountd/z:
                         nextg.node[i]['state']=0.5 
       g, nextg = nextg, g
 
# random drift effect in the cooperation-defection evolution 
    pre_neigh=[]  
    for i in g.nodes_iter():
        ncountc=0       # the counter counts cooperators in a neighbourhood with including focal itself
        ncountd=0       # the counter counts defectors in a neighbourhood with including focal itself
        if pre_neigh.count(i) == 0:  # 1st checking : the node is not allowable if it  participates  previously
           if g.node[i]['state']==1:
              ncountc +=1
           elif g.node[i]['state']==0:
              ncountd +=1
           nei=[g.node[j]['state'] for j in g.neighbors(i)]
           ncountc=nei.count(1)+ncountc
           ncountd=nei.count(0)+ncountd  
                            
           ratio1=ncountc/z
           ratio2=ncountd/z
           if 0<ratio2<ratio1:
              if random()<p*sig_x:      
                  j = find_out(i)
                  if j:
                       nextg.node[j]['state']=1
                       pre_neigh.append(j)          
                  else:
                    nextg.node[i]['state']=1     
                    pre_neigh.append(i)          #  i-th node may be vacant
              elif random()<q*sig_x: 
                 nextg.node[i]['state']=0.5          
           elif 0<ratio1<ratio2:
               if random()<p*sig_y: 
                    j = find_out(i)
                    if j:
                         nextg.node[j]['state']=0
                         pre_neigh.append(j)          
                    else:
                      nextg.node[i]['state']=0
                      pre_neigh.append(i)
               elif random()<q*sig_y:            
                  nextg.node[i]['state']=0.5         
    
    g, nextg = nextg, g  
            
    countc=0
    countd=0    
    
    nod=[g.node[i]['state'] for i in g.nodes_iter()]
    countc=nod.count(1)
    countd=nod.count(0)
    
    countv=N-countc-countd 

    print(countc,countd,countv)

#simulation of node movement
    g.pos = nx.spring_layout(g, pos = g.pos, iterations = 5)      
# --althernative of pycxsimulator-----
#t=0
#initialize()
#input("Press <Run[F5]> to continue")
#for t in range(1,20000):
#   update() 
#stop = timeit.default_timer()
#print('Time: ', stop - start)
#-----------------
import pycxsimulator
pycxsimulator.GUI(stepSize = 100).start(func=[initialize, observe, update])
