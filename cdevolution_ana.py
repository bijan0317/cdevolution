#Title:The cooperation-defection evolution on social networks [Analytical version]
#Author: Bijan Sarkar
#Reference: https://sourceforge.net/projects/pycx/files/   
import matplotlib
matplotlib.use('TkAgg')
from pylab import *

n = 100     # pattern of grid: n * n
Dh = 1./n   # spatial resolution, the space size is [0,1] * [0,1]
Dt = 0.002  # temporal resolution, another value 0.02 
# payoff elements
Pi_CC=2.2      
Pi_CD=0.3      
Pi_DC=1.5      
Pi_DD=0.5      
# network parameters
Gamma_max=1.145
z=5       # mean=average node dgree; another value z=50
v=0       # variance on the node dgree distribution
e_scale=0.01 # the factor is equal to (1/N), see the main article; another value e_scale=0.001

# birth-death part for cooperators
a_3=((Pi_CC-Pi_CD)/Gamma_max)*((v+z**2-z)/z)
a_2=((v+z**2-z)/z)*((Pi_CD/Gamma_max)-1)+(Pi_CC-Pi_CD)/Gamma_max 
a_1=(Pi_CD/Gamma_max)-1         

# birth-death part for defectors
b_3=((Pi_DD-Pi_DC)/Gamma_max)*((v+z**2-z)/z)
b_2=((v+z**2-z)/z)*((Pi_DC/Gamma_max)-1)+(Pi_DD-Pi_DC)/Gamma_max 
b_1=(Pi_DC/Gamma_max)-1        
        
        
def initialize():
    global x, y, nextx, nexty, cdata, ddata, ccount, dcount, ecount
    
    ccount=0    # counts the point at which cooperator group dominates
    dcount=0    # counts the point at which defector group dominates
    ecount=0    # counts the point at which both concentrations are same 
    
    cdata=[]    # measures the concentration x(p,q) over time variation
    ddata=[]    # measures the concentration y(p,q) over time variation
    x = zeros([n, n])
    y = zeros([n, n])
    for p in range(n):
        for q in range(n):
            x[p, q] = uniform(0, 0.3)   # initial concentration of cooperators
            y[p, q] = uniform(0, 0.3)   # initial concentration of defectors
            x[30,60]=0.3032             # <- the step at which the input values change
            y[30,60]=0.2728             # <- the step at which the input values change
            if x[p,q]>y[p,q]:
                 ccount +=1
            elif x[p,q]<y[p,q]:
                 dcount +=1
            elif x[p,q]==y[p,q]:
                 ecount +=1
    
    print(x[30,60],y[30,60],ccount,dcount,ecount)  # <- the step at which the input values change
                
    nextx = zeros([n, n])
    nexty = zeros([n, n])
    
    nextx =  x
    nexty =  y
    
def observe():
    global x, y, nextx, nexty, cdata, ddata, ccount, dcount, ecount
    clf()
    
    subplot(2, 2, 1)
    cla()
    imshow(x, vmin = 0, vmax = 1,cmap= 'Blues')
    cbar0=colorbar(orientation='vertical')
    title('Cooperator',fontsize=9)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    cbar0.ax.tick_params(labelsize=8)
    
    subplot(2, 2, 2)
    cla()
    imshow(y, vmin = 0, vmax = 1, cmap = 'Reds')
    cbar=colorbar(orientation='vertical')
    title('Defector',fontsize=9)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    cbar.ax.tick_params(labelsize=8)
    
    subplot(2, 1, 2)
    cla()
    cdata.append(x[30,60])  # <- the step at which the input values change
    ddata.append(y[30,60])  # <- the step at which the input values change
    grid()
    plot(cdata, label = 'Cooperator')
    plot(ddata, label = 'Defector')
    legend(loc='best', prop={'size': 9})
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    
def update():
    global x, y, nextx, nexty, cdata, ddata, ccount, dcount, ecount
    
    for p in range(n):
        for q in range(n):
            xC, xR, xL, xU, xD = x[p,q], x[(p+1)%n,q], x[(p-1)%n,q], \
                                 x[p,(q+1)%n], x[p,(q-1)%n]
            yC, yR, yL, yU, yD = y[p,q], y[(p+1)%n,q], y[(p-1)%n,q], \
                                 y[p,(q+1)%n], y[p,(q-1)%n]
            xrdo = 10**(-5)*((xR + xL - xU - xD) / (Dh**2))   # xrdo=sig_x x (random drift operator)
            yrdo = 10**(-5)*((yR + yL - yU - yD) / (Dh**2))   # yrdo=sig_y x (random drift operator)
            #xrdo=0  # <- the step at which the input values change
            #yrdo=0  # <- the step at which the input values change
            if xC+yC<1:
                if 0<xC<1 and 0<yC<1:
                   nextx[p,q] = xC+(e_scale*a_3*(xC/(xC+yC))**3 +e_scale*a_2*(xC/(xC+yC))**2+e_scale*a_1*(xC/(xC+yC))+xrdo)*Dt          
                   nexty[p,q] = yC+(e_scale*b_3*(yC/(xC+yC))**3 +e_scale*b_2*(yC/(xC+yC))**2+e_scale*b_1*(yC/(xC+yC))+yrdo)*Dt
                
                elif xC==0 and 0<yC<1:
                   nexty[p,q] = yC+(e_scale*b_3*(yC/(xC+yC))**3 +e_scale*b_2*(yC/(xC+yC))**2+e_scale*b_1*(yC/(xC+yC))+yrdo)*Dt
                
                elif 0<xC<1 and yC==0:
                   nextx[p,q] = xC+(e_scale*a_3*(xC/(xC+yC))**3 +e_scale*a_2*(xC/(xC+yC))**2+e_scale*a_1*(xC/(xC+yC))+xrdo)*Dt 
                
                elif xC==0 and yC==0:
                   nextx[p,q] = xC
                   nexty[p,q] = yC     
                        
            if nextx[p,q]<0:
               nextx[p,q]=0
            elif nextx[p,q]>1:
               nextx[p,q]=1
               nexty[p,q]=0
            if nexty[p,q]<0:
               nexty[p,q]=0
            elif nexty[p,q]>1:
               nexty[p,q]=1
               nextx[p,q]=0
            if xC+yC>1: 
               nextx[p,q] = xC
               nexty[p,q] = yC                              
   
    x, nextx = nextx, x
    y, nexty = nexty, y
    
    ccount=0
    dcount=0      
    ecount=0
    ecountr=0     # the refine version of ecount
    coexist=0
    extinction=0
    
    ccount = (x>y).sum()
    dcount = (x<y).sum()
    ecount = (x==y).sum()
    ecountr = ((y!=0) & (x==y)).sum()
    coexist=((y!=0) & (x!=0)).sum()
    extinction=((y==0) & (x==0)).sum()
   
    cmx=0      # the counter counts how many x=1 are noticeable in the system
    dmx=0      # the counter counts how many y=1 are noticeable in the system
    cmx = (x == 1).sum()
    dmx = (y == 1).sum()

    cmi=0      # the counter counts how many x=0 are noticeable in the system
    dmi=0      # the counter counts how many y=0 are noticeable in the system
    cmi = (x == 0).sum()
    dmi = (y == 0).sum()

    print("%.5f" %x[30,60],"%.5f" %y[30,60],ccount,cmx,cmi,dcount,dmx,dmi,ecount,ecountr,coexist,extinction)  # <- the step at which the input values change

import pycxsimulator
pycxsimulator.GUI(stepSize = 100).start(func=[initialize, observe, update])        
