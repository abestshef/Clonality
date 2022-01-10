### Import the necessary libraries

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from random import *

### Global constants

Rec=0.5 # Staph EC growth rate
Ric=0.25 # Staph IC growth rate
K=50000 # Staph EC carrying capacity
BETA=0.25 # Macrophage phagocytosis rate
CB=100 # Half saturation constant
MU=1 # Macrophage initial killing rate

REPS=20 #No. of replicate runs
S0=10000 # Initial no. of staph
M0=2000 # Initial no. of macrophages
MAX_TIME=25 # Max time

plt.rcParams.update({'font.size': 16})

bluestore=[]
redstore=[]
mastore=[]
bluemacs=[]
redmacs=[]
dualmacs=[]
masize=np.zeros((REPS,10))

# Main function
for reps in range(0,REPS):

    # Some local constants / lists
    tsteps=[0]
    SecBLUE=[int(S0/2)]
    SecRED=[int(S0/2)]
    SicBLUE=[0]
    SicRED=[0]
    MA=np.zeros((M0,2))
    current_t=0  
    muall=[MU for i in range(M0)]
    ros=np.zeros(M0)
    print('Replicate %i' % reps)
    
    # Main run
    while current_t < MAX_TIME:
        
        # After 2.5 hours, all EC bacteria killed
        if current_t > 2.5:
            SecBLUE[-1]=0
            SecRED[-1]=0
            if tsteps[-2]<2.5:
                for sizes in range(10):
                    masize[reps,sizes]=(len(np.where((MA[:,0]+MA[:,1])==sizes)[0]))
        
        # Calculate event rates
        event1 = Rec*(SecBLUE[-1]+SecRED[-1])*(1-(SecBLUE[-1]+SecRED[-1])/K) # EC staph replication
        event2 = np.sum([muall[i]*(MA[i,0]+MA[i,1]) for i in range(M0)]) # Staph killed
        event4 = BETA*(SecBLUE[-1]+SecRED[-1])/((SecBLUE[-1]+SecRED[-1])+CB)*(M0) # Phagocytosis
        event3 = Ric*(SicBLUE[-1]+SicRED[-1])

        # Find time to next event
        scale = event1+event2+event3+event4     
        if scale==0:
            break     
        dt = -np.log(np.random.uniform()) / scale
        current_t=tsteps[-1]+dt

        # Find event
        eventcheck = np.random.uniform()
        if eventcheck < event1/scale: #Event is staph replication  

            Seccheck=np.random.uniform()
            if Seccheck < SecBLUE[-1]/(SecBLUE[-1]+SecRED[-1]):
                SecBLUE.append(SecBLUE[-1]+1)
                SecRED.append(SecRED[-1])
            else:
                SecBLUE.append(SecBLUE[-1])
                SecRED.append(SecRED[-1]+1)
            
            SicBLUE.append(SicBLUE[-1])
            SicRED.append(SicRED[-1])
            
        elif eventcheck<(event1+event2)/scale: # Event is staph killed
            
            scheck=np.random.uniform()
            macheck=np.random.uniform()
            
            if scheck < SicBLUE[-1]/(SicBLUE[-1]+SicRED[-1]):
                colour=0
                sic=SicBLUE[-1]
            else:
                colour=1
                sic=SicRED[-1]
            
            for i in range(0,M0):
                if macheck<np.sum(MA[0:i+1,colour])/(sic):
                    MA[i,colour]-=1
                    if colour==0:
                        SicBLUE.append(SicBLUE[-1]-1)
                        SicRED.append(SicRED[-1])
                    else:
                        SicBLUE.append(SicBLUE[-1])
                        SicRED.append(SicRED[-1]-1)
                    break
            SecBLUE.append(SecBLUE[-1])
            SecRED.append(SecRED[-1])
            
        elif eventcheck<(event1+event2+event3)/scale: # Event is IC staph replicates
            
            scheck=np.random.uniform()
            macheck=np.random.uniform()
            
            if scheck < SicBLUE[-1]/(SicBLUE[-1]+SicRED[-1]):
                colour=0
                sic=SicBLUE[-1]
            else:
                colour=1
                sic=SicRED[-1]
            
            for i in range(0,M0):
                if macheck<np.sum(MA[0:i+1,colour])/(sic):
                    MA[i,colour]+=1
                    if colour==0:
                        SicBLUE.append(SicBLUE[-1]+1)
                        SicRED.append(SicRED[-1])
                    else:
                        SicBLUE.append(SicBLUE[-1])
                        SicRED.append(SicRED[-1]+1)
                    break
            SecBLUE.append(SecBLUE[-1])
            SecRED.append(SecRED[-1])

        else: #Event is phagocytosis
            
            Seccheck=np.random.uniform()
            choice=np.random.randint(M0)
            if ros[choice]==0:
                ros[choice]=current_t
            if Seccheck < SecBLUE[-1]/(SecBLUE[-1]+SecRED[-1]):
                SecBLUE.append(SecBLUE[-1]-1)
                SecRED.append(SecRED[-1])
                SicBLUE.append(SicBLUE[-1]+1)
                SicRED.append(SicRED[-1])
                MA[choice,0]+=1
            else:
                SecBLUE.append(SecBLUE[-1])
                SecRED.append(SecRED[-1]-1)
                SicBLUE.append(SicBLUE[-1])
                SicRED.append(SicRED[-1]+1)
                MA[choice,1]+=1
            
        # Update cell killing rates
        kchange=np.where(ros[:]>0)   
        if len(kchange[0])>0:
            for i in range(len(kchange[0])):
                muall[kchange[0][i]]=MU*(1-(current_t-ros[kchange[0][i]])**2/((current_t-ros[kchange[0][i]])**2+9**2))

        # Update lists
        tsteps.append(current_t)
        
    # Update variable stores
    bluestore.append(SicBLUE[-1])
    redstore.append(SicRED[-1])
    
    bm=0
    rm=0
    dm=0
    for i in range(M0):
        if MA[i,0]>0 and MA[i,1]==0:
            bm+=1
        elif MA[i,0]==0 and MA[i,1]>0:
            rm+=1
        elif MA[i,0]>0 and MA[i,1]>0:
            dm+=1
    bluemacs.append(bm)
    redmacs.append(rm)
    dualmacs.append(dm)
    
    if reps==0:
        mastore=MA
    else:
        if np.sum(MA)>np.sum(mastore):
            mastore=MA
        

# Values needed for plots
bmeans=np.mean(bluemacs)
bstd=np.std(bluemacs)
rmeans=np.mean(redmacs)
rstd=np.std(redmacs)
dmeans=np.mean(dualmacs)
dstd=np.std(dualmacs)
allmeans=[bmeans,rmeans,dmeans]
allstd=[bstd,rstd,dstd]
xax2=['Blue','Red','Dual']

# Build the main plot
fig, ax = plt.subplots()
ax.bar(xax2, allmeans, yerr=allstd, align='center', color=['blue','red','grey'], capsize=10)
ax.set_xticks(xax2)
plt.ylabel('Cells w/ IC bacteria at 25h')
plt.title('1:1')
plt.ylim([-0.5, 3])
plt.tight_layout()

# Second plot for each replicate
xax=np.linspace(1,REPS,REPS)
fig, ax = plt.subplots()
ax.bar(xax, bluemacs, 0.5, label='Blue',color='blue')
ax.bar(xax, redmacs, 0.5, bottom=bluemacs,
       label='Red',color='red')
ax.bar(xax, dualmacs, 0.5, bottom=redmacs,
       label='Dual',color='gray')
ax.set_ylabel('Macrophages w/ IC bacteria at 25h')
ax.set_xlabel('Simulation No.')
ax.set_xlim([0.5,20.5])
ax.set_ylim([0,5])
plt.legend()
plt.tight_layout()
