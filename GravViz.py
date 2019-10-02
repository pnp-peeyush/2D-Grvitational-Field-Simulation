import numpy as np
import matplotlib.pyplot as plt
import time

n = 25
m=15
e=0.5                #elasticity of collision (-1,1) factor by which momentum is lost after collision
G=15
t=34.0
dt=0.005
itn=int(t/dt)
skip=10
Y, X = np.mgrid[-m:m:0.10, -n:n:0.10]
#r=np.array([50.0,20.0])
#Xo,Yo=np.array([-50.0,50.0]),np.array([0.0,0.0])
#print(Yo)

def field(Xo,Yo,Mo,ro):
    R=(np.sqrt(((X-Xo)**2)+((Y-Yo)**2)))
    Rx=(X-Xo)/R
    Ry=(Y-Yo)/R
    Rx=np.where(np.isnan(Rx),0,Rx)
    Ry=np.where(np.isnan(Ry),0,Ry)
    g = -10*Mo/(R**2)
    g=np.where(np.isnan(g),0,g)
    g=np.where(g<-G*Mo/(ro**2),-G*Mo*R/(ro**3),g)
    U, V = g * Rx, g * Ry
    
    return U,V,g
                                                   #(% vel per sec) 
      #X loc  #Y loc   #Mass    #Radius  #Velocity    #decay        #color          (IFD==Initial Field Data)
IFD=[(-20.0,   00.0,   0150.0,  02.0,    ( 0,  3.75), -0.01,      '#9ccc3c'),
     ( 20.0,   00.0,   0150.0,  02.0,    ( 0, -3.75), -0.01,      '#80a9ff'),]
     #( 00.0,   00.0,   0150.0,  05.5,    ( 0,  0.00), -0.25,      '#f7347a')]
     
Simul_pos=[]    
All_min_dist=[]
t_stamp=[0.0]
Simul_field=[]
    
for i in range (0,len(IFD)):
    min_dist=[]
    for j in range(0,len(IFD)):
        if i!=j:
            min_dist.append(IFD[i][3]+IFD[j][3]) 
    All_min_dist.append(min_dist)
    Simul_pos.append([[IFD[i][0]],[IFD[i][1]]])


''''###############################################################################'''
All_pos=np.empty((0,2))
All_vel=np.empty((0,2))
All_mass=np.empty((0,2))
All_decay=np.empty((0,2))
for i in range(0,len(IFD)):
    All_pos=np.vstack((All_pos,np.array([IFD[i][0],IFD[i][1]])))
    All_vel=np.vstack((All_vel,np.array([IFD[i][4][0],IFD[i][4][1]])))
    All_mass=np.vstack((All_mass,np.array([IFD[i][2],IFD[i][2]])))
    All_decay=np.vstack((All_decay,np.array([1+(0.01*IFD[i][5]),1+(0.01*IFD[i][5])])))

''''###############################################################################'''
T1=time.time()
for s in range(0,itn):
    t1=time.time()
    All_fields=[]
    All_dist=[]
    All_dist_mag=[]
    for i in range (0,len(IFD)):
        #All_fields.append(field(All_pos[i][0],All_pos[i][1],IFD[i][2],IFD[i][3]))
        disti=[]
        r1=All_pos[i]
        
        for j in range(i,len(IFD)):
            r2=All_pos[j]
            if i!=j:
                ri=r2-r1
                disti.append(ri)
            elif i==j==len(IFD)-1:
                pass
            
        All_dist.append(disti)
        
    for i in range(0,len(IFD)):
        for j in range(0,i):
            All_dist[i]=[-1*All_dist[i-1-j][i-1]]+All_dist[i]
    
    if s%skip==0:        
        for i in range (0,len(IFD)):
            All_fields.append(field(All_pos[i][0],All_pos[i][1],IFD[i][2],IFD[i][3]))  
                
        final_field=[]
        U_f=0
        V_f=0
        g_f=0
        for i in range(0,len(All_fields)):
            U_f=U_f+All_fields[i][0]
            V_f=V_f+All_fields[i][1]
            g_f=g_f+All_fields[i][2]
            
        Simul_field.append([U_f,V_f,g_f])
        
    '''###############################################################################'''
    
    All_forces=np.empty((0,2))
    
    for i in range (0,len(IFD)):
        k=0
        fx=0
        fy=0
        dist_mag=[]
        for j in range(0,len(IFD)):
            if i!=j:
                r=np.hypot(All_dist[i][k][0],All_dist[i][k][1])
                f=(G*IFD[i][2]*IFD[j][2])/(r*r)
                fx=fx+((f*All_dist[i][k][0])/r)    
                fy=fy+((f*All_dist[i][k][1])/r)
                k=k+1
        All_forces=np.vstack((All_forces,np.array([fx,fy])))
        
    '''###############################################################################'''
    
    All_accl=All_forces/All_mass
        
    '''###############################################################################'''
    
    All_vel=All_decay*(All_vel+(All_accl*dt))

    for i in range (0,len(IFD)):
        for j in range(0,len(All_dist[i])):
            #print i,j
            if np.hypot(All_dist[i][j][0],All_dist[i][j][1])<=All_min_dist[i][j]:
                All_vel[i]=e*All_vel[i]
                if j<i:
                    All_vel[j]=e*All_vel[j]
                elif j>=i:
                    All_vel[j+1]=e*All_vel[j+1]
            #pass
        
    '''###############################################################################'''
    
    for i in range (0,len(IFD)):
        Simul_pos[i][0].append(All_pos[i][0])
        Simul_pos[i][1].append(All_pos[i][1])
    
    All_pos=(All_pos+(All_vel*dt))
    '''###############################################################################'''  
    t_stamp.append((s+1)*dt)  
    
    '''###############################################################################'''
    if s%100==0:
        #print All_dist
        t2=time.time()
        print(t2-t1)

T2=time.time()
print(T2-T1)
plt.figure(figsize=(13,9.8))
plt.axes([0.025, 0.025, 0.95, 0.95])
plt.axis('scaled')

plt.imshow(g_f,extent=[X.min(), X.max(), Y.min(), Y.max()],interpolation='none',cmap='viridis_r',origin='lower')
#plt.colorbar()
#plt.quiver(X, Y, U_f, V_f, g_f,linewidth=.5,scale=900,cmap='viridis_r')
plt.streamplot(X, Y, V_f, -U_f, color=g_f, linewidth=1, cmap='viridis_r',density=[1,1])
#plt.colorbar()

plt.xlim(-n, n)
plt.ylim(-m, m)
#plt.xticks(())
#plt.yticks(())
#plt.grid()
#plt.savefig("E:/BouncingBall/GravSim/render/fig"+str(s)+".png")

plt.show()
