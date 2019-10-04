import numpy as np
import matplotlib.pyplot as plt
import time

n    = 20                                       # 2*n is the x dimension of simulated space in meters
m    = 20                                       # 2*m is the y dimension of simulated space in meters
e    = 0.99                                     # coefficient of restitution varies from (-1,1) 
G    = 5                                        # Universal gravitational constant in m^3/(kg.s^2)
t    = 40.0                                     # simulated time in seconds
dt   = 0.01                                     # time step (interval between each iteration of the master loop)
itn  = int(t/dt)                                # number of iterations of the master loop          
skip = 20                                       # number of iterations to skip for gravitational potential field calculation
gs   = 0.2                                      # grid size
Y, X = np.mgrid[-m:m:gs, -n:n:gs]               # dividing simulated space into finer grid for gravitational potential field calculation

'''function to calculate magnitude of a vector'''    
def vmag(a):                               
    return np.sqrt(np.einsum('i,i',a,a))

'''function to normalize a vector'''    
def vnorm(a):
    return a/vmag(a)

'''function to calculate gravitational potential field by one object in the entire simulated space'''    
def field(Xo,Yo,Mo,ro):
    R=(np.sqrt(((X-Xo)**2)+((Y-Yo)**2)))        # Magnitude Position vector field (position vector from object location to each grid point)
    Rx=(X-Xo)/R                                 # X component of the position vector field
    Ry=(Y-Yo)/R                                 # Y component of the position vector field
    Rx=np.where(np.isnan(Rx),0,Rx)              
    Ry=np.where(np.isnan(Ry),0,Ry)              
    g=R                                         # Gravitational potential field initailised to position field just to get an array of correct size
    g=np.where(R>=ro,-G*Mo/(np.square(R)),g)    # field outside the sphere
    g=np.where(R<ro,-G*Mo*R/(ro**3),g)          # field inside the sphere
    U, V = g * Rx, g * Ry                       # x and y component of the gravitational potential vector field
    
    return U,V,g
                                                                
'''IFD (Initial Field Data) is a list of tuples where each tuple holds inital values of all necessary variables required for simulation for each object'''

#                                                         (% velocity lost per sec)  
#                                                                     \             
#                                                                     |            
#                                                                    \/            
      #X loc  #Y loc   #Mass(kg)    #Radius(m)  #Velocity(m/s)    #decay         #color          
IFD=[(-10.0,   00.0,   0165.0,      03.25,       ( 0.0,  3.2),      -0.005,        '#9ccc3c'),            #object 1 data
     ( 10.0,   00.0,   0135.0,      02.5,       ( 0.1, -3.3),      -0.005,        '#80a9ff'),]           #object 2 data
    #( 00.0,   00.0,   0020.0,      05.5,       ( 0.0,  2.00),     -0.001,        '#f7347a')]            ...and so on

'''###############################################################################'''     
               
Simul_pos=[]                                    #stores position of each object for each simlated time instance  
All_min_dist=[]                                 #stores minimum distance of between all possible pairs of objects (basically a threshold for collision)
Simul_field=[]                                  #stores potential field data of the system of objects for each simulated time instance

'''Storing initial values of Simul_pos and All_min_dist'''    
for i in range (0,len(IFD)):
    min_dist=[]
    for j in range(0,len(IFD)):
        if i!=j:
            min_dist.append(IFD[i][3]+IFD[j][3]) 
    All_min_dist.append(min_dist)
    Simul_pos.append([[IFD[i][0]],[IFD[i][1]]])

'''###############################################################################'''

'''
Initializing position, velocity, mass, and decay parameter for each objet, position and velocity arrays store values for single iteration of the master loop and 
are initialized in each iteration
'''
All_pos=np.empty((0,2))
All_vel=np.empty((0,2))
All_mass=np.empty((0,2))
All_decay=np.empty((0,2))
for i in range(0,len(IFD)):
    All_pos=np.vstack((All_pos,np.array([IFD[i][0],IFD[i][1]])))
    All_vel=np.vstack((All_vel,np.array([IFD[i][4][0],IFD[i][4][1]])))
    All_mass=np.vstack((All_mass,np.array([IFD[i][2],IFD[i][2]])))
    All_decay=np.vstack((All_decay,np.array([1+(0.01*IFD[i][5]),1+(0.01*IFD[i][5])])))

T1=time.time()
'''###############################################################################

   ................................THE MASTER LOOP................................

   ###############################################################################'''

for s in range(0,itn):
    t1=time.time()
    '''
    Calculation of position vector from a one object to every other object for all the objects (english is not my 1st language >_<)
    '''
    All_dist=[]                                 #stores position vectors
    for i in range (0,len(IFD)):
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
            
    '''###############################################################################'''
    
    '''
    Calculation of field data for the entire system of objects. Field data is calculated only for visualization purposes hence need not be calculated at each
    iteration of master loop, 'skip' value allows to reduce the number of times field is calculated which in turn reduces the time of execution of master loop
    pro tip: want to reduce time of calculation? then don't do calculation ;)
    '''
    
    All_fields=[]                               # stores field data for individual objects
    if s%skip==0:        
        for i in range (0,len(IFD)):
            All_fields.append(field(All_pos[i][0],All_pos[i][1],IFD[i][2],IFD[i][3]))  
                
        U_f=0
        V_f=0
        g_f=0
        
        '''Calculation of effective field data due to the system of objects'''
        
        for i in range(0,len(All_fields)):
            U_f=U_f+All_fields[i][0]
            V_f=V_f+All_fields[i][1]
            g_f=g_f+All_fields[i][2]
            
        Simul_field.append([U_f,V_f,g_f])
        
    '''###############################################################################'''
    
    '''Calculation of force acting on each object due to every other object'''
    
    All_forces=np.empty((0,2))                  #stores the exact thing mentioned above
    
    for i in range (0,len(IFD)):
        k=0
        fx=0
        fy=0
        for j in range(0,len(IFD)):
            if i!=j:
                r=vmag(All_dist[i][k])
                if r>=All_min_dist[i][k]:
                    f=(G*IFD[i][2]*IFD[j][2])/(r*r)
                    fx=fx+((f*All_dist[i][k][0])/r)    
                    fy=fy+((f*All_dist[i][k][1])/r)
                else:
                    f=(G*IFD[i][2]*IFD[j][2])/((IFD[i][3])**3)
                    fx=fx+(f*All_dist[i][k][0])    
                    fy=fy+(f*All_dist[i][k][1])
                k=k+1
        All_forces=np.vstack((All_forces,np.array([fx,fy])))
        
    '''###############################################################################'''
    
    '''Calculation of accleration of each object'''
    
    All_accl=All_forces/All_mass
        
    '''###############################################################################'''
    
    '''Calculation of accleration of each object'''
    
    All_vel=All_decay*(All_vel+(All_accl*dt))

    '''Check for collision and do post collision velocity updates'''

    for i in range (0,len(IFD)):
        for j in range(0,len(All_dist[i])):
            if vmag(All_dist[i][j])<=All_min_dist[i][j]:
                if j<i:
                    All_vel[i]=(((All_mass[i]+(e*All_mass[j]))*All_vel[i])+((1-e)*All_mass[j]*All_vel[j]))/(All_mass[i][0]+All_mass[j][0])
                    All_vel[j]=(((All_mass[j]+(e*All_mass[i]))*All_vel[j])+((1-e)*All_mass[i]*All_vel[i]))/(All_mass[i][0]+All_mass[j][0])
                elif j>=i:
                    All_vel[i]=(((All_mass[i]+(e*All_mass[j+1]))*All_vel[i])+((1-e)*All_mass[j+1]*All_vel[j+1]))/(All_mass[i][0]+All_mass[j+1][0])
                    All_vel[j+1]=(((All_mass[j+1]+(e*All_mass[i]))*All_vel[j+1])+((1-e)*All_mass[i]*All_vel[i]))/(All_mass[i][0]+All_mass[j+1][0])
        
    '''###############################################################################'''
    
    '''store position of each object for each simlated time instance'''
    
    for i in range (0,len(IFD)):
        Simul_pos[i][0].append(All_pos[i][0])
        Simul_pos[i][1].append(All_pos[i][1])
    
    '''Update position of all objects'''
    
    All_pos=(All_pos+(All_vel*dt))
    
    '''###############################################################################'''  
    
    '''(just for debugging) print time for execution of every '500th' iteration (although it is hardcoded, doesn't have to be 500 could be anything)'''

    if s%500==0:
        t2=time.time()
        print (t2-t1), 'sec. : time for iteration '+str(s)


T2=time.time()
print '\n','......................................................','\n\n',(T2-T1), 'sec. : time for execution of master loop'
