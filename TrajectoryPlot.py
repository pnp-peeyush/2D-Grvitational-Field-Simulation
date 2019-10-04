import matplotlib.pyplot as plt
import matplotlib.animation
from GravitationalFieldandKinematicsSimulator import Simul_pos,IFD,dt,m,n

size_scale      = 5                             #multiplication factor for size of each representative 'dot' of an object
catchup_space   = 1                             #to contain the simulation inside the final figure
skip            = 2                            #similar concept as in GravitationalFieldandKinematicsSimulator file


fig=plt.figure(figsize=(13,9.8))
fig.subplots_adjust(left=0.03, right=0.9975, top=0.975, bottom=0.025)
ax1=fig.add_subplot(1,1,1)
ax1.axis('scaled')
plt.grid()

size=[]
col=[]
lines=[]
dots=[]

for i in range(0,len(IFD)):
    size.append(size_scale*IFD[i][3])
    col.append(IFD[i][-1])


for i in range(0,len(Simul_pos)):
    lines.append(ax1.plot(Simul_pos[i][0][:0],Simul_pos[i][1][:0],color=col[i]))
    dots.append(ax1.plot(Simul_pos[i][0][0],Simul_pos[i][1][0],marker='.',markersize=size[i],color=col[i],alpha=0.5))
    


a=[]
b=[]
d=[]
e=[]


ax1.set_xlim(-n,n)
ax1.set_ylim(-m,m)

def update(i):
    ax1.set_title('t+{0} sec'.format(i*dt))
    for c in range(0,len(Simul_pos)):
        line,=lines[c]
        line.set_data(Simul_pos[c][0][:i],Simul_pos[c][1][:i])
        dot,=dots[c]
        dot.set_data(Simul_pos[c][0][i],Simul_pos[c][1][i])
          
        if i >0:
            a.append(min(Simul_pos[c][0][:i]))
            b.append(max(Simul_pos[c][0][:i]))
            d.append(min(Simul_pos[c][1][:i]))
            e.append(max(Simul_pos[c][1][:i]))
    
    if i > 0:
        pass
        ax1.set_xlim(min(a)-catchup_space,max(b)+catchup_space)
        ax1.set_ylim(min(d)-catchup_space,max(e)+catchup_space)
    #plt.savefig('E:/BouncingBall/GravSim/render/foo'+str(i)+'.png', transparent=True)
        
        
ani = matplotlib.animation.FuncAnimation(fig, update, frames=range(0,len(Simul_pos[0][0]),skip),interval=10, repeat=False)   
mng = plt.get_current_fig_manager()
mng.window.showMaximized()
plt.show()