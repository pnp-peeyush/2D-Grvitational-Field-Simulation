import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation
from GravitationalFieldandKinematicsSimulator import Simul_field,X,Y,dt,n,m,skip

fig = plt.figure(figsize=(13,9.8))
ax = fig.add_subplot(111) 

im = ax.imshow(Simul_field[0][2],extent=[X.min(), X.max(), Y.min(), Y.max()],cmap='viridis_r',origin='lower',interpolation='gaussian')
cb = fig.colorbar(im,cax=fig.add_axes([0.925, 0.025, 0.025, 0.95]))
tx = ax.set_title('t+ 0 sec')
fig.subplots_adjust(left=0.03, right=0.9975, top=0.975, bottom=0.025)   

def update(i):
    vmax     = np.max(Simul_field[i][2])
    vmin     = np.min(Simul_field[i][2])
    im.set_data(Simul_field[i][2])
    im.set_clim(vmin, vmax)
    tx.set_text('t+{0} sec'.format(i*dt*skip))
    plt.savefig('E:/BouncingBall/GravSim/render/fld'+str(i)+'.png')
    return [im]
    
        
ani = matplotlib.animation.FuncAnimation(fig, update, frames=len(Simul_field),interval=10,repeat=False,blit=False)   
#ani.save('test_anim.mp4', fps=25, extra_args=['-vcodec', 'libx264'])
mng = plt.get_current_fig_manager()
mng.window.showMaximized()
plt.show()
