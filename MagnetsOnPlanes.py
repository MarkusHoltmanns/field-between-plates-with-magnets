import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import colors

#https://matplotlib.org/stable/api/asgen/matplotlib.pyplot.quiver.html#matplotlib.pyplot.quiver

fig = plt.figure(dpi=200)
ax = plt.figure(num=1, figsize=(12, 8)).add_subplot(projection='3d')
#ax.set_title('Magnetic Field Between Sets Of Magnets On Two Plates')

# Making the grid
x, y, z = np.meshgrid(np.arange(-0.4, 0.41, 0.03),
                      np.arange(-0.4, 0.41, 0.03),
                      np.arange(-0.4, 0.41, 0.1))


# magnetic moment
mx=0.422
my=0.0 #A m^2
mz=0.0
# location of the centre
sx = 0.
sy = 0.
sz = 0.




#https://stackoverflow.com/questions/49277753/plotting-cuboids
def cuboid_data2(o, size=(1,1,1)):
    X = [[[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0]],
         [[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0]],
         [[1, 0, 1], [1, 0, 0], [1, 1, 0], [1, 1, 1]],
         [[0, 0, 1], [0, 0, 0], [0, 1, 0], [0, 1, 1]],
         [[0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0]],
         [[0, 1, 1], [0, 0, 1], [1, 0, 1], [1, 1, 1]]]
    X = np.array(X).astype(float)
    for i in range(3):
        X[:,:,i] *= size[i]
    X += np.array(o)
    return X

def plotCubeAt2(positions,sizes=None,colors=None, **kwargs):
    if not isinstance(colors,(list,np.ndarray)): colors=["C0"]*len(positions)
    if not isinstance(sizes,(list,np.ndarray)): sizes=[(1,1,1)]*len(positions)
    g = []
    for p,s,c in zip(positions,sizes,colors):
        g.append( cuboid_data2(p, size=s) )
    return Poly3DCollection(np.concatenate(g), facecolors=np.repeat(colors,6), **kwargs)
    

positions = []
sizes = []
colors = []

start=0
length=13 # is equal to the amount of magnets in a line on the plate
ssx=0.007
ssy=0.015
ssz=0.015
dipole_positions=[]
distancetoplanes = 0.3
for i in range (start,start+length):
    for j in range (start,start+length):        
        for k in range (0,2):
            #print(i,j,k)
            planes = np.where(k<1,-1,1)
            radius = 0.07
            dx0 = sx + planes * distancetoplanes
            #dy0 = sy + radius * (i-(length-1)/2)
            #dz0 = sz + radius * (j-(length-1)/2)
            ii=(i-(length-1)/2)
            jj=(j-(length-1)/2)
            dy0 = sy + radius * ii *np.exp(-pow(1.2*ii/length,2))
            dz0 = sz + radius * jj *np.exp(-pow(1.2*jj/length,2))
            
            dx11 = dx0 -ssx 
            dx12 = dx0
            dy1 = dy0 -ssy/2 
            dz1 = dz0 -ssz/2 

            dipole_positions.append((dx0,dy0,dz0))

            positions.append((dx11,dy1,dz1))
            positions.append((dx12,dy1,dz1))

            sizes.append((ssx,ssy,ssz))
            sizes.append((ssx,ssy,ssz))
            colors.append("crimson")
            colors.append("limegreen")
            #print(2*np.pi*i/length)
    
dipole_positions = list(dipole_positions)
positions = list(positions) 
sizes = list(sizes)
pc = plotCubeAt2(positions,sizes,colors=colors, edgecolor="k", linewidth=0.1)
ax.add_collection3d(pc)    

# Magnetic field formula
# b_field = 1/(4*Pi)*(3*r*(r*m)/r^5-m/r^3)

# useful factors and variables
mu = 1.25663706212E-6 # H/m
factor = mu/(4*np.pi)
tempfieldx = 0.
tempfieldy = 0.
tempfieldz = 0.

#calculating the B-field at point (x,y,z)
for i in range (start,start+length*length*2):
    xx = dipole_positions[i][0]
    yy = dipole_positions[i][1]
    zz = dipole_positions[i][2]
    r=np.sqrt((x-xx)**2+(y-yy)**2+(z-zz)**2)
    rx=x-xx
    ry=y-yy
    rz=z-zz
    rm=rx*mx+ry*my+rz*mz
    
    tempfieldx = tempfieldx + factor*(3*rx*rm/r**5-mx/r**3)
    tempfieldy = tempfieldy + factor*(3*ry*rm/r**5-my/r**3)
    tempfieldz = tempfieldz + factor*(3*rz*rm/r**5-mz/r**3)
    
# calculating the field vectors in the 3d grid
multiplyfactor=1E2 #for the length of the arrows
u = tempfieldx*multiplyfactor
v = tempfieldy*multiplyfactor
w = tempfieldz*multiplyfactor
laenge=np.sqrt(u*u+v*v+w*w) #length of the printed b-field (arrow)

maxlength=0.03 #maximum length of the arrow in the 3d plot
uu = maxlength*u/laenge
vv = maxlength*v/laenge
ww = maxlength*w/laenge
#check if any arrow is longer than maxlength and resize it if necessary
u=np.where(laenge<maxlength,u,uu)
v=np.where(laenge<maxlength,v,vv)
w=np.where(laenge<maxlength,w,ww)
laenge=np.sqrt(u*u+v*v+w*w)

# Color by azimuthal angle
c = laenge #np.arctan2(v, u)
# Flatten and normalize
c = (c.ravel() - c.min()) / c.ptp()
# Repeat for each body line and two head lines
c = np.concatenate((c, np.repeat(c, 2)))
# Colormap
c = plt.cm.hsv(c)

ax.quiver(x, y, z, u, v, w, linewidths=0.2, color=c)#, length=0.1, normalize=True)

#set the limits of the shown axis in the 3d plot
ax.set_xlim([-0.4,0.41])
ax.set_ylim([-0.4,0.41])
ax.set_zlim([-0.4,0.41])

#labels of the axis
xLabel = ax.set_xlabel('x [m]', fontsize=8, linespacing=3.1)
yLabel = ax.set_ylabel('y [m]', fontsize=8, linespacing=3.1)
zLabel = ax.set_zlabel('z [m]', fontsize=8, linespacing=3.1)

ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)
ax.zaxis.set_tick_params(labelsize=8)

#plt.title('Magnetic Field Between Sets Of Magnets On Two Plates')

#fig.suptitle("Figure title", fontsize=16)
# Figure subtitle
#fig.text(0.5, 0.9, "Figure subtitle", horizontalalignment="center")
# Figure footer title
fig.text(0.5, 0.015, "Magnetic Field Between Sets Of Magnets On Two Plates",fontsize=8, horizontalalignment="center")

# If you wanna save your image file in a specific directory
#filename = "/Users/Holtmanns/Desktop/Plot.png"
#plt.savefig(filename, dpi=1000)

plt.show()
