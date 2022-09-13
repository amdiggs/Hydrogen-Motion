import matplotlib.pyplot as plt
import numpy as np

def sphere(at,ax):
    rads=[0.5, 1.2]
    clrs=["m","r","b"]
    r=rads[at.type -1]
    theta=np.arange(0,3.24,0.1)
    phi=np.arange(0,6.38,0.1)
    T,P=np.meshgrid(theta,phi)
    X=r*np.sin(T)*np.cos(P) + at.x
    Y=r*np.sin(T)*np.sin(P) + at.y
    Z=r*np.cos(T) + at.z
    ax.plot_surface(X,Y,Z, color=clrs[at.type])
    return

def bond(at1, at2, ax):
    clrs=["m","r","b"]
    mid=(at1.coords + at2.coords)/2
    x, y, z = [[at1.x,mid[0]], [at1.y, mid[1]], [at1.z, mid[2]]]
    x2, y2, z2 = [[at2.x,mid[0]], [at2.y, mid[1]], [at2.z, mid[2]]]
    ax.plot(x,y,z, linewidth=10, color=clrs[at1.type])
    ax.plot(x2,y2,z2, linewidth=10, color=clrs[at2.type])
    return

def draw_config(at):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    sphere(at,ax)
    for n in at.neighbors:
        sphere(n,ax)
        bond(at,n,ax)
    plt.show()