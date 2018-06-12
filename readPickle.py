import cPickle as pickle
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.colors as col

def readPickle():
    [zsL,t1sL,citL,fname]=pickle.load(open('/home/grecu/Citation/nov23.pklz','rb'))
    return np.array(zsL),t1sL,citL,fname

def plotret(t1sL,h,zwC2D,vmax,title,fname):
    plt.figure()
    z2=np.ma.array(zwC2D,mask=zwC2D<1)
    plt.pcolormesh(t1sL,h,z2.T,vmin=1,vmax=vmax,cmap="jet")
    plt.xlabel('Time')
    plt.ylabel('Height')
    plt.title(title)
    plt.ylim(0.,8.)
    plt.savefig(fname)
    plt.colorbar()

def plotWC(qKu,qW):
    fig=plt.figure()
    ax=plt.subplot(111)
    ax.set_aspect('equal')
    ax.scatter(qKu,qW,s=2)
    plt.xlabel("Ku-Ka PWC (g/m^3)")
    plt.ylabel("W-band PWC (g/m^3)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.xlim(0.03,1)
    plt.ylim(0.03,1)
