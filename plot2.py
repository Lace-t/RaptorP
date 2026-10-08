import numpy as np
from rapplot import Raptor
import matplotlib.pyplot as plt
import matplotlib.colors as colo
import matplotlib.ticker as ticker
from matplotlib import rcParams, rc
rcParams['font.family'] = 'Times New Roman'
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Times New Roman'
rcParams['mathtext.bf'] = 'Times New Roman:bold'
rcParams['mathtext.it'] = 'Times New Roman:italic'
rcParams['mathtext.default'] = 'it'
rcParams['mathtext.cal'] = 'cursive'
rc('text', usetex=True)
plt.rcParams.update({'font.size': 15})

c = 2.99792458e10          # cm/s
me = 9.10938356e-28        # g
e = 4.80320425e-10         # statcoulomb
sigmaT = 6.6524587158e-25 

##print("%e"%((6*np.pi*me*c)/sigmaT))
##print("%e"%((27*np.pi*me*e*c)/(sigmaT)**2))

angle=1.4
mode='show'
for i in [40]:
    jet=Raptor(5.2e5,9.4,f'img_data_{i}_{angle}0.h5',offset=0,unit='rg',
            root='output')
    
##    jet.plot(('tau',1.5),figsize=(5,8),scale='log',cmap='Greys')
##    plt.tight_layout()
##    if mode=='save':
##        plt.savefig(f'figure/tau_{i}.png')
##    elif mode=='show':
##        plt.show()
##    plt.close()

    # jet.plot(('F',1.5),figsize=(5,8),scale='log',cmap='Greys')
    # plt.tight_layout()
    # plt.show()

    
    jet.plot_evpa(4.6e5,figsize=(7,6))
    plt.tight_layout()
    if mode=='save':
        plt.savefig(f'figure/flux_{i}.png')
    elif mode=='show':
        plt.show()
    plt.close()

##    jet.plot_poldeg(1.5,figsize=(5,8))
##    plt.tight_layout()
##    if mode=='save':
##        plt.savefig(f'figure/fp_{i}.png')
##    elif mode=='show':
##        plt.show()
##    plt.close()

##    d=np.loadtxt(f'output/spectrum_{i}_{angle}.00.dat')
##    plt.figure(dpi=150)
##    plt.plot(d[:,0]/1e9,1e3*d[:,1])
##    plt.plot(d[:,0]/1e9,1e3*d[:,1],'r.')
##    plt.xscale('log')
##    plt.yscale('log')
##    plt.xlabel('Freq [GHz]')
##    plt.ylabel('FLux [mJy]')
##    plt.grid()
##    plt.tight_layout()
##    if mode=='save':
##        plt.savefig(f'figure/sed_{i}.png')
##    elif mode=='show':
##        plt.show()
##    plt.close()
##
##    plt.figure(dpi=150)
##    plt.plot(d[:,0]/1e9,np.sqrt(np.square(d[:,2])+np.square(d[:,3]))/d[:,1])
##    plt.plot(d[:,0]/1e9,np.sqrt(np.square(d[:,2])+np.square(d[:,3]))/d[:,1],'r.')
##    plt.xscale('log')
##    plt.xlabel('Freq [GHz]')
##    plt.ylabel('Degree of Polarization')
##    plt.grid()
##    plt.tight_layout()
##    if mode=='save':
##        plt.savefig(f'figure/sfp_{i}.png')
##    elif mode=='show':
##        plt.show()
##    plt.close()
