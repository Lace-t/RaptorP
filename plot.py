from rapplot import Raptor
import matplotlib.pyplot as plt

angle=90
jet=Raptor(6.5e9,16.8e3,f'img_data_40_{angle}.00.h5',offset=30,unit='rg',
           root='output')
#print(jet.mas)
print(jet.keys())

if angle==1:
    jet.plot_stokes(230,figsize=(12,10))
    plt.suptitle(r'$1^\circ$ PowerLaw a=4',y=0.99,fontsize=15)
    plt.tight_layout()
    plt.savefig('1.png')
elif angle==10:
    jet.plot_evpa(86,size=((-3,3),(-8,8)),figsize=(10,16))
    plt.show()
elif angle==9:
    jet.plot_stokes(230,size=((-3,3),(-7,7)),figsize=(10,16))
    plt.suptitle(r'$9^\circ$ Thermal',y=0.99,fontsize=15)
    plt.tight_layout()
    plt.savefig('figure/9.png')
elif angle==17:
    jet.plot_stokes(230,size=((-200,200),(-300,200)),figsize=(16,16))
    plt.suptitle(r'$17^\circ$ PowerLaw a=4',y=0.99,fontsize=15)
    plt.tight_layout()
    plt.savefig('17.png',dpi=150)
elif angle=='17d':
    jet.plot_stokes(230,size=((-4,4),(-22,22)),figsize=(6,14),scale='log',vlim=([1e16,1e19],[0,0],[0,0],[0,0]))
    plt.suptitle(r'$17^\circ$ Thermal',y=0.99,fontsize=15)
    plt.tight_layout()
    plt.savefig('figure/17d.png',dpi=150)
elif angle==45:
    jet.plot_stokes(1)
    #plt.suptitle(r'$45^\circ$',y=0.99,fontsize=15)
    plt.tight_layout()
    #plt.savefig('45.png',dpi=150)
    plt.show()
if angle==90:
    jet.plot_stokes(1.5)
    #jet.plot_poldeg(1)
    #plt.suptitle(r'$90^\circ$ Thermal',y=0.99,fontsize=15)
    #plt.tight_layout()
    #plt.savefig('90.png',dpi=150)
    plt.show()
elif angle==135:
    jet.plot_stokes(230,size=((-3,3),(-22,23)),figsize=(6,16))
    plt.suptitle(r'$135^\circ$ Thermal',y=0.99,fontsize=15)
    plt.tight_layout()
    plt.savefig('figure/135.png',dpi=150)
elif angle==163:
    jet.plot_stokes(230,size=((-150,150),(-200,100)),figsize=(16,16))
    plt.suptitle(r'$163^\circ$ PowerLaw a=4',y=0.99,fontsize=15)
    plt.tight_layout()
    plt.savefig('163.png',dpi=150)
'''
if angle==1:
    jet.plot_evpa(230,figsize=(6,5),npix=10,arrow=4,width=0.005)
    plt.suptitle(r'$1^\circ$ Thermal',y=0.99,fontsize=15)
    plt.tight_layout()
    plt.savefig('figure/1p_evpa.png',dpi=150)
elif angle==17:
    jet.plot_evpa(230,size=((-3,3),(-10,12)),figsize=(4,8),
                  npix=10,arrow=1.5,width=0.01)
    plt.suptitle(r'$17^\circ$ tracer',y=0.99,fontsize=15)
    plt.tight_layout()
    plt.savefig('figure/17t_evpa.png',dpi=150)
elif angle==45:
    jet.plot_evpa(230,size=((-3,3),(-22,23)),figsize=(3,8),
                  npix=10,arrow=0.4,width=0.02)
    plt.suptitle(r'$45^\circ$ Thermal',y=0.99,fontsize=15)
    plt.tight_layout()
    plt.savefig('figure/45p_evpa.png',dpi=150)
'''
