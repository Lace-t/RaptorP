import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
import h5py

#font = {'size': 16}
#matplotlib.rc('font', **font)

G = 6.674e-8
MSUN = 1.989e33
SPEED_OF_LIGHT = 2.998e10
KPC = 3.086e21
SEC2DAY = 1/86400.
DEG2MAS= 206264.806*1000

class Raptor(object):
    '''>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    Raptor class is a reconsitution of rapplot.py
    Initialization
    --------------
    M              : double
                     Mass of black hole(unit: Msun)
    d              : double
                     distance(unit:kpc)
    name           : str
                     number of datafile
    root(optional) : str
                     location of datafile
    --------------        
    Written by Xufan Hu 2024.11.30
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    '''
    def __init__(self,M,d,name,root=''):
        #reading .h5 file
        self.__name=root+'/'+name
        print("Reading keys from: ",self.__name)
        
        self.__data=h5py.File(self.__name,'r')
        self.__keys=[key for key in self.__data.keys()]
        #consider the possibility of multi-frequencies
        self.__freq=np.array([])
        for i in self.__keys:
            if i[0]=='I':
                self.__freq=np.append(self.__freq,eval(i[1:])/1e9)
            else:
                break
        print("frequencies: ",self.__freq)
        #calculate the scale
        self.__M=M*MSUN
        self.__d=d*KPC
        self.rg=G*self.__M/SPEED_OF_LIGHT**2#show it by Raptor.rg
        self.mas=self.rg/self.__d*DEG2MAS#show it by Raptor.mas
        self.__tunit =self.rg/SPEED_OF_LIGHT

    def keys(self):
        '''Return the index of the datafile

        Returns
        -------
        keys : list
        '''
        return self.__keys

    def data(self):
        '''Return the data in the datafile

        Returns
        -------
        data : h5py._hl.files.File
        '''
        return self.__data

    def freq(self):
        '''Return the frequencies of data

        Returns
        -------
        freq : np.array
               unit: GHz
        '''
        return self.__freq

    def order(self,ind):
        '''Return the index of a certain Stokes parameter at a certain frequency

        parameters
        ----------
        ind : any combination data type
              ('stokes parameter', frequency) e.g. ('I',230)

        Returns
        -------
        order : int
        '''   
        dic={'I':0,'Q':1,'U':2,'V':3,'tau':4,'F':5}
        if dic[ind[0]]<4:
            return dic[ind[0]]*len(self.__freq)+np.where(self.__freq==ind[1])[0][0]
        else:
            return dic[ind[0]]*len(self.__freq)+np.where(self.__freq==ind[1])[0][0]+2

    def close(self):
        '''close the h5 file
        '''
        self.__data.close()
        gc.collect()
        self.__data = None
    
    def plot(self,ind,label='',size=None,vmin=0,vmax=0,figsize=(6,8),scale='linear',cmap="afmhot",ax=None,fig=None):
        '''Plot a parameter of any index

        parameters
        ----------
        ind               : int or any combination data type
                            int                       : index
                            any combination data type : ('stokes parameter', frequency) e.g. ('I',230)
        label             : str
                            label of colorbar
        size(optional)    : None or tuple with two tuples
                            x and y span of figure(unit: rg)
                            None  : use the limit of canvas
                            tuple : set it manually ((xmin,xmax),(ymin,ymax))
                            (default: None)
        vmin(optional)    : double
                            minimum value in the colorbar
                            (default: 0)
        vmax(optional)    : double
                            maximum value in the colorbar
                            (default: 0)
        figsize(optional) : tuple with two elements
                            figure size
                            (default: (6,8))
        scale(optional)   : "log" or "linear"
                            scale of parameter
                            (default: 'linear')
        cmap(optional)    : str
                            cmap in the colormap
                            (default: 'afmhot')
        ax(optional)      : matplotlib.ax
                            subplot, if not defined, the function will generate one automatically
                            (default: None)
        fig(optional)     : matplotlib.fig
                            canvas, if not defined, the function will generate one automatically
                            (default: None)                          
        '''
        #determine ind
        if type(ind)==int:
            index=ind
        else:
            index=self.order(ind)
            
        #If not specify, find extreme value
        if vmin==vmax:
            vmin=np.min(self.__data[self.__keys[index]][0])
            vmax=np.max(self.__data[self.__keys[index]][0])
            for i in range(1,len(self.__data[self.__keys[index]])):
                vmin=min(vmin,np.min(self.__data[self.__keys[index]][i]))
                vmax=max(vmax,np.max(self.__data[self.__keys[index]][i]))
            print(label,"vmin=%.3e"%vmin,"vmax=%.3e"%vmax)
        #single panel or multiple panels
        if ax==None and fig==None:
            fig,ax= plt.subplots(figsize=figsize)
        #plot scale
        if scale == 'log':
            for i in range(0,len(self.__data[self.__keys[index]])):
                self.__pixels=int(np.sqrt(len(self.__data[self.__keys[index]][i])))
                array=((np.reshape(self.__data[self.__keys[index]][i],(self.__pixels,self.__pixels))))
                alpha=((np.reshape(self.__data['alpha'][i],(self.__pixels,self.__pixels))))*self.mas
                beta=((np.reshape(self.__data['beta'][i],(self.__pixels,self.__pixels))))*self.mas
                
                #ax.set_aspect('equal')                
                figure=ax.pcolormesh(alpha,beta,np.log10(array+1),vmin=np.log10(vmin+1),
                                     vmax=np.log10(vmax),cmap=cmap,shading='auto')
        elif scale=='linear':
            for i in range(0,len(self.__data[self.__keys[index]])):  
                pixels=int(np.sqrt(len(self.__data[self.__keys[index]][i])))
                array=((np.reshape(self.__data[self.__keys[index]][i],(pixels,pixels))))
                alpha=((np.reshape(self.__data['alpha'][i],(pixels,pixels))))*self.mas
                beta=((np.reshape(self.__data['beta'][i],(pixels,pixels))))*self.mas
                
                #ax.set_aspect('equal')
                figure=ax.pcolormesh(alpha,beta,array,vmin=vmin,vmax=vmax,cmap=cmap,shading='auto')
        else :
            raise ValueError("only 'log' or 'linear' is valid for scale.")
        
        fig.colorbar(figure,label=label,ax=ax)
        if size!=None:
            ax.set_xlim(size[0][0]*self.mas,size[0][1]*self.mas)
            ax.set_ylim(size[1][0]*self.mas,size[1][1]*self.mas)
        #define scientific notation
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True) 
        formatter.set_powerlimits((0,0)) 
        ax.yaxis.set_major_formatter(formatter)
        ax.xaxis.set_major_formatter(formatter)

        ax.tick_params(axis='both',which='both',direction='out')
        ax.set_xlabel(r"x [mas]",loc='left')
        ax.set_ylabel(r"y [mas]")
        plt.tight_layout()

    def plot_pol(self,freq,vmin=0,vmax=0,size=None,scale='linear',figsize=(6,8),cmap="BuPu"):
        '''Plot the linear polarization at a certain frequency

        parameters
        ----------
        freq              : double
                            frequency
        size(optional)    : None or tuple with two tuples
                            x and y span of figure(unit: rg)
                            None  : use the limit of canvas
                            tuple : set it manually ((xmin,xmax),(ymin,ymax))
                            (default: None)
        vmin(optional)    : double
                            minimum value in the colorbar
                            (default: 0)
        vmax(optional)    : double
                            maximum value in the colorbar
                            (default: 0)
        scale(optional)   : "log" or "linear"
                            scale of parameter
                            (default: 'linear')
        figsize(optional) : tuple with two elements
                            figure size
                            (default: (6,8))
        cmap(optional)    : str
                            cmap in the colormap
                            (default: 'BuPu')
        '''
        index=np.where(self.__freq==freq)[0][0]
        #Find vmin & vmax
        if vmin==vmax:
            block=np.sqrt(np.square(self.__data[self.__keys[len(self.__freq)+index]][0])+
                          np.square(self.__data[self.__keys[2*len(self.__freq)+index]][0]))
            vmin=np.min(block)
            vmax=np.max(block)
            for i in range(1,len(self.__data[self.__keys[index]])):
                block=np.sqrt(np.square(self.__data[self.__keys[len(self.__freq)+index]][i])+
                          np.square(self.__data[self.__keys[2*len(self.__freq)+index]][i]))
                vmax=max(vmax,np.max(block))
                vmin=max(vmin,np.min(block))
            print("Linear Polarization : vmin=%.3e"%vmin,"vmax=%.3e"%vmax)
        
        fig,ax= plt.subplots(figsize=figsize)
        
        if scale=='log':
            for i in range(0,len(self.__data[self.__keys[index]])):  
                pixels=int(np.sqrt(len(self.__data[self.__keys[index]][i])))
                array_Q=((np.reshape(self.__data[self.__keys[len(self.__freq)+index]][i],(pixels,pixels))))
                array_U=((np.reshape(self.__data[self.__keys[2*len(self.__freq)+index]][i],(pixels,pixels))))
                array = np.sqrt(np.power(array_Q,2)+np.power(array_U,2))
                
                alpha=((np.reshape(self.__data['alpha'][i],(pixels,pixels))))*self.mas
                beta=((np.reshape(self.__data['beta'][i],(pixels,pixels))))*self.mas
                
                #ax.set_aspect('equal')
                figure=ax.pcolormesh(alpha,beta,np.log10(array+1e-40),vmin=np.log10(vmin+1e-40),
                                     vmax=np.log10(vmax),cmap=cmap,shading='auto')
        if scale=='linear':
            for i in range(0,len(self.__data[self.__keys[index]])):  
                pixels=int(np.sqrt(len(self.__data[self.__keys[index]][i])))
                array_Q=((np.reshape(self.__data[self.__keys[len(self.__freq)+index]][i],(pixels,pixels))))
                array_U=((np.reshape(self.__data[self.__keys[2*len(self.__freq)+index]][i],(pixels,pixels))))
                array = np.sqrt(np.power(array_Q,2)+np.power(array_U,2))
                
                alpha=((np.reshape(self.__data['alpha'][i],(pixels,pixels))))*self.mas
                beta=((np.reshape(self.__data['beta'][i],(pixels,pixels))))*self.mas
                
                #ax.set_aspect('equal')
                figure=ax.pcolormesh(alpha,beta,array,vmin=0,vmax=vmax,cmap=cmap,shading='auto')

        fig.colorbar(figure,label='Linear Pol (%.2eGHz)'%freq,ax=ax)
        if size!=None:
            ax.set_xlim(size[0][0]*self.mas,size[0][1]*self.mas)
            ax.set_ylim(size[1][0]*self.mas,size[1][1]*self.mas)
        #define scientific notation
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True) 
        formatter.set_powerlimits((0,0)) 
        ax.yaxis.set_major_formatter(formatter)
        ax.xaxis.set_major_formatter(formatter)

        ax.tick_params(axis='both',which='both',direction='out')
        ax.set_xlabel(r"x",loc='left')
        ax.set_ylabel(r"y")
        plt.tight_layout()

    def plot_poldeg(self,freq,size=None,label='|m|',figsize=(6,8),cmap="afmhot"):
        '''Plot the fraction of polarization at a certain frequency

        parameters
        ----------
        freq              : double
                            frequency
        label             : str
                            label of colorbar
                            (default: '|m|')
        size(optional)    : None or tuple with two tuples
                            x and y span of figure(unit: rg)
                            None  : use the limit of canvas
                            tuple : set it manually ((xmin,xmax),(ymin,ymax))
                            (default: None)
        vmin(optional)    : double
                            minimum value in the colorbar
                            (default: 0)
        vmax(optional)    : double
                            maximum value in the colorbar
                            (default: 0)
        figsize(optional) : tuple with two elements
                            figure size
                            (default: (6,8))
        cmap(optional)    : str
                            cmap in the colormap
                            (default: 'afmhot')
        '''
        index=np.where(self.__freq==freq)[0][0]
        #Find vmax
        vmax=np.max(self.__data[self.__keys[index]][0])
        for i in range(1,len(self.__data[self.__keys[index]])):
            vmax=max(vmax,np.max(self.__data[self.__keys[index]][i]))
        print("Stokes I : vmax=%.3e"%vmax)
        
        fig,ax= plt.subplots(figsize=figsize)
        for i in range(0,len(self.__data[self.__keys[index]])):  
            pixels=int(np.sqrt(len(self.__data[self.__keys[index]][i])))
            array_I=((np.reshape(self.__data[self.__keys[index]][i],(pixels,pixels))))
            array_Q=((np.reshape(self.__data[self.__keys[len(self.__freq)+index]][i],(pixels,pixels))))
            array_U=((np.reshape(self.__data[self.__keys[2*len(self.__freq)+index]][i],(pixels,pixels))))
            array = np.sqrt(np.power(array_Q,2)+np.power(array_U,2))/(array_I+1e-40)
            array[array_I/vmax<1e-7]=0
            
            alpha=((np.reshape(self.__data['alpha'][i],(pixels,pixels))))*self.mas
            beta=((np.reshape(self.__data['beta'][i],(pixels,pixels))))*self.mas
            
            #ax.set_aspect('equal')
            figure=ax.pcolormesh(alpha,beta,array,vmin=0,vmax=1,cmap=cmap,shading='auto')

        fig.colorbar(figure,label=label+'(%.2eGHz)'%freq,ax=ax)
        if size!=None:
            ax.set_xlim(size[0][0]*self.mas,size[0][1]*self.mas)
            ax.set_ylim(size[1][0]*self.mas,size[1][1]*self.mas)
        #define scientific notation
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True) 
        formatter.set_powerlimits((0,0)) 
        ax.yaxis.set_major_formatter(formatter)
        ax.xaxis.set_major_formatter(formatter)

        ax.tick_params(axis='both',which='both',direction='out')
        ax.set_xlabel(r"x [mas]",loc='left')
        ax.set_ylabel(r"y [mas]")
        plt.tight_layout()

    def plot_evpa(self,freq,size=None,mask=0.1,floor=0.1,num=1,arrow=0.5,width=0.02,vmin=0,vmax=0,figsize=(6,8),scale='linear',cmap="afmhot",ax=None,fig=None):
        '''Plot the EVPA at a certain frequency

        parameters
        ----------
        freq              : double
                            frequency
        size(optional)    : None or tuple with two tuples
                            x and y span of figure(unit: rg)
                            None  : use the limit of canvas
                            tuple : set it manually ((xmin,xmax),(ymin,ymax))
                            (default: None)
        mask(optional)    : double
                            threshold of I/I_max for plotting the EVPAs
                            (default: 0.1)
        floor(optional)   : double
                            threshold of degree of polarization for plotting the EVPAs
                            (default: 0.1)
        num(optional)     : int
                            the number of EVPA arrow(s) in one block
                            (default: 1)
        arrow(optional)   : double
                            scale of arrow, the less this value is, the longer the length of the arrow
                            (default: 0.5)
        width(optional)   : double
                            width of arrow
                            (default: 0.02)
        vmin(optional)    : double
                            minimum value in the colorbar
                            (default: 0)
        vmax(optional)    : double
                            maximum value in the colorbar
                            (default: 0)
        figsize(optional) : tuple with two elements
                            figure size
                            (default: (6,8))
        scale(optional)   : "log" or "linear"
                            scale of parameter
                            (default: 'linear')
        cmap(optional)    : str
                            cmap in the colormap
                            (default: 'afmhot')
        ax(optional)      : matplotlib.ax
                            subplot, if not defined, the function will generate one automatically
                            (default: None)
        fig(optional)     : matplotlib.fig
                            canvas, if not defined, the function will generate one automatically
                            (default: None)
        '''
        index=np.where(self.__freq==freq)[0][0]
        #If not specify, find extreme value
        if vmin==vmax:
            vmin=np.min(self.__data[self.__keys[index]][0])
            vmax=np.max(self.__data[self.__keys[index]][0])
            for i in range(1,len(self.__data[self.__keys[index]])):
                vmin=min(vmin,np.min(self.__data[self.__keys[index]][i]))
                vmax=max(vmax,np.max(self.__data[self.__keys[index]][i]))
            print('Stokes I(%.2eGHz):'%freq," vmin=%.3e"%vmin," vmax=%.3e"%vmax)
        
        norm=np.max(np.sqrt(np.power(self.__data[self.__keys[index+len(self.__freq)]][0],2)+
                            np.power(self.__data[self.__keys[index+2*len(self.__freq)]][0],2)))                                                             
        for i in range(1,len(self.__data[self.__keys[1]])):
            norm=max(norm,np.max(np.sqrt(np.power(self.__data[self.__keys[index+len(self.__freq)]][i],2)+
                                         np.power(self.__data[self.__keys[index+2*len(self.__freq)]][i],2))))
        #print("sqrt(Q^2+U^2)=%.3e"%norm)

        dpix=np.ceil(20/num).astype(int)
        #single panel or multiple panels
        if ax==None and fig==None:
            fig,ax= plt.subplots(figsize=figsize)
        if scale == 'log':
            for i in range(0,len(self.__data[self.__keys[index]])):  
                pixels=int(np.sqrt(len(self.__data[self.__keys[index]][i])))
                I=((np.reshape(self.__data[self.__keys[index]][i],(pixels,pixels))))
                Q=((np.reshape(self.__data[self.__keys[index+len(self.__freq)]][i],(pixels,pixels))))
                U=((np.reshape(self.__data[self.__keys[index+2*len(self.__freq)]][i],(pixels,pixels))))
                alpha=((np.reshape(self.__data['alpha'][i],(pixels,pixels))))*self.mas
                beta=((np.reshape(self.__data['beta'][i],(pixels,pixels))))*self.mas
                #calculate EVPA
                evpa = (180/np.pi)*0.5*np.arctan2(U,Q+1e-40)
                evpa[evpa>90]-=180            
                vxp=np.sqrt(Q*Q+U*U)*np.sin(evpa*np.pi/180)/norm
                vyp=-np.sqrt(Q*Q+U*U)*np.cos(evpa*np.pi/180)/norm
    
                #ax.set_aspect('equal')
                figure=ax.pcolormesh(alpha,beta,np.log10(I+1e-40),vmin=np.log10(vmin+1),vmax=np.log10(vmax),cmap=cmap,shading='auto')
                #use a threshold to determine whether plot EVPA
                if np.mean(I)>mask*vmax and np.mean(np.sqrt(np.square(vxp)+np.square(vyp)))>floor:
                    ax.quiver(alpha[dpix//2::dpix,dpix//2::dpix], beta[dpix//2::dpix,dpix//2::dpix],
                              vxp[dpix//2::dpix,dpix//2::dpix],vyp[dpix//2::dpix,dpix//2::dpix],pivot='mid',angles='xy',
                              headwidth=0,headlength=0,headaxislength=0,scale=arrow,width=width,color='green',scale_units='xy')
            fig.colorbar(figure,label=r'$\log_{10}$ Stokes I(%.2eGHz)'%freq,ax=ax)
        elif scale=='linear':
            for i in range(0,len(self.__data[self.__keys[1]])):  
                pixels=int(np.sqrt(len(self.__data[self.__keys[1]][i])))
                I=((np.reshape(self.__data[self.__keys[index]][i],(pixels,pixels))))
                Q=((np.reshape(self.__data[self.__keys[index+len(self.__freq)]][i],(pixels,pixels))))
                U=((np.reshape(self.__data[self.__keys[index+2*len(self.__freq)]][i],(pixels,pixels))))
                alpha=((np.reshape(self.__data['alpha'][i],(pixels,pixels))))*self.mas
                beta=((np.reshape(self.__data['beta'][i],(pixels,pixels))))*self.mas
                #calculate EVPA
                evpa = (180/np.pi)*0.5*np.arctan2(U,Q+1e-40)
                evpa[evpa>90]-=180            
                vxp=np.sqrt(Q*Q+U*U)*np.sin(evpa*np.pi/180)/norm*self.mas
                vyp=-np.sqrt(Q*Q+U*U)*np.cos(evpa*np.pi/180)/norm*self.mas
    
                #ax.set_aspect('equal')
                figure=ax.pcolormesh(alpha,beta,I,vmin=vmin,vmax=vmax,cmap=cmap,shading='auto')
                #use a threshold to determine whether plot EVPA
                if np.mean(I)>mask*vmax and np.mean(np.sqrt(np.square(vxp)+np.square(vyp)))>floor*self.mas:
                    ax.quiver(alpha[dpix//2::dpix,dpix//2::dpix], beta[dpix//2::dpix,dpix//2::dpix],
                              vxp[dpix//2::dpix,dpix//2::dpix],vyp[dpix//2::dpix,dpix//2::dpix],pivot='mid',angles='xy',
                              headwidth=0,headlength=0,headaxislength=0,scale=arrow,width=width,color='green',scale_units='xy')
            fig.colorbar(figure,label='Stokes I(%.2eGHz)'%freq,ax=ax)
        else :
            raise ValueError("only 'log' or 'linear' is valid for scale.")
        
        if size!=None:
            ax.set_xlim(size[0][0]*self.mas,size[0][1]*self.mas)
            ax.set_ylim(size[1][0]*self.mas,size[1][1]*self.mas)
        #define scientific notation
        formatter = ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True) 
        formatter.set_powerlimits((0,0)) 
        ax.yaxis.set_major_formatter(formatter)
        ax.xaxis.set_major_formatter(formatter)

        ax.tick_params(axis='both',which='both',direction='out')
        ax.set_xlabel(r"x [mas]",loc='left')
        ax.set_ylabel(r"y [mas]")
        plt.tight_layout()
        
    def plot_stokes(self,freq,size=None,figsize=(12,16),evpa=None,scale='linear',cmaps=("afmhot","RdBu","RdBu","RdBu"),vlim=([0,0],[0,0],[0,0],[0,0])):
        '''Plot 2*2 Stokes panels at a certain frequency

        parameters
        ----------
        freq              : double
                            frequency
        size(optional)    : None or tuple with two tuples
                            x and y span of figure(unit: rg)
                            None  : use the limit of canvas
                            tuple : set it manually ((xmin,xmax),(ymin,ymax))
                            (default: None)
        figsize(optional) : tuple with two elements
                            figure size
                            (default: (6,8))
        evpa(optional)    : True or any combination data type
                            None                  : plot "I" parameter only
                            True                  : plot EVPA with default setup(mask=0.1,floor=0.1,num=1,arrow=0.5,width=0.02,see plot_evpa())
                            combination data type : plot EVPA with user-defined (mask,floor,num,arrow,width)
        scale(optional)   : "log" or "linear"
                            scale of parameter
                            (default: 'linear')
        cmaps(optional)   : tuple with four elements
                            cmap in every panel
                            (default: ("afmhot","RdBu","RdBu","RdBu"))
        vlim(optional)    : tuple with 4 lists
                            [vmin,vmax] for each panel
                            (default: ([0,0],[0,0],[0,0],[0,0]))
        '''
        dic=('I','Q','U','V')
        labels=('Stokes I(%.2eGHz)'%freq,'Stokes Q(%.2eGHz)'%freq,'Stokes U(%.2eGHz)'%freq,'Stokes V(%.2eGHz)'%freq)
        
        fig,ax= plt.subplots(2,2,figsize=figsize)
        axs=(ax[0][0],ax[0][1],ax[1][0],ax[1][1])
        #plot Stokes I with or without EVPA
        if evpa==None:
            self.plot((dic[0],freq),labels[0],size,vmin=vlim[0][0],vmax=vlim[0][1],figsize=(figsize[0]/2,figsize[1]/2),
                      scale=scale,cmap=cmaps[0],ax=axs[0],fig=fig)
        elif evpa==True:
            self.plot_evpa(freq,size,vmin=vlim[0][0],vmax=vlim[0][1],figsize=(figsize[0]/2,figsize[1]/2),
                           scale=scale,cmap=cmaps[0],ax=axs[0],fig=fig)
        else:
            self.plot_evpa(freq,size,mask=evpa[0],floor=evpa[1],num=evpa[2],arrow=evpa[3],width=evpa[4],
                           vmin=vlim[0][0],vmax=vlim[0][1],figsize=(figsize[0]/2,figsize[1]/2),
                           scale=scale,cmap=cmaps[0],ax=axs[0],fig=fig)        
        #plot Stokes Q,U,V
        for i in range(1,4):
            self.plot((dic[i],freq),labels[i],size,vmin=vlim[i][0],vmax=vlim[i][1],figsize=(figsize[0]/2,figsize[1]/2),
                      cmap=cmaps[i],ax=axs[i],fig=fig)
