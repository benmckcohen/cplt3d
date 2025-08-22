import cplt3d.generator_funcs
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as colormaps
import numpy as np
import os
import tqdm.notebook as tqdm
# import sklearn.neighbors as neighbors
import time
import matplotlib.style as mplstyle
# import multiprocess
import copy
import numpy as np
# N_cpu = multiprocess.cpu_count()
# print(N_cpu)
# mplstyle.use('fast')

def in_box(pts,x,y,z,dx,dy,dz):
    in_x = (pts[:,0] > x) & (pts[:,0] < x + dx)
    in_y = (pts[:,1] > y) & (pts[:,1] < y + dy)
    in_z = (pts[:,2] > z) & (pts[:,2] < z + dz)
    return in_x & in_y & in_z

class GadgetSimulation:
    def __init__(self,name,coords = None,verbose = False,particle_type = 'PartType1'):
        '''Input the gadget_result directory path'''
        l = os.listdir(name)
        hdfs = []
        file_names = []
        file_nums = []
        for f in l:
            if "snapshot" in f:
                file_name = name + '/' + f
                file_num = int(''.join([i for i in f if str.isnumeric(i)][:-1]))
                file_names.append(file_name)
                file_nums.append(file_num)
        file_names = [x for _,x in sorted(zip(file_nums, file_names))]
        for file_name in file_names:
            try:
                file = h5py.File(file_name, 'r')
                hdfs.append(file)
            except:
                print(file_name)
        self.hdfs = hdfs
        self.cs = coords
        self.identifications = None
        self.vels = None
        self.verbose = verbose
        self.particle_type = particle_type
        # print(self.hdfs)
    def __iter__(self):
        return iter(self.hdfs)
    def __getitem__(self,i):
        return self.hdfs[i]
    def __len__(self):
        return len(self.hdfs)
    def _fetch(self,name,index = None):
        vals = []
        for tab in self:
            val = np.array(tab[self.particle_type][name]).T
            if not index is None:
                val = val[:,index]
            vals.append(val)
        return np.array(vals)
    def find_particles(self,ids,i,name):
        tab = self[i]
        val = np.array(tab[self.particle_type][name])
        identifications = tab[self.particle_type]['ParticleIDs']
        mask = np.in1d(identifications,ids,assume_unique=True)
        return val[mask]
        
    def coords(self):
        if self.cs is None:
            self.cs = self._fetch("Coordinates")
        return self.cs
    def vs(self):
        if self.vels is None:
            self.vels = self._fetch("Velocities")
        return self.vels
    def ids(self):
        if self.identifications is None:
            self.identifications = self._fetch("ParticleIDs")
        return self.identifications
    def x(self):
        return self._fetch("Coordinates",index = 0)
    def y(self):
        return self._fetch("Coordinates",index = 1)
    def z(self):
        return self._fetch("Coordinates",index = 2)
    def vx(self):
        return self._fetch("Velocities",index = 0)
    def vy(self):
        return self._fetch("Velocities",index = 1)
    def vz(self):
        return self._fetch("Velocities",index = 2)
 
    def KernelDensity(self,i,shape,**kwargs):
        coords = self.coords()[i]
        x,y,z = coords
        min_x,max_x = np.min(x),np.max(x)
        min_y,max_y = np.min(y),np.max(y)
        min_z,max_z = np.min(z),np.max(z)
        if type(shape) != int:
            N_x,N_y,N_z = shape
        else:
            N_x,N_y,N_z = shape,shape,shape
        X,Y,Z = np.linspace(min_x,max_x,N_x),np.linspace(min_y,max_y,N_y),np.linspace(min_z,max_z,N_z)
        X,Y,Z = np.meshgrid(X,Y,Z)
        X_score = np.array([X.flatten(),Y.flatten(),Z.flatten()]).T
        kde = neighbors.KernelDensity(**kwargs)
        kde.fit(coords.T)
        Y = kde.score_samples(X_score)
        Y = np.exp(Y)
        return kde,X_score,Y
    def _particles_in(self,i,particle_func):
        coords = self.coords()[i]
        return particle_func(coords)
    def particles_in_square(self,bounds,i):
        def particle_func(coords):
            x = bounds[0][0]
            dx = bounds[0][1] - bounds[0][0]
            y = bounds[1][0]
            dy = bounds[1][1] - bounds[1][0]
            z = bounds[2][0]
            dz = bounds[2][1] - bounds[2][0]
            return in_box(coords.T,x,y,z,dx,dy,dz)
        return self._particles_in(i,particle_func)

    def plotting(f,cbar = True):
        def wrap(self,i,savevals=None,size = (6,8),axlab=['X','Y','Z'],ax = None,fig = None,**kwargs):
            if ax is None or fig is None:
                print(f'obtained {ax} and {fig}, need both to use input')
                if self.verbose:
                    print('Making new axis')
                fig = plt.figure(figsize = size)
                ax = fig.add_subplot(projection='3d')
                
            if not axlab is None:
                ax.set_xlabel(axlab[0])
                ax.set_ylabel(axlab[1])
                ax.set_zlabel(axlab[2])
            
            res = f(self,i,ax = ax,**kwargs)
            
            if savevals!=None:
                fig.savefig(savevals['name'],dpi = savevals['dpi'])
                plt.close()
            return fig,ax,res
        return wrap
    def animated(autoparallel = False):
        def wrap_wrap(f):
            def wrap(self,i,arr = None,anim_dir=None,ft='.png',forceparallel = False,**kwargs):
                if i=='a':
                    reses = []
                    kwargs_iter = {}
                    for key in copy.copy(list(kwargs.keys())):
                        if "iter_"+key in kwargs.keys() and kwargs["iter_"+key]:
                            kwargs_iter[key] = kwargs[key]
                            kwargs.pop("iter_"+key)
                    for j in tqdm.trange(len(self)):
                        for key in kwargs_iter.keys():
                            kwargs[key] = kwargs_iter[key][j]
                        res = f(self,j,savevals = {"name":anim_dir+'/'+str(j)+ft,"dpi":300},**kwargs)
                        reses.append(res)
                    return reses
                else:
                    return f(self,i,**kwargs)
            return wrap
        return wrap_wrap
    @animated(autoparallel=True)
    @plotting
    def scatter3d(self,i,coords = None,N=10,ms = 0.01,ax = None):
        coords = self.coords()[i]
        ax.plot(coords[0][::N],coords[1][::N],coords[2][::N],'ro',ms = ms)
    @animated(autoparallel=False)
    @plotting
    def kdescatter3d(self,i,shape=(100,100),bandwidth=0.1,ms=0.01,ax=None):
        if self.verbose:
            t = time.time()
            print("KDE...")
        kde,X_score,Y = self.KernelDensity(i,shape,bandwidth = bandwidth,kernel='tophat')
        if self.verbose:
            print("finished kde",time.time()-t,'seconds')
            print("plotting...")
        X,Y,Z = X_score.T[0],X_score.T[1],X_score.T[2]
        return ax.scatter(X,Y,Z,c = Y)
    @animated(autoparallel=False)
    @plotting
    def hist3d(self,i,ax=None,N=10,add_colorbar = True,cmap = colormaps.get_cmap('viridis'),shrink = 0.8,pad = -0.2,**kwargs):
        print('coords')
        coords = self.coords()[i][:,::N]
        print('res')
        res = cplt3d.generator_funcs.uniform_histogram(ax,coords.T,np.ones(len(coords[0])),cmap = cmap,verbose = self.verbose,**kwargs)
        if add_colorbar:
            plt.colorbar(ax = ax,mappable = colormaps.ScalarMappable(norm=res[0], cmap=cmap),shrink=shrink,pad = pad)
        return res
    
    @animated(autoparallel=False)
    @plotting
    def hist3d_tree(self,i,ax=None,N=10,add_colorbar = True,cmap = colormaps.get_cmap('viridis'),shrink = 0.8,pad = -0.2,**kwargs):
        print('coords')
        coords = self.coords()[i][:,::N]
        print('res')
        res = cplt3d.generator_funcs.tree_histogram(ax,coords.T,np.ones(len(coords[0])),cmap = cmap,verbose = self.verbose,**kwargs)
        if add_colorbar:
            plt.colorbar(ax = ax,mappable = colormaps.ScalarMappable(norm=res[0], cmap=cmap),shrink=shrink,pad = pad)
        return res
    
    def to_csv(self,i):
        df = pd.DataFrame(columns = ['x','y','z'])
        df['x'] = self.coords()[i][0]
        df['y'] = self.coords()[i][1]
        df['z'] = self.coords()[i][2]
        df.to_parquet('cosmological_simulation.parquet',index = False)