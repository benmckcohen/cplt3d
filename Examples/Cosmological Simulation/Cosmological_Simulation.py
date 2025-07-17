if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import matplotlib.cm as colormaps
    import matplotlib.colors as colors
    import numpy as np
    import matplotlib.style as mplstyle
    import multiprocess
    import numpy as np
    import platform

    animate = True

    N_cpu = multiprocess.cpu_count()
    mplstyle.use('fast')

    
    from gadget_utils import GadgetSimulation
    from cplt3d.generator_funcs import spin_3d_plot



    sim = GadgetSimulation('./',verbose = True)

    

    
    linear_cmap = colormaps.get_cmap('inferno')
    linear_cmap_l = linear_cmap(np.linspace(0,1,linear_cmap.N))
    linear_cmap_l[:,-1] = np.linspace(0,1,linear_cmap.N)**1.1 # 1
    linear_cmap = colors.ListedColormap(linear_cmap_l)

    rainbow_cmap = colormaps.get_cmap("gist_rainbow")
    rainbow_cmap_l = rainbow_cmap(np.linspace(0,1,rainbow_cmap.N))
    rainbow_cmap_l[:,-1] = np.linspace(0,1,rainbow_cmap.N)**1.1 # 1
    rainbow_cmap = colors.ListedColormap(rainbow_cmap_l)

    
    # fig = plt.figure(figsize = (8,8),constrained_layout=True)
    fig = plt.figure(figsize = (5,5),constrained_layout=True)
    ax2 = fig.add_subplot(projection = '3d')

    # first plot a wireframe example to show how the method automatically hones in on high density/slope regions

    _,_,poly = sim.hist3d_tree(0,fig = fig, ax=ax2,cmap = rainbow_cmap,
                        add_colorbar = not animate, vmin = 1, vmax = 100,shrink=0.7,pad = 0.05,
                        min_resolution = 1,max_resolution = 4,focus='magnitude',
                        filled=None,edgecolor_function = lambda x:(0,0,0,0.6))#bins=0.00025,dist=dist

    ax2.set_xlim(0,100)
    ax2.set_ylim(0,100)
    ax2.set_zlim(0,100)

    

    # Animate the wireframe example
    animate = True
    already_set = False
    if animate:
        if platform.system() == "Darwin":
            multiprocess.set_start_method('spawn')
            already_set = True
        folder_in = '/Users/benjamincohen/Documents/Research/Inflation/Project/cplt3d/Examples/Cosmological Simulation/Wireframe'
        spin_3d_plot(fig,ax2,'rotating_wireframe',step=1,times = 1,fps = 30,
                    Animation_Generation_Folder=folder_in,
                    verbose = True,delete = False,merge=True,dpi = 300,constrained_layout = True)

    
    # run the real renderer!
    # fig = plt.figure(figsize = (8,8),constrained_layout=True)
    fig = plt.figure(figsize = (5,5),constrained_layout=True)
    ax2 = fig.add_subplot(projection = '3d')

    # next plot the real deal
    max_resolution = 9

    adj = 6.75000007e-01#0.1
    adj_frac = np.array([0,0,0.1,0.5,2,3,2,1,-10])/10
    dist = np.array([5.21540646e-08,
            4.17232517e-07,
            3.33786013e-06,
            2.67028811e-05,
            2.13623048e-04,
            1.70898439e-03,
            1.36718751e-02,
            1.09375001e-01,
            8.75000007e-01])
    dist = list(np.array(dist) + adj * adj_frac)

    _,_,poly = sim.hist3d_tree(0,fig = fig, ax=ax2,cmap = linear_cmap,
                        add_colorbar = not animate, vmin = 1, vmax = 3000,shrink=0.7,pad = 0.05,
                        min_resolution = 1,max_resolution = 9,focus='magnitude',
                        filled=0.2,bins=0.00025,dist=dist)

    ax2.set_xlim(0,100)
    ax2.set_ylim(0,100)
    ax2.set_zlim(0,100)


    ax2.set_proj_type('ortho')

    print('loading plot...')

    if not animate:
        plt.show()

        plt.close()

    # Animate the full plot
    
    if animate:
        if platform.system() == "Darwin" and not already_set:
            multiprocess.set_start_method('spawn')
        folder_in = '/Users/benjamincohen/Documents/Research/Inflation/Project/cplt3d/Examples/Cosmological Simulation/Simulation'
        spin_3d_plot(fig,ax2,'rotating_density',step=1,fps = 30,
                    parallel = True, Animation_Generation_Folder=folder_in,
                    verbose = True,delete = False,merge=True,dpi = 300,constrained_layout = True)


    