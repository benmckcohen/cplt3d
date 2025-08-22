# Example python file to generate tree_histogram of cosmological simulation
# Set animate = True to animate both the wireframe and the histogram that are generated

animate = False

if __name__ == '__main__':
    # Import stuff
    import matplotlib.pyplot as plt
    import matplotlib.cm as colormaps
    import matplotlib.colors as colors
    import numpy as np
    import matplotlib.style as mplstyle
    import multiprocess
    import numpy as np
    import platform
    import pandas as pd
    from matplotlib.colors import LogNorm
    from cplt3d.generator_funcs import tree_histogram,spin_3d_plot
    N_cpu = multiprocess.cpu_count()
    mplstyle.use('fast')

    # ===============================================================================================

    # Prepare colormaps
    linear_cmap = colormaps.get_cmap('inferno')
    linear_cmap_l = linear_cmap(np.linspace(0,1,linear_cmap.N))
    linear_cmap_l[:,-1] = np.linspace(0,1,linear_cmap.N)**1.1 # 1
    linear_cmap = colors.ListedColormap(linear_cmap_l)

    rainbow_cmap = colormaps.get_cmap("gist_rainbow")
    rainbow_cmap_l = rainbow_cmap(np.linspace(0,1,rainbow_cmap.N))
    rainbow_cmap_l[:,-1] = np.linspace(0,1,rainbow_cmap.N)**1.1 # 1
    rainbow_cmap = colors.ListedColormap(rainbow_cmap_l)
    #===============================================================================================

    # Load Data
    print('loading data...')
    df = pd.read_parquet('./cosmological_simulation.parquet')
    # df = df.iloc[::5]
    # df.to_parquet("cosmological_simulation.parquet")
    coords = np.array([list(df['x']),list(df['y']),list(df['z'])]).T
    vals = np.ones(len(coords))
    print('finished loading data')
    #===============================================================================================

    # Prepare figure for the wireframe
    fig = plt.figure(figsize = (5,5),constrained_layout=True)
    ax2 = fig.add_subplot(projection = '3d')

    # Plot a wireframe example to show how the method automatically hones in on high density/slope regions
    res = tree_histogram(ax2,coords,vals,cmap = rainbow_cmap,verbose = True, 
                        vmin = 1, vmax = 100,
                        min_resolution = 1,max_resolution = 4,focus='magnitude',
                        filled=None,edgecolor_function = lambda x:(0,0,0,0.6))
    norm = res[0]
    poly = res[2]
    # if not animate:
    #     plt.colorbar(ax = ax2,mappable = colormaps.ScalarMappable(norm=norm, cmap=rainbow_cmap),shrink=0.7,pad = 0.05)
    

    # Set axis limits for the wireframe example
    ax2.set_xlim(0,100)
    ax2.set_ylim(0,100)
    ax2.set_zlim(0,100)

    # Animate the wireframe example
    already_set = False
    if not animate:
        plt.savefig("./wireframe.png",dpi = 300)
    else:
        if platform.system() == "Darwin":
            multiprocess.set_start_method('spawn')
            already_set = True
        folder_in = './Wireframe'
        spin_3d_plot(fig,ax2,'rotating_wireframe',step=1,times = 1,fps = 30,
                    Animation_Generation_Folder=folder_in,
                    verbose = True,delete = False,merge=True,dpi = 300,constrained_layout = True)

    print('')
    print('')
    print('=' * 30)
    print('='*7,'Finished Wireframe, Starting Density','=' * 7)
    print('=' * 30)
    print('')
    print('')
    
    # Prepare for the full histogram
    fig = plt.figure(figsize = (5,5),constrained_layout=True)
    ax2 = fig.add_subplot(projection = '3d')

    # Set the maximum resolution
    max_resolution = 9

    # adj = 6.75000007e-01#0.1
    # adj_frac = np.array([0,0,0.1,0.5,2,3,2,1,-10])/10
    # dist = np.array([5.21540646e-08,
    #         4.17232517e-07,
    #         3.33786013e-06,
    #         2.67028811e-05,
    #         2.13623048e-04,
    #         1.70898439e-03,
    #         1.36718751e-02,
    #         1.09375001e-01,
    #         8.75000007e-01])
    # dist = list(np.array(dist) + adj * adj_frac)
    
    # Plot the full histogram
    res = tree_histogram(ax2,coords,vals,cmap = linear_cmap,verbose = True, #norm = LogNorm,
                                        vmin = 1, vmax = 6000,dist = 'sigmoid',
                                        min_resolution = 1,max_resolution = 9,focus='magnitude',
                                        filled=0.2,bins=50000)
    norm = res[0]
    poly = res[2]
    # if not animate:
    #     plt.colorbar(ax = ax2,mappable = colormaps.ScalarMappable(norm=norm, cmap=linear_cmap),shrink=0.7,pad = 0.05)


    # Update the axis limits and the projection
    ax2.set_xlim(0,100)
    ax2.set_ylim(0,100)
    ax2.set_zlim(0,100)
    ax2.set_proj_type('ortho')

    print('loading plot...')

    # Show or Animate the plot
    if not animate:
        # plt.show()
        plt.savefig("./density.png",dpi = 300)
        plt.close()
    
    else:
        if platform.system() == "Darwin" and not already_set:
            multiprocess.set_start_method('spawn')
        folder_in = './Simulation'
        spin_3d_plot(fig,ax2,'rotating_density',step=1,fps = 30,
                    parallel = True, Animation_Generation_Folder=folder_in,
                    verbose = True,delete = False,merge=True,dpi = 300,constrained_layout = True)


    