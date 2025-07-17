import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import time
import cplt3d.cube_plotter as cube_plotter
import cplt3d.helper as helper
from inspect import signature
import cplt3d.edgecolor_funcs as edgecolor_funcs

def volume_plot(fetch_coordinates):
    '''General plotting decorator. Needs to encapsulate something that will transform the points and values into box coordinates and their values.
    Function
    ----------
    fetch_coordinates : function: pts, vals -> x, y, z, dx, dy, dz, results
        A function that takes in pts, vals and returns the coordinates (x,y,z) and extents (dx,dy,dz) of the boxes in addition to their result (results).
    
    '''
    @helper.verbose_print    
    def wrap(ax,pts,vals,
            cmap = plt.get_cmap('viridis'),norm = colors.Normalize,vmin = None,vmax = None,
            _range = None,
            filled = 0.2,filled_invert = False,
            edgecolor_function = edgecolor_funcs.none_func,
            **kwargs):
        '''A general function that actually plots the voxels. Needs to wrap a way to generate the values and box positions from the input pts and vals.
        Parameters
        ----------
        ax : Axis
            The axis on which to plot the points
        pts : array of shape (N,3)
            The points to plot
        vals: array of shape (N,)
            The values to plot that are associated with each point
        cmap: Colormap
            The colormap to use while plotting, strongly recommended to have a gradient in alpha.
        norm: Function
            The normalization function to use. Colors are computed using cmap(NORM(x)) where NORM is the instantiation of norm with the correct vmax and vmin.
        vmin: float
            The minimum of the dynamic range. Set automatically if None.
        vmax: float
            The maximum of the dynamic range. Set automatically if None.
        _range: array of shape (3,2)
            The minima and maxima to generate bins within. Set automatically if None
        filled: float or function
            The method of making boxes invisible. If float the code removes the bottom(top) % of the bins. If function, must be of form f(x,y,z,dx,dy,dz,result,color), be vectorized, and return a np array of bools which are True if the box is plotted and False if not.
        filled_invert: bool
            If true, it inverts filled (so removes the top % of the data instead of bottom %). Has no impact if filled is a function.
        edgecolor_func: Function
            A function that takes in the facecolors and returns the edgecolors array. Can be useful if you want to shade your edge colors differently from your face colors (or change alpha)
        **kwargs:
            Other arguments for voxelize (and the polygon collection). facecolor and edgecolor are overriden.

        Returns
        -------
        Function
            The instantiation of the normalization function used.
        Poly3dCollection
            The polygons that are plotted.
        Array
            The x values of the boxes.
        Array
            The y values of the boxes.
        Array
            The z values of the boxes.
        Array
            The dx widths of the boxes.
        Array
            The dy widths of the boxes.
        Array
            The dz widths of the boxes.
        Array
            The results of the transformation into bins.
        '''

        assert pts.shape[1] == 3 and len(pts.shape) == 2, "Points must have shape (N,3)"
        # assert _range is None or (np.array(_range).shape == np.array([3,2])), "_range must be either None or an array of shape (3,2)"

        # prepare minimum and maximum values of bins
        if _range is None:
            min_x,max_x = np.min(pts[:,0]),np.max(pts[:,0])
            min_y,max_y = np.min(pts[:,1]),np.max(pts[:,1])
            min_z,max_z = np.min(pts[:,2]),np.max(pts[:,2])
        else:
            min_x,max_x = _range[0]
            min_y,max_y = _range[1]
            min_z,max_z = _range[2]

        print(f"finding bins within region\nx in ({min_x,max_x})\ny in ({min_y,max_y})\nz in ({min_z,max_z})")

        # getting the keywords for the coordinate fetching tool
        sig = signature(fetch_coordinates)
        keywords = [p.name for p in sig.parameters.values() if p.kind == p.POSITIONAL_OR_KEYWORD]
        kwargs_coords = {}
        for key in keywords:
            if key in kwargs.keys():
                kwargs_coords[key] = kwargs[key]
                kwargs.pop(key)

        # running the coordinate fetcher
        t = time.time()
        print(f"starting coordinate fetch using {fetch_coordinates.__name__}...")
        X,Y,Z,dX,dY,dZ,results = fetch_coordinates(pts,vals,min_x,max_x,min_y,max_y,min_z,max_z,**kwargs_coords)
        print('finished coordinate fetch in',time.time() - t,'seconds')
        print(f'identified {np.prod(X.shape),X.shape} boxes')

        print("computing things to use in plotter...")
        # computing vmin and vmax
        if vmin is None:
            vmin = np.min(results)
        if vmax is None:
            vmax = np.max(results)
        print('using vmin of',vmin)
        print('using vmax of',vmax)

        # Compute the color
        normalizer = norm(vmin = vmin,vmax = vmax)
        C = cmap(normalizer(results.flatten()))
        
        # Kill boxes we don't want visible in the final run        
        if filled is None:
            print("filled is None, so no boxes removed")
        elif type(filled) is float or type(filled) is int:
            if filled_invert:
                filled = 1-filled
            threshold = np.quantile(results,filled)
            if filled_invert:
                mask = results < threshold
                print(f'filled is float and filled_invert=True, so removing all boxes with result over {filled*100}% or result={threshold}')
                print(f'this eliminates {np.sum(~mask)} voxels, which is {np.sum(~mask)/len(X)}% of the total')
            else:
                mask = results > threshold
                print(f'filled is float and filled_invert=False, so removing all boxes with result under {filled*100}% or result={threshold}')
                print(f'this eliminates {np.sum(~mask)} voxels, which is {np.sum(~mask)/np.prod(X.shape)*100}% of the total')


            X,Y,Z = X[mask],Y[mask],Z[mask]
            dX,dY,dZ = dX[mask],dY[mask],dZ[mask]
            C = C[mask.flatten()]
        
        else:
            mask = filled(X.flatten(),Y.flatten(),Z.flatten(),dX.flatten(),dY.flatten(),dZ.flatten(),C)
        
            X,Y,Z = X[mask],Y[mask],Z[mask]
            dX,dY,dZ = dX[mask],dY[mask],dZ[mask]
            C = C[mask.flatten()]

        dX,dY,dZ = dX.flatten(),dY.flatten(),dZ.flatten()
        X,Y,Z = X.flatten(),Y.flatten(),Z.flatten()
        

        print('finished, ready to plot!')
        print('plotting...')
    
        res = cube_plotter.voxelize(ax,X,Y,Z,dX,dY,dZ,facecolors = C,edgecolors = edgecolor_function(C),**kwargs)
    
        print('finished plotting')

        # ax.plot(X,Y,Z,'bo')
        return normalizer,res,X,Y,Z,dX,dY,dZ,results

    return wrap