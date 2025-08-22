# cplt3d - A Package for Nifty 3d Plotting in Matplotlib

![til](/Examples/Gaussian/Images/7_Gaussian_Histogram.gif)

An example of the code's tree-mesh histogram generator applied to a sum of two 3d Gaussian distributions. Note how the code automatically increases the resolution in regions of higher density.

## Installation

This package is in development and does not have a Python Package Index yet. To install the package, close the repository and pip install manually. For example, in the directory where you wish to store `cplt3d`, run:

```
git clone https://github.com/benmckcohen/cplt3d.git
pip install ./cplt3d
rm -r cplt3d
```

## Use

### 3d Plotting

`cplt3d` has four functions to allow 3d voxel plotting.

#### cplt3d.uniform_histogram

`cplt3d.uniform_histogram` allows users to plot 3d histograms of datasets at a given resolution. 
**Parameters**
- ax : Axis
  - The axis on which to plot the points
- pts : array of shape (N,3)
  - The points to plot
- vals: array of shape (N,)
  - The values to plot that are associated with each point
- cmap: Colormap
  - The colormap to use while plotting, strongly recommended to have a gradient in alpha.
- norm: Function
  - The normalization function to use. Colors are computed using cmap(NORM(x)) where NORM is the instantiation of norm with the correct vmax and vmin.
- vmin: float
  - The minimum of the dynamic range. Set automatically if None.
- vmax: float
  - The maximum of the dynamic range. Set automatically if None.
- _range: array of shape (3,2)
  - The minima and maxima to generate bins within. Set automatically if None
- filled: float or function
  - The method of making boxes invisible. If float the code removes the bottom(top) % of the bins. If function, must be of form f(x,y,z,dx,dy,dz,result,color), be vectorized, and return a np array of bools which are True if the box is plotted and False if not
- filled_invert: bool
  - If true, it inverts filled (so removes the top % of the data instead of bottom %). Has no impact if filled is a function.
- edgecolor_function: Function
  - A function that takes in the facecolors and returns the edgecolors array. Can be useful if you want to shade your edge colors differently from your face colors (or change alpha)
- bins: int or list of ints of shape (3,)
  - The number of bins to use. If a list, sets [X bins, Y bins, Z bins]. If an int, generates that number on each axis. 
- statistic: string
  - The way to combine vals to generate the histogram. Input to scipy `binned_statistic_dd`. Note that cplt histograms generate densities, so whatever for statistic is computed, the number plotted is statistic/volume_of_bin.
- **kwargs:
  - Other arguments for voxelize (and the polygon collection). facecolor and edgecolor are overriden.

        Returns
        -------
        Function
            The instantiation of the normalization function used.
        Poly3dCollection
            The polygons that are plotted.
        Array
            The x values of the bins.
        Array
            The y values of the bins.
        Array
            The z values of the bins.
        Array
            The dx widths of the bins.
        Array
            The dy widths of the bins.
        Array
            The dz widths of the bins.
        Array
            Values of the bins.



### Animating 3d Plots
