# cplt3d - A Package for Nifty 3d Plotting in Matplotlib

![til](/Examples/Gaussian/Images/7_Gaussian_Histogram.gif)

`cplt3d` is a python package for creating nifty 3d plots in matplotlib. It was created for the purpose of creating 3d histograms, specifically with cosmological N-body simulations in mind. It is still a work in progress. 

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

#### `cplt3d.uniform_histogram`

`cplt3d.uniform_histogram` allows users to plot 3d histograms of datasets at a given resolution. Note that the rendering of high resolution 3d histograms can be very computationally expensive. Therefore, it is highly recommended to set filled to at least $0.1-0.2$ for non-trivial bin sizes.

**Parameters**
        -------------------------------------
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
- kwargs:
  - Other arguments for voxelize (and the polygon collection). facecolor and edgecolor are overriden.
**Returns**
- Function
  - The instantiation of the normalization function used.
- Poly3dCollection
  - The polygons that are plotted.
- Array
  - The x values of the bins.
- Array
  - The y values of the bins.
- Array
  - The z values of the bins.
- Array
  - The dx widths of the bins.
- Array
  - The dy widths of the bins.
- Array
  - The dz widths of the bins.
- Array
  - Values of the bins.

**Example**
        -------------------------------------
As an example, we can generate a histogram of a Gaussian distribution. First we can generate the distribution:
```python
# take some samples from a unimodal distribution
N = 10000
std = 0.3
np.random.seed(123456)
pts_samples = np.random.multivariate_normal(mean = [0,0,0],cov = [[std,0,0],[0,std,0],[0,0,std]],size = N)
GAUSSIAN_samples = np.ones(len(pts_samples))
```
Then we can plot
```python
# Import
from cplt3d.generator_funcs import uniform_histogram
```
```python
# Prepare the colormap
cmap = colormaps.get_cmap('viridis')
use_cmap = cmap(np.arange(cmap.N))
use_cmap[:,-1] = np.linspace(0,1,cmap.N)**1.7
use_cmap = colors.ListedColormap(use_cmap)
# Prepare the axis
fig = plt.figure(figsize = (3,3))
ax = fig.add_subplot(projection = '3d')
```
```python
# Actually plot
N_bins = 2**5
uniform_histogram(ax,pts_samples,GAUSSIAN_samples,bins = N_bins,filled = 0.3,cmap = use_cmap,verbose = False)
```
```python
# Save the figure
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.zaxis.labelpad=-3.5
fig.savefig('Images/3_Histogram-Uniform_Gaussian.png',dpi = 300)
```
Which leads to the figure:

![til](/Examples/Gaussian/Images/3_Histogram-Uniform_Gaussian.png)

#### `cplt3d.uniform_nearest_interpolator` and `cplt3d.uniform_nearest_interpolator`

`cplt3d.uniform_nearest_interpolator` interpolates a function using nearest neighbor interpolation to the voxels. This is faster than linear interpolation, which can be done through `cplt3d.linear_nearest_interpolator`. Both have the same input and output structure.

**Parameters**
        -------------------------------------
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
      - The method of making boxes invisible. If float the code removes the bottom(top) % of the bins. If function, must be of form f(x,y,z,dx,dy,dz,result,color), be vectorized, and return a np array of bools which are True if the box is plotted and False if not.
  - filled_invert: bool
      - If true, it inverts filled (so removes the top % of the data instead of bottom %). Has no impact if filled is a function.
  - edgecolor_function: Function
      - A function that takes in the facecolors and returns the edgecolors array. Can be useful if you want to shade your edge colors differently from your face colors (or change alpha)
  - bins: int or list of ints of shape (3,)
      - The number of bins to use. If a list, sets [X bins, Y bins, Z bins]. If an int, generates that number on each axis. 
  - kwargs:
      - Other arguments for voxelize (and the polygon collection). facecolor and edgecolor are overriden.

**Returns**
        -------------------------------------
        
  - Function
    - The instantiation of the normalization function used.
  - Poly3dCollection
    - The polygons that are plotted.
  - Array
    - The x values of the bins.
  - Array
    - The y values of the bins.
  - Array
    - The z values of the bins.
  - Array
    - The dx widths of the bins.
  - Array
    - The dy widths of the bins.
  - Array
    - The dz widths of the bins.
  - Array
    - Values of the bins.

**Example**
        -------------------------------------

As an example, we can generate a plot of a Gaussian distribution. Like before, we first generate the function:

```python
N = 32
X,Y,Z = np.linspace(-4,4,N),np.linspace(-4,4,N),np.linspace(-4,4,N)
X,Y,Z = np.meshgrid(X,Y,Z)
X = X.flatten()
Y = Y.flatten()
Z = Z.flatten()
pts = np.array([X,Y,Z]).T
GAUSSIAN = scipy.stats.multivariate_normal.pdf(pts,mean = [0,0,0],cov = [[1,0,0],[0,2,0],[0,0,3]])
```

Then we can plot it
```python
# Import
from cplt3d.generator_funcs import uniform_nearest_interpolator,uniform_linear_interpolator
```
```python
# Prepare colormap
cmap = colormaps.get_cmap('viridis')
use_cmap = cmap(np.arange(cmap.N))
use_cmap[:,-1] = np.linspace(0,1,cmap.N)**1.7
use_cmap = colors.ListedColormap(use_cmap)
# Prepare the axis
fig = plt.figure(figsize = (3,3))
ax = fig.add_subplot(projection = '3d')
```

```python
# Actually plot
uniform_nearest_interpolator(ax,pts,GAUSSIAN,bins = N_bins,filled = 0.3,cmap = use_cmap,verbose = True)
# or
uniform_linear_interpolator(ax,pts,GAUSSIAN,bins = N_bins,filled = 0.3,cmap = use_cmap,verbose = False)
```

```python
# Save the plots
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.zaxis.labelpad=-3.5
fig.savefig('Images/1_Nearest-Uniform_Gaussian.png',dpi = 300)
# or
fig.savefig('Images/2_Linear-Uniform_Gaussian.png',dpi = 300)
```

This produces plots that look like

![til](/Examples/Gaussian/Images/1_Nearest-Uniform_Gaussian.png)

#### `cplt3d.tree_histogram`

This is an experimental histogram plotter that uses different bin sizes to "zoom-in" on structure in the histogram. This allows very small bins to be plotted without massive cost in displaying the histogram when much of the structure is present in small regions of the plot. Note that this is an **experimental** feature and it is probably better to use `uniform_histogram` in most applications. Ways to improve computing `dist` automatically for the plots to look good in more situations are welcome!

  **Parameters**
        -------------------------------------
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
    - The method of making boxes invisible. If float the code removes the bottom(top) % of the bins. If function, must be of form f(x,y,z,dx,dy,dz,result,color), be vectorized, and return a np array of bools which are True if the box is plotted and False if not.
  - filled_invert: bool
    - If true, it inverts filled (so removes the top % of the data instead of bottom %). Has no impact if filled is a function.
  - edgecolor_function: Function
    - A function that takes in the facecolors and returns the edgecolors array. Can be useful if you want to shade your edge colors differently from your face colors (or change alpha)
  - min_resolution: int or None
    - The minimum resolution to use. Note this is log2(bins) at minimal bin size. If `None`, will compute the min_resolution as one higher resolution step than the minimal resolution to have an empty bin (or max_resolution + 1 if max_resolution is specified and smaller). 
  - max_resolution: int or None
    - The maximum resolution to use. Note that this is log2(bins) at maximal **bin** size. If `None`, will compute the resolution such that the bin size is as close to the average volume per particle (or min_resolution + 1 if min_resolution < that).
  - dist: list of floats or None
    - The distribution of percentages of volumes for each. Must sum to 1 and have a length = max_resolution - min_resolution + 1. For example [0.5,0.5] would cause the code to make as close to half the bins as possible level min_resolution bins and half the bins level max_resolution bins. 
    - If `"equivolume"` or `None` then will automatically generate the distribution as
    $$
    dist(L) = 1/N
    $$
    So each level takes up the same volume.
    - If `"sigmoid"` then will use a sigmoid distribution which is then normalized by the equivolume distribution. i.e. the un-normalized distribution takes the form:

    $$
    dist(L) = \frac{2}{1 + e^{2.1(L - \sqrt{\frac{1}{n}\sum_{\min}^{\max}n^2})}}\left (\frac{1}{2^{3L}}\right )^1.25
    $$

    The code distributes bins according to these volumes percentages as closely as possible. When there is overflow, it adds the overflow percentage to the next bin size.  

  <!-- - bins: int or float -->
    <!-- - int: The total number of bins. This will not be the exact number of bins plotted, but will be the ballpark number which is aimed for while the program allocates bin sizes using `dist`. -->
    - float: The percentage of volume taken up by the smallest bin size. Should be very small, especially for small bins, to prevent rendering issues. 
  - focus: string
    - How the code determines where to make smaller bins. The current two foci are 'slope' which extracts a percent with the largest slopes and 'magnitude' which extracts a percent with the largest magnitude.
  - kwargs:
    - Other arguments for voxelize (and the polygon collection). facecolor and edgecolor are overriden.

**Returns**
        -------------------------------------
        
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

**Example**
        -------------------------------------
As an example, we can once again plot a histogram of the Gaussian samples. Using the same setup as for `uniform_histogram`, we can plot simply with

```python
poly = tree_histogram(ax,pts_samples,GAUSSIAN_samples,cmap = use_cmap,
               filled=None,verbose = True,dist = 'sigmoid',
               min_resolution = None,max_resolution = None,edgecolor_function = lambda color:(0,0,0,0.01))
```
which gives us

![til](/Examples/Gaussian/Images/4_Histogram-Tree_Gaussian_1.png) 


### Animating 3d Plots

cplt3d comes equipped with two functions to make it easy to animate plots, especially 3d plots. The first is a general animation framework. 

#### `parallel_animate`

`parallel_animate` animates functions in parallel. Note that during animation creation this method creates a temperary directory and saves the frames as individual images. It then combines those frames and removes the file. The `merge` keyword controls whether or not the method actually combines the images. The `delete` keyword controls whether the method deletes the temperary folder.
  
  **Parameters**
          -------------------------------------
  - fig : Figure
      - The figure which contains the axes to animate.
  - func : Function
      - A function that takes in a frame from frames and performs any modifications necessary to the axes.
  - frames : list
      - A list of inputs to func.
  - result_name: str
      - The filename that the code saves the animation to.
  - out: str
      - The filetype of the animation. Default is `.gif`.
  - fps: int
      - The frames per second of the animation. 
  - parallel: bool
      - Whether or not to parallelize saving the frames of the animation.    
  - Animation_Generation_Folder: str
      - The folder to generate the animation frames without. If set to None will create folder with random integer name in the working directory.
  - delete: bool
      - Whether or not to delete the Animation_Generation_Folder and all of the frames within upon completion of the animation.
  - merge: bool
      - Whether or not to merge the frames into an actual animation. 
  - kwargs:
      - Arguments for the `savefig` operation on `fig`.

  **Returns**
        -------------------------------------
  - None


#### `spin_3d_plot`

`spin_3d_plot` spins a 3d plot or a set of 3d plots in a circle. This allows easier visualization. The method can run in parallel and uses `parallel_animate` under the hood, so the discussion above applies for this method too. 

  **Parameters**
  ----------
  - fig : Figure
      - The figure which contains the axes to rotate
  - axs : Axis or list
      - The axis on which to plot the points. If given a list it will animate all of the axes in the list (useful for animating a multi-cell plot). Make sure that all the axes are in the 3d projection!
  - result_name : str
      - The name of the file to save the animation to
  - times: float
      - The number of times to rotate. i.e. the function will rotate the plot int(360 * times) degrees.
  - step: int
      - The step size of the animation. i.e. the function will rotate the plot in steps of `step` degrees
  - parallel: bool
      - Whether or not to parallelize saving the frames of the animation.
  - axis: str (or 3d vector)
      - The axis to rotate around. Currently, only 'x', 'y', and 'z' are implemented. TODO implement 3d vector to rotate about arbitrary axis. 
  - fps: int
      - The frames per second of the animation. 
  - Animation_Generation_Folder: str
      - The folder to generate the animation frames without. If set to None will create folder with random integer name in the working directory.
  - delete: bool
      - Whether or not to delete the Animation_Generation_Folder and all of the frames within upon completion of the animation
  - merge: bool
      - Whether or not to merge the frames into an actual animation. 
  - kwargs:
      - Arguments for the general parallel animation generation function.
          - out: str
              The filetype of the animation. Default is `.gif`
          - kwargs
              What to input into the `fig` `savefig` method when saving frames (e.g. dpi)
  **Returns**
          -------------------------------------
  - None

**Example**
We can easily rotate subplots using this function, as seen in this example. We first create a bimodal distribution by combining two Gaussian distributions. 
```python
# take some samples from a bimodal distribution
N = 10000
std = 0.3
pos = 1
pts_samples_bi_1 = np.random.multivariate_normal(mean = [pos,pos,pos],cov = [[std,0,0],[0,std,0],[0,0,std]],size = N)
pts_samples_bi_2 = np.random.multivariate_normal(mean = [-pos,-pos,-pos],cov = [[std,0,0],[0,std,0],[0,0,std]],size = N)
pts_samples_bi = list(pts_samples_bi_1) + list(pts_samples_bi_2)
GAUSSIAN_samples_bi = np.ones(len(pts_samples_bi))
```
Then we can plot it as we have before
```python
fig = plt.figure(figsize = (6*2,6))
ax = fig.add_subplot(121,projection = '3d')
ax2 = fig.add_subplot(122,projection = '3d')

cmap = colormaps.get_cmap('inferno')
use_cmap = cmap(np.arange(cmap.N))
use_cmap[:,-1] = np.linspace(0,1,cmap.N)**2
use_cmap = colors.ListedColormap(use_cmap)
d = 3

poly = tree_histogram(ax,np.array(pts_samples_bi),np.array(GAUSSIAN_samples_bi),cmap = use_cmap,
               filled=None,verbose = False,
               _range = [[-d,d],[-d,d],[-d,d]],
               min_resolution = 1,max_resolution = 5,
               vmin = 1,norm = LogNorm,edgecolor_function = lambda x:(0,0,0,0.2),linewidths = 0.4)

poly = tree_histogram(ax2,np.array(pts_samples_bi),np.array(GAUSSIAN_samples_bi),cmap = use_cmap,
               filled=None,verbose = False,
               _range = [[-d,d],[-d,d],[-d,d]],
               min_resolution = 1,max_resolution = 5,
               vmin = 1,norm = LogNorm)


ax.set_xlim(-d,d)
ax.set_ylim(-d,d)
ax.set_zlim(-d,d)
fig.tight_layout()
fig.savefig('Images/6_Histogram-Tree Dual-Gaussian.png',dpi = 500)
```
Now animating this plot is as simple as passing it to `spin_3d_plot` as follows
```python
if __name__ == '__main__':
    folder_in = './Images'
    spin_3d_plot(fig,[ax,ax2],'Images/7_Gaussian_Histogram',step=1,merge=True,delete=True,fps = 15,
                parallel = True,verbose = True,Animation_Generation_Folder=folder_in,dpi = 500)
```
