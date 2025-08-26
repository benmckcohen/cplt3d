import numpy as np
import mpl_toolkits.mplot3d.art3d as mpt3d
import time
import cplt3d.helper as helper


_vertices = np.array([
    [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],  # Bottom face
    [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]   # Top face
])



def get_face(v,a,b,c,d):
    '''Prepares a cube face from the vertices list.
    
    Given an array of vertices of the form 
        _vertices = np.array([
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],  # Bottom face
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]   # Top face
        ])
    and indices a, b, c, and d, this method returns an array of points corresponding to the polygon faces of the cube for plotting. 
    
    '''

    return np.swapaxes(np.array([v[:,a],v[:,b],v[:,c],v[:,d]]),0,1)

def create_cubes(x,y,z,dx,dy,dz):
    '''Given positions (in 1d arrays) and sizes (in 1d arrays), this method produces the polygons of the cube faces. The cube has extent (x,y,z) to (x+dx,y+dy,z+dz).

    Parameters
    ----------
    x : float or array with shape (N,)
        The x values of the lower corners of the cubes
    y : float or array with shape (N,)
        The y values of the lower corners of the cubes
    z : float or array with shape (N,)
        The z values of the lower corners of the cubes

    Returns
    -------
    top
        The top faces of the cubes
    bottom
        The bottom faces of the cubes
    front
        The front faces of the cubes
    back
        The back faces of the cubes
    left
        The left faces of the cubes
    right
        The right faces of the cubes
    '''

    # Produce the base vertices
    v = np.repeat(_vertices[np.newaxis,:],len(x),axis = 0)
    
    # Update the vertex positions
    s = np.array([dx,dy,dz]).T
    s = s[:,np.newaxis,:]
    t = np.array([x,y,z]).T
    t = t[:,np.newaxis,:]
    v = v*s + t
    
    # Getting the faces
    top = get_face(v,4,5,6,7)
    bottom = get_face(v,0,1,2,3)
    front = get_face(v,0,1,5,4)
    back = get_face(v,3,2,6,7)
    left = get_face(v,0,3,7,4)
    right = get_face(v,1,2,6,5)
    
    return top,bottom,front,back,left,right

@helper.verbose_print
def voxelize(ax,x,y,z,dx,dy,dz,**kwargs):
    '''This method adds voxels to the axis using a Poly3dCollection which is faster than the baked in matplotlib voxel objects. 
    
    Parameters
    ----------
    ax : matplotlib.Axis
        The axis on which to plot the voxels. **Important** make sure the axis is configured for 3d.
    x : float or array of shape (N,)
        The x positions of the voxels
    y : float or array of shape (N,)
        The y positions of the voxels
    z : float or array of shape (N,)
        The z positions of the voxels
    dx : float or array of shape (N,)
        The x extents of the voxels
    dy : float or array of shape (N,)
        The y extents of the voxels
    dz : float or array of shape (N,)
        The z extents of the voxels
    verbose : bool
        Flag whether to print run information

    Returns
    -------
    Poly3DCollection
        A collection of the polygons that make up the voxels.
    
    '''    
    
    assert (len(x) == len(y)) and (len(x) == len(z)) and (len(dx) == len(dy)) and (len(dx) == len(dz)) and (len(dx) == len(x)), \
    f"Input arrays must be the same length, yet I see lengths of\nx={len(x)}\ny={len(y)}\nz={len(z)}\ndx={len(dx)}\ndy={len(dy)}\ndz={len(dz)}"

    # Process facecolors and edgecolors as we will end up using 6 data points (faces) for each input value
    facecolors_included = 'facecolors' in kwargs.keys()
    edgecolors_included = 'edge' in kwargs.keys()

    if facecolors_included:
        try:
            iter(kwargs['facecolors'])
            facecolors_isiter = True
        except:
            facecolors_isiter = False
    if edgecolors_included:
        try:
            iter(kwargs['facecolors'])
            edgecolors_isiter = True
        except:
            edgecolors_isiter = False

    t = time.time()
    print('creating cubes...')
    
    top,bottom,front,back,left,right = create_cubes(x,y,z,dx,dy,dz)

    print("finished cubes in",time.time()-t,'seconds')
    
    t = time.time()
    print('creating polygons...')
    polygons = np.array(list(top) + list(bottom) + list(front) + list(back) + list(left) + list(right))
    
    # Adjusting matplotlib properties that are pointwise because we actually have 6 times the number of points
    if facecolors_included and facecolors_isiter:
        print("multiplying facecolors...")
        kwargs['facecolors'] = list(kwargs['facecolors']) * 6
    if edgecolors_included and edgecolors_isiter:
        print("multiplying edgecolors..")
        kwargs['edgecolors'] = list(kwargs['edgecolors']) * 6

    
    polygons = mpt3d.Poly3DCollection(polygons,**kwargs)
    
    print("finished polygons in",time.time()-t,'seconds')
    t = time.time()
    print('adding polygons to axis...')
    
    ax.add_collection3d(polygons)
    
    print('finished adding polies in',time.time()-t,'seconds')

    return polygons