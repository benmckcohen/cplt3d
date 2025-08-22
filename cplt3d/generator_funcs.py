import numpy as np
import scipy.interpolate as interpolate
import multiprocess
import scipy
import cplt3d.plotting as plotting
import time
import __main__ as main
interactive = not hasattr(main, '__file__')
if interactive:
    import tqdm.notebook as tqdm
else:
    import tqdm
import os
import cplt3d.helper as helper
import multiprocess
import matplotlib.animation 
import subprocess
import sys
import imageio

# print("cplt3d using",N_cpu,'cpus if necessary.')

def uniform_generator(interpolator):
    '''A decorator of interpolation functions to use a uniform distribution of bins
    interpolator: function that takes in (pts,vals,X,Y,Z,dX,dY,dZ) and returns results as the thing to be plotted for each bin
    '''
    
    def wrap(pts,vals,min_x,max_x,min_y,max_y,min_z,max_z,bins=25,**kwargs):
        if type(bins) == int:
            N_x,N_y,N_z = bins,bins,bins
        else:
            N_x,N_y,N_z = bins[0],bins[1],bins[2]
        N_x,N_y,N_z = N_x+1,N_y+1,N_z+1
        X,Y,Z = np.linspace(min_x,max_x,N_x),np.linspace(min_y,max_y,N_y),np.linspace(min_z,max_z,N_z)
        dX,dY,dZ = X[1]-X[0],Y[1]-Y[0],Z[1]-Z[0]
        
        X,Y,Z = np.meshgrid(X[:-1],Y[:-1],Z[:-1])
        dX,dY,dZ = dX * np.ones(X.shape),dY * np.ones(Y.shape),dZ * np.ones(Z.shape)
        return X,Y,Z,dX,dY,dZ,interpolator(pts,vals,X,Y,Z,dX,dY,dZ,**kwargs)
    
    return wrap

def _convolve_one_level(X,Y,Z,dX,dY,dZ,convolving_function):
    new_X,new_Y,new_Z = X[0::2,0::2,0::2],Y[0::2,0::2,0::2],Z[0::2,0::2,0::2]
    new_dX,new_dY,new_dZ = dX[0::2,0::2,0::2]*2,dY[0::2,0::2,0::2]*2,dZ[0::2,0::2,0::2]*2
    new_results = convolving_function(new_X,new_Y,new_Z,new_dX,new_dY,new_dZ)
    return new_X,new_Y,new_Z,new_dX,new_dY,new_dZ,new_results

def expand(array):
    new_array = np.zeros(np.array(array.shape) * 2)
    
    new_array[0:-1:2,0:-1:2,0:-1:2] = array
    new_array[0:-1:2,1::2,0:-1:2] = array
    new_array[0:-1:2,0:-1:2,1::2] = array
    new_array[0:-1:2,1::2,1::2] = array
    
    new_array[1::2,0:-1:2,0:-1:2] = array
    new_array[1::2,1::2,0:-1:2] = array
    new_array[1::2,0:-1:2,1::2] = array
    new_array[1::2,1::2,1::2] = array

    return new_array

def and_contract(array):
    return array[0:-1:2,0:-1:2,0:-1:2] &\
    array[0:-1:2,1::2,0:-1:2] &\
    array[0:-1:2,0:-1:2,1::2] &\
    array[0:-1:2,1::2,1::2] &\
    array[1::2,0:-1:2,0:-1:2] &\
    array[1::2,1::2,0:-1:2] &\
    array[1::2,0:-1:2,1::2] &\
    array[1::2,1::2,1::2]

def or_contract(array):
    return array[0:-1:2,0:-1:2,0:-1:2] |\
    array[0:-1:2,1::2,0:-1:2] |\
    array[0:-1:2,0:-1:2,1::2] |\
    array[0:-1:2,1::2,1::2] |\
    array[1::2,0:-1:2,0:-1:2] |\
    array[1::2,1::2,0:-1:2] |\
    array[1::2,0:-1:2,1::2] |\
    array[1::2,1::2,1::2]

def do_I_convolve(dX,dY,dZ,results,bound,depth):
    focus = bound[0]
    bound = bound[1]

    print('checking convolution')
    print(f'using focus "{focus}"')
    # print(f'at depth {depth}')

    if bound[depth - 1] == '-':
        print('below minimum depth, exiting')
        return np.zeros(results.shape).astype(bool)
    
    if focus == 'slope':
        results_x_1,results_x_2 = results[0:-1:2,0:-1:2,0:-1:2],results[0:-1:2,1::2,0:-1:2]
        results_y_2 = results[1::2,0:-1:2,0:-1:2]
        results_z_2 = results[0:-1:2,0:-1:2,1::2]
        ddx,ddy,ddz = np.abs(results_x_2-results_x_1), np.abs(results_y_2-results_x_1), np.abs(results_z_2-results_x_1)
        dd = np.max(np.array([ddx,ddy,ddz]),axis=0)
    elif focus == 'magnitude':
        results_x_1= results[0:-1:2,0:-1:2,0:-1:2]
        dd = np.abs(results_x_1)
    else:
        raise NotImplementedError(f"There is no current implementation for focus = {focus}")

    # r_flat = dd
    r_flat_ind = (np.argsort(np.argsort(dd.flatten()))/len(dd)**3).reshape(dd.shape)
    # print(r_flat_ind)
    in_need_of_convolve_c = r_flat_ind < 1-bound[depth-1]
    
    in_need_of_convolve = expand(in_need_of_convolve_c)

    return in_need_of_convolve.astype(bool)

def _convolve_up_dict(X,Y,Z,dX,dY,dZ,results,convolving_function,bound,depth,min_depth = 0):
    '''Takes in the current level's boxes and returns lists of X,Y,Z,dX,dY,dZ for plotting.'''
    print('-' * 10 + f' Depth = {depth} ' + '-' * 10)
    in_need_of_convolve = do_I_convolve(dX,dY,dZ,results,bound,depth)
    above = {str(len(in_need_of_convolve)):[X,Y,Z,dX,dY,dZ,results,in_need_of_convolve]}
    if len(X) == 1+min_depth:
        print("depth=min reached, terminating")
        return above
    elif not np.sum(in_need_of_convolve):
        print("there is no more need to convolve, terminating")
        return above
    else:
        print(f'convolving {np.sum(in_need_of_convolve)} voxels, {np.sum(in_need_of_convolve)/len(in_need_of_convolve)**3}% of total, {np.sum(~in_need_of_convolve)} left over')
        X2,Y2,Z2,dX2,dY2,dZ2,results2 = _convolve_one_level(X,Y,Z,dX,dY,dZ,convolving_function)
        print(f'{len(X)}^3 -> {len(X2)}^3')
        prev = _convolve_up_dict(X2,Y2,Z2,dX2,dY2,dZ2,results2,convolving_function,bound,depth - 1)
        prev.update(above)
        return prev

def slice_levels(level_index,i,j,k,conv,levels):
    cur_level = levels[level_index]
    if level_index == len(levels)-1:
        return [[conv[cur_level][a][i,j,k] for a in range(7)]]
    # print(level_index,levels)
    next_level = levels[level_index+1]
    next_bools = conv[next_level][-1]
    next_bools = next_bools[i:i+2,j:j+2,k:k+2]
    all_true_on_next = np.prod(next_bools)
    if all_true_on_next:
        return [[conv[cur_level][a][i,j,k] for a in range(7)]]
    else:
        return  slice_levels(level_index+1,i*2,j*2,k*2,conv,levels)\
                + slice_levels(level_index+1,i*2,j*2+1,k*2,conv,levels)\
                + slice_levels(level_index+1,i*2,j*2,k*2+1,conv,levels)\
                + slice_levels(level_index+1,i*2,j*2+1,k*2+1,conv,levels)\
                + slice_levels(level_index+1,i*2+1,j*2,k*2,conv,levels)\
                + slice_levels(level_index+1,i*2+1,j*2+1,k*2,conv,levels)\
                + slice_levels(level_index+1,i*2+1,j*2,k*2+1,conv,levels)\
                + slice_levels(level_index+1,i*2+1,j*2+1,k*2+1,conv,levels)

def _convolve_up(X,Y,Z,dX,dY,dZ,results,convolving_function,bound,depth,min_depth = 0):
    conv = _convolve_up_dict(X,Y,Z,dX,dY,dZ,results,convolving_function,bound,depth,min_depth = min_depth)
    print('-' * 33)
    
    print("loading boolean maps...")
    levels = list(conv.keys())
    lowest_level = levels[0]
    highest_level = levels[-1]
    
    voxels = []
    bool_maps = []
    for level in levels:
        bool_maps.append(conv[level][-1])
    print('finished loading boolean maps')
    bool_maps.reverse()
    print(f'processing boolean maps... ({[bool_map.shape[0] for bool_map in bool_maps]})')
    bool_maps_processed = []
    convolved = None
    i = 0
    for bool_map in bool_maps:
        print(f'processing {len(bool_map)} level bool_map')
        if i==0:
            # at the lowest level we just want to keep anything that isn't scheduled to be convolved
            keep = ~bool_map
            used = keep
        else:
            # now we have to worry about a few effects. 
            # generally, we want to keep any component of boxes that have had a component be kept
            # as well as any boxes where all of the components are not scheduled for convolution
            # this translates to our only wanting to convolve sets of 8 boxes where NONE of the boxes are either
            #  (1) not scheduled for convolution on this level 
            #  (2) have been included before

            # the boxes scheduled for convolution are bool_map, 
            # so the regions where >=one of the 8 boxes is not scheduled are
            # not_all_scheduled = or_contract(~bool_map)
            # and thus the boxes are
            # not_all_scheduled = expand(not_all_scheduled)
            # but this is of course the same as the original bool_map because the convolution decision is computed for sets of 8 and then itself expanded
            # assert np.array_equal(not_all_scheduled,bool_map)

            # if any of a set of 8 have been included before, we want do not want to convolve. used = True if a bin has been included before, so the regions that have a bin which has been included before are
            any_included = or_contract(used).astype(bool)
            # so the boxes are
            any_included = expand(any_included).astype(bool)

            # thus the boxes to convolve are
            # to_convolve = bool_map & ~any_included
            # so the boxes to keep are
            keep = ~bool_map | any_included
            # note that keep is expanded because both bool_map and any_included are expanded, so we may freely contract it
            
            # however, we must be careful not to include boxes which have already been included!
            keep = keep & ~used

            # then we must update the boxes we have used
            used = used | keep

        # print(f'actually keeping {np.sum(keep)}')
        bool_maps_processed.append(keep)
        
        
        used = and_contract(used)
        print(f'    - plotted {np.sum(used)/len(used)**3*100}% of region')
        # print(used)
        i+=1
        

    bool_maps_processed.reverse()
    
    print('finished processing boolean maps')
    print(f'bin distribution: {[np.sum(boolean_map) for boolean_map in bool_maps_processed]}')
    
    print('preparing voxels...')
    voxels = []
    r = len(levels)
    
    for i in range(r):
        level = levels[i]
        if level not in list(conv.keys()):
            print(f"didn't find {level} in keys, probably skipped due to no need for convolution.")
        level_X,level_Y,level_Z,level_dX,level_dY,level_dZ,level_results,_ = conv[level]
        bool_map = bool_maps_processed[i]
        X_include,Y_include,Z_include,dX_include,dY_include,dZ_include,results_include = level_X[bool_map],level_Y[bool_map],level_Z[bool_map],level_dX[bool_map],level_dY[bool_map],level_dZ[bool_map],level_results[bool_map]
        voxels+=list(np.array([X_include,Y_include,Z_include,dX_include,dY_include,dZ_include,results_include]).T)

    voxels = np.array(voxels).T
    voxels = [list(x) for x in voxels]
    if len(voxels) == 0:
        voxels = [[],[],[],[],[],[],[]]
    
    print('finished preparing voxels')

    return voxels + [conv[highest_level][-1]]
    

def compute_dist(pts,vals,min_x,max_x,min_y,max_y,min_z,max_z,
             min_resolution = 0,max_resolution = 10,dist = None):
    if dist == 'equivolume' or dist is None:
        print('using equal area dist calculation...')
        N = max_resolution - min_resolution + 1
        return [1/N] * N,False
    elif dist == 'sigmoid':
        print('using sigmoid dist calculation')
        dist = np.arange(min_resolution,max_resolution+1,1)
        center = np.sqrt(np.sum(dist**2)/len(dist))
        print(f'centering distribution at {center}')
        sigmoid = 2/(1+np.exp(-2.1 * (dist - center)))
        equivol = 1*2**(3 * dist)
        dist = sigmoid/equivol**1.25
        dist = dist/np.sum(dist)
        return dist,True
    else:
        raise NotImplementedError(f"dist={dist} not implemented")
    

def compute_min_resolution(pts,vals,min_x,max_x,min_y,max_y,min_z,max_z,
             min_resolution = 0,max_resolution = 10,dist = None, bins = 0.1,focus = 'slope',
             **kwargs):
    print('finding largest size with one region having zero particles...')
    finished = False
    min_resolution = 0
    while not finished:
        min_resolution += 1
        N = 2**min_resolution
        X,Y,Z,dX,dY,dZ,res = uniform_generator(_histogram)(pts,vals,min_x,max_x,min_y,max_y,min_z,max_z,bins=N,**kwargs)
        finished = (0 in res)
        if not finished:
            print(f'min_resolution={min_resolution} doesn\'t do it, continuing')
    print(f'found minimal gets-zero min_resolution={min_resolution}')
    # print(f'adding one to make sure we capture everything... Now is {min_resolution+1}')
    # min_resolution += 1
    if not max_resolution is None:
        assert min_resolution < max_resolution,'Your max resolution is smaller than the generated minimum resolution, your data are too uniform!'
    return min_resolution

def compute_max_resolution(pts,vals,min_x,max_x,min_y,max_y,min_z,max_z,
             min_resolution = 0,max_resolution = 10,dist = None, bins = 0.1,focus = 'slope',
             **kwargs):
    print('finding average volume per particle...')
    number_of_particles = len(vals) # or number of max bins so if perfectly uniform get 1 particle per bin
    max_resolution = np.log2(number_of_particles)/3
    max_resolution = int(max_resolution) + 1
    volume = (max_x - min_x) * (max_y - min_y) * (max_z - min_z)
    print(f'found max_resolution={max_resolution} - volume={volume/2**(3 * max_resolution)} units^3')
    if max_resolution <= min_resolution:
        print(f'this is the same as or less than the min_resolution... setting to 1 more than min_resolution')
        max_resolution = min_resolution+1
    return max_resolution

def compute_bound_from_dist(dist,optimize,bins,max_resolution):
    bound = []
    N_remains = bins
    total_boxes = 0
    percent_to_makeup = 0
    for i in range(0,len(dist)):
        # Percent of volume to keep
        volume_percent = dist[len(dist)-1-i]
        # Level and number of bins
        level = max_resolution - i
        level_num = 2**(3*level)
        # The number of bins to ideally keep, may be decimal
        number_making_volume_percent = level_num * min(volume_percent + percent_to_makeup,1)
        # The number to actually keep, must be int
        number_keep = int(number_making_volume_percent)
        total_boxes+=number_keep
        to_bind = number_keep/level_num
        print(f'at level {level} we ideally keep {number_keep} ({volume_percent*100}% for this level, {percent_to_makeup*100}% from previous levels)')
        print(f'    -> so of {level_num} actually keep {number_keep}, or {to_bind*100}%')
        # The percent that is left over. e.g. if level 1, then 8 bins. If percent was 2, then 0.16 ideal, 0 actual, and 0.16/8 = 0.02 (2%) to makeup
        percent_to_makeup = (number_making_volume_percent - number_keep)/level_num
        print(f'    -> so we need to makeup {percent_to_makeup * 100}% of the volume')
        if total_boxes < bins:
            print(f'    -> currently at {total_boxes} boxes < upper limit {bins}')
        else:
            print(f'    -> currently at {total_boxes} boxes >= upper limit {bins}, forcing this layer to be 100% bound')
            to_bind = 1
        # then bound represents the percent that we keep
        bound.append(to_bind)
    bound.reverse()
    return bound

def tree_generator(interpolator):
    '''A decorator of interpolation functions to use a tree-based distribution of bins
    interpolator: function that takes in (pts,vals,X,Y,Z,dX,dY,dZ) and returns results as the thing to be plotted for each bin
    '''

    def wrap(pts,vals,min_x,max_x,min_y,max_y,min_z,max_z,
             min_resolution = None,max_resolution = None,dist = None, bins = 10000,focus = 'slope',
             **kwargs):
        
        print(f"running from Resolution {min_resolution} to {max_resolution}")
        
        if min_resolution is None:
            print('Found min_resolution is None, recomputing min_resolution...')
            min_resolution = compute_min_resolution(pts,vals,min_x,max_x,min_y,max_y,min_z,max_z,
                            min_resolution = min_resolution,max_resolution = max_resolution,dist = dist, bins = bins,focus = focus,
                            **kwargs)
            print(f'using min_resolution={min_resolution}')

        if max_resolution is None:
            print('Found max_resolution is None, recomputing max_resolution...')
            max_resolution = compute_max_resolution(pts,vals,min_x,max_x,min_y,max_y,min_z,max_z,
                            min_resolution = min_resolution,max_resolution = max_resolution,dist = dist, bins = bins,focus = focus,
                            **kwargs)
            print(f'using max_resolution={max_resolution}')
        print('dist',dist,type(dist))
        if dist is None or type(dist) is str:
            print(f'Found dist is {dist}, recomputing dist...')
            dist,optimize = compute_dist(pts,vals,min_x,max_x,min_y,max_y,min_z,max_z,
             min_resolution = min_resolution,max_resolution = max_resolution,dist = dist)
            print(f'using dist={dist}, sum = {np.sum(dist)}')
        else:
            optimize = False

        if type(bins) == int:
            print("set bins using integer, implying (approximate) total bin number")
            print(f'effective nside is {bins**(1/3)}')
            print(f"    - l={max_resolution}")
            print(f"    - highest resolution n={2**(max_resolution)}")
        if type(bins) == float:
            print(f"set bins using float, implying total fraction of volume occupied by smallest bin size")
            print(f"    - l={max_resolution}")
            print(f"    - highest resolution n={2**(max_resolution)}")
            assert 0<bins and bins < 1, 'set bins using integer, needs to be fraction of volume occupied by smallest bin size so between 0 and 1'
            bins = int(bins * 2**(3 * max_resolution)/dist[-1])
            print(f"found {bins} total voxels")

        assert len(dist) == (max_resolution - min_resolution) + 1,f'Distribution is too long, len(dist) should = max_resolution-min_resolution + 1. Instead, len(dist) = {len(dist)} != {max_resolution - min_resolution + 1} (max_resolution={max_resolution},min_resolution={min_resolution})'

        bound = compute_bound_from_dist(dist,optimize,bins,max_resolution)
        
        # for i in range(len(dist)):
        #     percent = dist[i]
        #     level = min_resolution + i
        #     level_num = 2**(3*level)
        #     # compute how many we ideally would like to keep in this bin
        #     num_ideal = N_remains * percent\
        #                 + extra_bins
        #     # but of course we can't go above the level_num
        #     print(f'at level {level} we ideally keep {num_ideal} > {N_remains * percent}')
        #     num_actual = min(num_ideal,level_num)
        #     print(f'-> so of {level_num} actually keep {num_actual}, or {num_actual/level_num*100}%')
        #     # then bound represents the percent that we keep
        #     bound.append(num_actual/level_num)
        #     extra_bins = (num_ideal - num_actual)
        #     N_remains = N_remains - num_actual

        
        # for i in range(min_resolution):
        #     bound[i] = '-'
        
        bound = ['-'] * min_resolution + bound

        bound = [focus,bound]

        print(f'convolving specifications:')
        print(f'    - bin target: {bins}')
        print(f'    - focus: {bound[0]}')
        print(f'    - distribution: {dist}')
        print(f'    - % to keep: {bound[1]}')
        
        # generate tree
    
        level_generator = lambda min_x,max_x,min_y,max_y,min_z,max_z,bins:uniform_generator(interpolator)(pts,vals,min_x,max_x,min_y,max_y,min_z,max_z,bins=bins,**kwargs)
        
        convolving_function = lambda X,Y,Z,dX,dY,dZ:interpolator(pts,vals,X,Y,Z,dX,dY,dZ,**kwargs)

        X_all,Y_all,Z_all,dX_all,dY_all,dZ_all,results_all = [],[],[],[],[],[],[]

        bins = 2**max_resolution
        boxes_to_expand = [
            [min_x,min_y,min_z,max_x-min_x,max_y-min_y,max_z-min_z,bins]
            ]

        def step(box_to_expand):
            # print(highest)
            min_x,min_y,min_z,dx,dy,dz,bins = box_to_expand
            max_x,max_y,max_z = min_x + dx, min_y + dy, min_z + dz
            
            X,Y,Z,dX,dY,dZ,results = level_generator(min_x,max_x,min_y,max_y,min_z,max_z,bins)
            
            X_new,Y_new,Z_new,dX_new,dY_new,dZ_new,results_new,convolved = _convolve_up(X,Y,Z,dX,dY,dZ,results,convolving_function,bound,max_resolution,min_depth=min_resolution)
            if np.any(~convolved):
                X_todo,Y_todo,Z_todo,dX_todo,dY_todo,dZ_todo = X[~convolved],Y[~convolved],Z[~convolved],dX[~convolved],dY[~convolved],dZ[~convolved]
                X_todo,Y_todo,Z_todo,dX_todo,dY_todo,dZ_todo = X_todo.flatten(),Y_todo.flatten(),Z_todo.flatten(),dX_todo.flatten(),dY_todo.flatten(),dZ_todo.flatten()
                new_boxes_to_expand = [X_todo,Y_todo,Z_todo,dX_todo,dY_todo,dZ_todo]
            else:
                new_boxes_to_expand = []
            
            new_voxels = [X_new,Y_new,Z_new,dX_new,dY_new,dZ_new,results_new]
            
            return new_voxels,new_boxes_to_expand

        print(f'expanding {len(boxes_to_expand)} regions...')
        
        # generate next step
        processing_result = []
        
        s = lambda boxes_to_expand:step(boxes_to_expand)
        # if parallel:
        #     print('computing...')
        #     with multiprocess.Pool(N_cpu) as p:
        #         processing_result = p.map(s,boxes_to_expand)
        # else:
        processing_result = [s(x) for x in boxes_to_expand]
        
        for voxels,_ in processing_result:
            # add the new voxels to the list of things to plot
            X_new,Y_new,Z_new,dX_new,dY_new,dZ_new,results_new = voxels
            
            X_all += X_new
            Y_all += Y_new
            Z_all += Z_new
            dX_all += dX_new
            dY_all += dY_new
            dZ_all += dZ_new
            results_all += results_new

        print(f"completed processing.")
        print(f"current total boxes is {len(results_all)}, decreased to {len(results_all)/2**(3*max_resolution)*100}%")
        
        return np.array(X_all),np.array(Y_all),np.array(Z_all),np.array(dX_all),np.array(dY_all),np.array(dZ_all),np.array(results_all)

    return wrap
            

def _histogram(pts,vals,X,Y,Z,dX,dY,dZ,statistic = 'sum'):
    _X,_Y,_Z = list(X[0,:,0]),list(Y[:,0,0]),list(Z[0,0,:])
    _X.append(_X[-1] + dX[0,0,0])
    _Y.append(_Y[-1] + dY[0,0,0])
    _Z.append(_Z[-1] + dZ[0,0,0])
    results,_,_ = scipy.stats.binned_statistic_dd(pts,vals,statistic = statistic,bins = [_X,_Y,_Z])
    results.reshape(X.shape)
    results = np.swapaxes(results,0,1)
    return results

@plotting.volume_plot
@uniform_generator
def uniform_histogram(pts,vals,X,Y,Z,dX,dY,dZ,statistic = 'sum'):
    '''A Method for Plotting a 3D Histogram
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
        edgecolor_function: Function
            A function that takes in the facecolors and returns the edgecolors array. Can be useful if you want to shade your edge colors differently from your face colors (or change alpha)
        bins: int or list of ints of shape (3,)
            The number of bins to use. If a list, sets [X bins, Y bins, Z bins]. If an int, generates that number on each axis. 
        statistic: string
            The way to combine vals to generate the histogram. Input to scipy `binned_statistic_dd`. 
            Note that cplt histograms generate densities, so whatever for statistic is computed, the number plotted is statistic/volume_of_bin
        **kwargs:
            Other arguments for voxelize (and the polygon collection). facecolor and edgecolor are overriden.

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
    '''
    return _histogram(pts,vals,X,Y,Z,dX,dY,dZ,statistic = statistic)

@plotting.volume_plot
@uniform_generator
def uniform_nearest_interpolator(pts,vals,X,Y,Z,dX,dY,dZ):
    '''A Method for Plotting a 3d Function in the cubic style. Interpolates to the centers of each cube using nearest neighbor interpolation.
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
        edgecolor_function: Function
            A function that takes in the facecolors and returns the edgecolors array. Can be useful if you want to shade your edge colors differently from your face colors (or change alpha)
        bins: int or list of ints of shape (3,)
            The number of bins to use. If a list, sets [X bins, Y bins, Z bins]. If an int, generates that number on each axis. 
        **kwargs:
            Other arguments for voxelize (and the polygon collection). facecolor and edgecolor are overriden.

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
    '''
    Xc,Yc,Zc = X+dX/2,Y+dY/2,Z+dZ/2
    to_plot = np.array([Xc.flatten(),Yc.flatten(),Zc.flatten()]).T
    
    t = time.time()
    print("interpolating...")
    I = interpolate.NearestNDInterpolator(pts, vals)
    print("done interpolating",time.time()-t,'seconds')
    
    t = time.time()
    print("evaluating...")
    results = I(to_plot)
    print("done evaluating",time.time()-t,'seconds')

    results = results.reshape(X.shape)
    results = np.swapaxes(results,0,1)
    return results

@plotting.volume_plot
@uniform_generator
def uniform_linear_interpolator(pts,vals,X,Y,Z,dX,dY,dZ):
    '''a 3d Function in the cubic style. Interpolates to the centers of each cube using linear interpolation. Note that this is much slower than nearest neighbor interpolation...
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
        edgecolor_function: Function
            A function that takes in the facecolors and returns the edgecolors array. Can be useful if you want to shade your edge colors differently from your face colors (or change alpha)
        bins: int or list of ints of shape (3,)
            The number of bins to use. If a list, sets [X bins, Y bins, Z bins]. If an int, generates that number on each axis. 
        **kwargs:
            Other arguments for voxelize (and the polygon collection). facecolor and edgecolor are overriden.

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
    '''
    Xc,Yc,Zc = X+dX/2,Y+dY/2,Z+dZ/2
    to_plot = np.array([Xc.flatten(),Yc.flatten(),Zc.flatten()]).T
    
    t = time.time()
    print("interpolating...")
    I = interpolate.LinearNDInterpolator(pts, vals)
    print("done interpolating",time.time()-t,'seconds')
    
    t = time.time()
    print("evaluating...")
    results = I(to_plot)
    print("done evaluating",time.time()-t,'seconds')

    results = results.reshape(X.shape)
    results = np.swapaxes(results,0,1)
    return results

@plotting.volume_plot
@tree_generator
def tree_histogram(pts,vals,X,Y,Z,dX,dY,dZ,statistic = 'sum'):
    '''Plots a histogram and assigns bin sizes using a tree method. Effectively, this attempts to place smaller bins in regions of the space that require a higher resolution. This attempts to plot fewer overall voxels on the matplotlib plot which is faster. Note that this is an **experimental** feature.
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
        edgecolor_function: Function
            A function that takes in the facecolors and returns the edgecolors array. Can be useful if you want to shade your edge colors differently from your face colors (or change alpha)
        min_resolution: int
            The minimum resolution to use. Note this is log2(bins) at minimal bin size
        max_resolution: int
            The maximum resolution to use. Note that this is log2(bins) at maximal bin size
        dist: list of floats, String, or None
            The distribution of the number of each bin size. Must sum to 1. For example [0.5,0.5] would cause the code to make as close to half the bins as possible level 1 bins and half the bins level 2 bins. If string will generate according to the method related to the string. Currently, there are two methods:
                - "equivolume" or None: computes a distribution assuming equal volumes.
                - "sigmoid": computes a distribution using a sigmoid distribution which is then adjusted by volume.
        bins: int or float
            int:  The total number of bins. This will not be the exact number of bins plotted, but will be the ballpark number which is aimed for while the program allocates bin sizes using `dist`.
            float: The percentage of volume taken up by the smallest bin size. Should be very small, especially for small bins, to prevent rendering issues. 
        focus: string
            How the code determines where to make smaller bins. The current two foci are 'slope' which extracts a percent with the largest slopes and 'magnitude' which extracts a percent with the largest magnitude.
        **kwargs:
            Other arguments for voxelize (and the polygon collection). facecolor and edgecolor are overriden.

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
    '''
    return _histogram(pts,vals,X,Y,Z,dX,dY,dZ,statistic = statistic)/(dX*dY*dZ)

def parallel_animate(fig,func,frames,result_name,
                     out = 'gif',fps = 50,parallel = True,Animation_Generation_Folder = None,
                     delete = True,merge = True,**kwargs):
    '''A function that generates an animation based on an input function. Can be run in parallel for complex plots.
        Parameters
        ----------
        fig : Figure
            The figure which contains the axes to animate
        func : Function
            A function that takes in a frame from frames and performs any modifications necessary to the axes.
        frames : list
            A list of inputs to func
        result_name: str
            The filename that the code saves the animation to
        out: str
            The filetype of the animation. Default is `.gif`
        fps: int
            The frames per second of the animation. 
        parallel: bool
            Whether or not to parallelize saving the frames of the animation.    
        Animation_Generation_Folder: str
            The folder to generate the animation frames without. If set to None will create folder with random integer name in the working directory.
        delete: bool
            Whether or not to delete the Animation_Generation_Folder and all of the frames within upon completion of the animation
        merge: bool
            Whether or not to merge the frames into an actual animation. 
        **kwargs:
            Arguments for the `savefig` operation on `fig`

        Returns
        -------
        None
    '''
    
    
    if not merge:
        delete = False
    
    N_cpu = multiprocess.cpu_count()
    if parallel:
        print('setup for parallel run')
        print(f"Using {N_cpu} cpus")
        N_per = int(len(frames)/N_cpu)+1
        frames_todo = [list(range(i,i+N_per)) for i in range(0,len(frames)-N_per,N_per)] 
        frames_todo += [list(range(frames_todo[-1][-1]+1,len(frames)))]
        # print(frames_todo)
    else:
        print('setup for linear run')
        frames_todo = None

    print('preparing file structure...')
    random_int = str(np.random.randint(0,10000)) + str(np.random.randint(0,10000))
    if Animation_Generation_Folder is None:
        if interactive:
            raise Exception("Requires a file path in interactive mode!")
        Animation_Generation_Folder = f'{os.path.dirname(os.path.realpath(sys.argv[0]))}/{random_int}'
    
    try:
        os.mkdir(Animation_Generation_Folder)
    except:
        print(f'{Animation_Generation_Folder} already exists, adding temp inside...')
        Animation_Generation_Folder = Animation_Generation_Folder + f"/{random_int}"
        os.mkdir(Animation_Generation_Folder)
    print('completed filestructure preparation')
    
    print(f'saving frames into {Animation_Generation_Folder}')
    print(f'{delete} - this folder and its contents will be deleted at the end of the run')
    if not delete:
        print(f'{merge} - this code will merge the frames into a single animation')
    if merge:
        print(f'the code will save into a(n) {out} file')

    print("preparing saving function...")
    def save(func,fig,frames,**kwargs):
        def wrapper(i):
            frame = frames[i]
            func(frame)
            fig.savefig(f'{Animation_Generation_Folder}/{i}.png',**kwargs)
        return wrapper
    
    func = save(func,fig,frames)
    
    
    def save_on_array(i):
        f = frames_todo[i]
        if interactive:
            print(' ', end='', flush=True) #issue with this turning off when verbose == False. Also need to fix the tqdm turning on upon verbose == False
        for i in tqdm.tqdm(f,desc = f'Animation {i}',total = len(frames_todo[i]),leave=True,position=i):
            func(i)


    print('saving functions prepared')

    print(f"generating {len(frames)} animation frames...")

    if parallel:
        with multiprocess.Pool(N_cpu) as p:
            p.map(save_on_array,list(range(len(frames_todo))))
                
    else:
        [func(i) for i in tqdm.trange(len(frames),desc = "Generating Animation",leave = True)]

    print('finished generating animation')
    
    if merge:
        print("combining animation...")
        
        lines = os.listdir(Animation_Generation_Folder)
        # lines = ['file '+f'"{Animation_Generation_Folder}/{i}"'+'\n' for i in lines]
        lines.sort(key = lambda u:int(u.split('.')[0]))
        # lines = [line + '.png' for line in lines]
        lines = [f'{Animation_Generation_Folder}/{line}' for line in lines]
        print(f"{len(lines)} frames found:")
    

        frames_png = []
        
        

        for line in tqdm.tqdm(lines,desc = "Loading Frames"):
            image = imageio.v2.imread(line)
            frames_png.append(image)
        
        print(f"saving at {fps} fps...")
        if out == 'gif':
            imageio.mimsave(f'{result_name}.gif',frames_png,fps = fps,loop = 0)
        else:
            imageio.mimsave(f'{result_name}.{out}',frames_png,fps = fps)
        
    if delete:
        print(f'removing animation generation folder: {Animation_Generation_Folder}')
        subprocess.call(['rm','-r', Animation_Generation_Folder])
    

@helper.verbose_print
def spin_3d_plot(fig,axs,result_name,times = 1,step = 1,parallel=True,
                 axis = 'z',fps = 50,Animation_Generation_Folder = None,
                 delete = True,merge = True,**kwargs):
    '''A function that generates an animation rotating the input 3d plot. Can run in parallel to speed up animation for complex plots.
        Parameters
        ----------
        fig : Figure
            The figure which contains the axes to rotate
        axs : Axis or list
            The axis on which to plot the points. If given a list it will animate all of the axes in the list (useful for animating a multi-cell plot). Make sure that all the axes are in the 3d projection!
        result_name : str
            The name of the file to save the animation to
        times: float
            The number of times to rotate. i.e. the function will rotate the plot int(360 * times) degrees.
        step: int
            The step size of the animation. i.e. the function will rotate the plot in steps of `step` degrees
        parallel: bool
            Whether or not to parallelize saving the frames of the animation.
        axis: str (or 3d vector)
            The axis to rotate around. Currently, only 'x', 'y', and 'z' are implemented. TODO implement 3d vector to rotate about arbitrary axis. 
        fps: int
            The frames per second of the animation. 
        Animation_Generation_Folder: str
            The folder to generate the animation frames without. If set to None will create folder with random integer name in the working directory.
        delete: bool
            Whether or not to delete the Animation_Generation_Folder and all of the frames within upon completion of the animation
        merge: bool
            Whether or not to merge the frames into an actual animation. 
        **kwargs:
            Arguments for the general parallel animation generation function.
                out: str
                    The filetype of the animation. Default is `.gif`
                **kwargs
                    What to input into the `fig` `savefig` method when saving frames (e.g. dpi)
        Returns
        -------
        None
    '''

    if not merge:
        delete = False
    
    if type(axs) != list:
        axs = [axs]

    # define the animation function
    def update(frame):
        for axi in axs:
            if axis == 'z':
                axi.view_init(elev=30, azim=frame) # Rotate around the z-axis
            elif axis == 'x':
                axi.view_init(elev=frame, azim=0) # Rotate around the x-axis
            elif axis == 'y':
                axi.view_init(elev = frame,azim = 90) # Rotate around the y-axis
            else:
                raise NotImplementedError(f"Axis {axis} not implemented")
    
    
    frames = np.arange(0,int(360 * times),step)
    
    print('animating...')
    parallel_animate(fig = fig,
                     func = update,
                     frames = frames,
                     result_name=result_name,
                     fps = fps,
                     parallel=parallel,
                     Animation_Generation_Folder=Animation_Generation_Folder,
                     delete = delete,merge = merge,**kwargs)