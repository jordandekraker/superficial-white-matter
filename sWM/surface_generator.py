# code from https://github.com/khanlab/hippunfold/blob/master/hippunfold/workflow/scripts/create_warps.py

# shifts a wm surface inward along a Laplace field

import copy
import nibabel as nib
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import sys

print('starting surface shift')

# here we will set some parameters
in_surf = sys.argv[1]
in_laplace = sys.argv[2]
out_surf_prefix = sys.argv[3]
def arg2float_list(arg):
    return list(map(float, arg.split(',')))
if len(sys.argv)>4:
    depth = arg2float_list(sys.argv[4])
else:
    depth = [1,2,3] # default depths

convergence_threshold = 1e-4
step_size = 0.1 # mm
max_iters = int(depth[-1]/step_size)


# load data
surf = nib.load(in_surf)
V = surf.get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0].data
F = surf.get_arrays_from_intent('NIFTI_INTENT_TRIANGLE')[0].data
laplace = nib.load(in_laplace)
lp = laplace.get_fdata()
print('loaded data and parameters')

# laplace to gradient
dx,dy,dz = np.gradient(lp)

# make interpolator of gradients
points = (range(lp.shape[0]), range(1,lp.shape[1]+1), range(1,lp.shape[2]+1))
interp_x = RegularGridInterpolator(points, dx)
interp_y = RegularGridInterpolator(points, dy)
interp_z = RegularGridInterpolator(points, dz)
print('gradient interpolator ready')


distance_travelled = np.zeros((len(V)))
n=0
for d in depth:
    # apply inverse affine to surface to get to matrix space
    print(laplace.affine)
    V[:,:] = V - laplace.affine[:3,3].T
    for xyz in range(3):
        V[:,xyz] = V[:,xyz]*(1/laplace.affine[xyz,xyz])
    for i in range(max_iters):
        Vnew = copy.deepcopy(V)
        pts = distance_travelled < d
        stepx = interp_x(V[pts,:])
        stepy = interp_y(V[pts,:])
        stepz = interp_z(V[pts,:])
        magnitude = np.sqrt(stepx**2 + stepy**2 + stepz**2)
        if len(magnitude)>0:
            for m in range(len(pts)):
                if magnitude[m]>0:
                    stepx[m] = stepx[m] * (step_size/magnitude[m])
                    stepy[m] = stepy[m] * (step_size/magnitude[m])
                    stepy[m] = stepz[m] * (step_size/magnitude[m])
        Vnew[pts,0] += stepx
        Vnew[pts,1] += stepy
        Vnew[pts,2] += stepz
        distance_travelled[pts] += step_size
        ssd = np.sum((V-Vnew)**2,axis=None)
        print(f'itaration {i}, convergence: {ssd}, still moving: {np.sum(pts)}')
        if ssd < convergence_threshold:
            break
        V[:,:] = Vnew[:,:]
    # return to world coords
    for xyz in range(3):
        V[:,xyz] = V[:,xyz]*(laplace.affine[xyz,xyz])
    V[:,:] = V + laplace.affine[:3,3].T

    nib.save(surf, out_surf_prefix + str(d) + 'mm.surf.gii')
