# MarchingCubes with Color

This is a python wrapper to get colored mesh using marching cubes.

prerequisites:

- python 3.8
- cmake version 3.4


set up:

- git clone --recurse-submodules https://github.com/Sangluisme/ColoredMarchingCubes.git



build python wrapper

```
python setup.py install
```

after build the c++ python wrapper, you can use it as a python package.

here is a small example:

```
import colored_marching_cubes as mc

grid_dim = np.array([128,128,128], dtype=np.int)
grid_size = (0.02 * grid_dim).astype(np.float32)
origin = np.array([1, 1, 1], dtype=np.float32)

# build the marching cube object
mesh_mc = mc.MarchingCubes(grid_dim, grid_size, origin)

# volumn is sdf value of the voxels
# red, green, blue is the color of the voxels
# volumn, and red, green, blue should be numpy array with size of grid_dim
# volumn should be float, red, green, blue is int8
vertices, faces = mc.computeIsoSurface(mesh_mc, volumn, red, green, blue)

# save mesh return a bool to indicate if saved successfully
if_saved = mc.savePly(mesh_mc, filename) 
```

trouble shooting:

- wrong mesh: c++ code take pointer as input and compute voxels in (i,j,k) order, so try  
```
volumn.permute(2,1,0)

# permute color as well 
```
if your volumn is computed in python, especially using neural implicit functions.

(optional) compile c++ code

```
mkdir build
cd build
cmake ..
make -j

```