

# SPRiemann

SPRiemann is a C++ implementation of the [Cube Sphere](https://www.github.com/SPGEOMESH/CubeSphere) Fortran library. 
This reconstruction of the Cube Sphere code moves communication and other features over to the C++ side, and extends support from single-level cube sphere meshes to a state-of-the-art spherical geodesic adaptive mesh refinement (AMR) mesh. 

# Steps to run

Navigate to a directory where you would like to install. 
Then:
```bash
pip install -U dgcalc
dgcalc get DGCALC SPGEOMESH
```
Next, compile the CubeSphere library:
```bash
dgcalc cd CubeSphere
dgcalc config
dgcalc make
```
Now you can start working with the CSRiemann project. 
For example:
```
dgcalc cd CSRiemann
dgcalc config
dgcalc make
# run a test (edit `main.cpp` to select a target)
dgcalc run
```

Some tests will run without first compiling CubeSphere.

Will add more documentation later. 


