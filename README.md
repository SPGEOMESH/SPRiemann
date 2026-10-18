

# SPRiemann

![Development Status](https://img.shields.io/badge/status-experimental-orange)

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

> *Will add more documentation.*


## Acknowledgements

This project is supported by the National Science Foundation:

| Award Number                                                                  | Principal Investigator                 |
|:------------------------------------------------------------------------------|:---------------------------------------|
| **[NSF-2514037](https://www.nsf.gov/awardsearch/show-award/?AWD_ID=2514037)** | Dr. Vladimir Florinski (UAH)           |
| **[NSF-2514038](https://www.nsf.gov/awardsearch/show-award/?AWD_ID=2514038)** | Dr. Dinshaw Balsara (Univ. Notre Dame) |

*Disclaimer: Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.*


