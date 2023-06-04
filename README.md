This is the `OpenFoam`&reg; solver for transient flow and reactive transport with dissolution. 

This solver is designed for transport limit dissoluiton, assuming time scale in dissolution is much larger than time scale in fluid flow, no meshing motion involved here.

PIMPLE algorithm is used to deal with pressure velocity coupling in the Navier-Stokes flow. Mass transport is simulated using the resulted flux (at the reactive surface) from flow simualtion. 

This solver added a pusedo freesurface boundary condition, which can be traced to `potentialFreeSurfaceFoam`&reg; in OpenFOAM repositories. The boundary condition enables simulations for pesodu two-phase flow with much lower computational cost.

This version was developed with OpenFOAM-v1912.

Additional documentation about the libraries, solver and cases can be found in the .zip attachment to the each release.

If you use these codes for published work, please cite the following references:

Liu, Z., & Ladd, A. (2022). Onset of instabilities in rotating flows by direct numerical simulation. Journal of Fluid Mechanics, 945, A31. doi:10.1017/jfm.2022.566
