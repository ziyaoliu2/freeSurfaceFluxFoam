/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    freeSurfaceFluxFoam

Group
    grpMultiphaseSolvers

Description
    Incompressible Navier-Stokes solver with inclusion of a wave height field
    to enable single-phase free-surface approximations

    Wave height field, zeta, used by pressure boundary conditions

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Incompressible Navier-Stokes solver with inclusion of a wave height"
        " field to enable single-phase free-surface approximations."
    );

    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
//        #include "CourantNo.H"
//        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }
    runTime.printExecutionTime(Info);
   
    Info << "Begin cycle: Time = " << runTime.timeName()
         << nl << endl;

    int  iter = 0;
    while ( true )
    {
        iter++;

        fvScalarMatrix CEqn
        (
            fvm::ddt(C) + fvm::div(phi, C) - fvm::laplacian(D, C)
        );

        CEqn.relax();
        double residual = solve( CEqn ).initialResidual();

        if( residual < tolerance )
        {
            Info << "Convection-diffusion: "
                 << "ExecutionTime = " << runTime.elapsedCpuTime() << "s"
                 << "ClockTime = " << runTime.elapsedClockTime() << "s"
                 << nl << "Converged in " << iter << " steps.  Residual = "
                 << residual << nl << endl;

            if(iter >= maxIter)
            {
                Info << nl << "dissolFoam Runtime WARNING:"
                     << "Convection-diffusion solver did not converge."<< nl
                     << "Maximum number of iterations"
                     << "  iter: "<< iter << endl;
            }
        break;
        }
        else
        {
            Info << " Step " << iter
                 << " residual: "<< residual << " > " << tolerance << endl;
        }
    }

        volFluxC == U*C - D*fvc::grad(C);

        fluxC == phi/mag(mesh.Sf())*fvc::interpolate(C) - D*fvc::snGrad(C);

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
