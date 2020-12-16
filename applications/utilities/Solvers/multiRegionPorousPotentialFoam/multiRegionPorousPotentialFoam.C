/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is NOT part of OpenFOAM.

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
    multiRegionPorousPotentialFoam

Description
    Steady-state solver for porous electrodes. It solves the potential field
    with conjugate potential transfer between solid, porous and fluid regions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"

#include "fvOptions.H"
#include "coordinateSystem.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"

    // Residual control
    bool allRegionsConverged = false;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        forAll(fluidRegions, i)
        {
            Info<< "\nSolving for fluid region "
                << fluidRegions[i].name() << endl;
            #include "setRegionFluidFields.H"
            #include "readFluidMultiRegionSIMPLEControls.H"
	    #include "initFluidConvergenceCheck.H"
            #include "solveFluid.H"
	    #include "convergenceFluidCheck.H"
        }

	forAll(porousRegions, i)
        {
            Info<< "\nSolving for porous region "
                << porousRegions[i].name() << endl;
            #include "setRegionPorousFields.H"
            #include "readPorousMultiRegionSIMPLEControls.H"
	    #include "initPorousConvergenceCheck.H"
            #include "solvePorous.H"
	    #include "convergencePorousCheck.H"
        }

        forAll(solidRegions, i)
        {
            Info<< "\nSolving for solid region "
                << solidRegions[i].name() << endl;
            #include "setRegionSolidFields.H"
            #include "readSolidMultiRegionSIMPLEControls.H"
	    #include "initSolidConvergenceCheck.H"
            #include "solveSolid.H"
	    #include "convergenceSolidCheck.H"
        }

	#include "checkResidualControls.H"

        runTime.write();

        Info<< " ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
