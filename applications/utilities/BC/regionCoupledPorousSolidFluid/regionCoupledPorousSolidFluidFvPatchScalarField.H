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

Class
    Foam::regionCoupledPorousSolidFluidFvPatchScalarField

Description
    Mixed boundary condition for potential in a solid, to be used for secondary 
    distribution problems.

    See paper description: https://doi.org/10.1016/j.electacta.2018.09.121


Usage
    \table
        Property     | Description             | Required    | Default value
        side         | name of the side    | yes | solid or fluid
        nbrField     | name of the field    | no | fi
        kappaMethod  | do not touch | inherited |
        kappa        | electrical conductivity | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {

	    type	    regionCoupledPorousSolidFluid;
	    nbrField        fi;
            kappaMethod     lookup;
            kappa           ks or kf;
            value           $internalField;
        
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.

See also
    Foam::compressible::turbulentTemperatureCoupledBaffleMixed

SourceFiles
   regionCoupledPorousSolidFluidFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledPorousSolidFluidFvPatchScalarField_H
#define regionCoupledPorousSolidFluidFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class regionCoupledPorousSolidFluidFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledPorousSolidFluidFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase//potentialCoupledBase
{
    // Private data

        //- Name of field on the neighbour region
        const word neighbourFieldName_;
        //- Name of the region side
        //const word sideName_; 
        //- Kinetic parameters
        //scalarList exchangeCurrentDensity_;
        //scalarList beta_;
        //scalarList equilibriumPotential_;

public:

    //- Runtime type information
    TypeName("regionCoupledPorousSolidFluid");

    // Constructors

        //- Construct from patch and internal field
        regionCoupledPorousSolidFluidFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupledPorousSolidFluidFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  regionCoupled onto a
        //  new patch
        regionCoupledPorousSolidFluidFvPatchScalarField
        (
            const regionCoupledPorousSolidFluidFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new regionCoupledPorousSolidFluidFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        regionCoupledPorousSolidFluidFvPatchScalarField
        (
            const regionCoupledPorousSolidFluidFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new regionCoupledPorousSolidFluidFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
