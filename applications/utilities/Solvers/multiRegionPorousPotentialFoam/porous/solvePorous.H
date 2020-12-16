#include "solid_potential.H"
#include "fluid_potential.H"
 
#include "CEqn.H"

{
// Total current
    scalar sumField = 0.0;
    scalar sumField_bipolar = 0.0;
    {
        sumField = gSum((mesh.V()*Jvol.internalField())());
	sumField_bipolar = gSum((mesh.V()*mag(Jvol.internalField()))());
    }

    if (mag(sumField) > mag(sumField_bipolar/2.0))
    {
	Info<< "Total current at the monopolar electrode = " << sumField << " A "<< endl;
    }
    else
    {
	Info<< "Total current at the bipolar electrode = " << sumField_bipolar/2.0 << " A "<< endl;
    }
}
