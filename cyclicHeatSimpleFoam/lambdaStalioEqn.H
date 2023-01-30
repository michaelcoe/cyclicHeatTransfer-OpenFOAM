// Calculate decay constant lambda

volScalarField alphaEff("alphaEff",turbulence->nu()/Pr + turbulence->nut()/Prt);

label inletID = mesh.boundaryMesh().findPatchID("inlet");
label wallsID = mesh.boundaryMesh().findPatchID("walls");

if(runTime.value() >= 20)
{
    // Calculate decay constant lambda
    dimensionedScalar b0
    (
        "b0",
        dimensionSet(0, 3, -1, 1, 0, 0, 0),
        gSum(alphaEff.boundaryField()[wallsID]*theta.boundaryField()[wallsID].snGrad()*mesh.boundary()[wallsID].magSf())
    );

    dimensionedScalar b1
    (
        "b1",
        dimensionSet(0, 4, -1, 1, 0, 0, 0),
        0.02*gSum(0.5*(flowDir & U.boundaryField()[inletID])*)
    );

    dimensionedScalar b2
    (
        "b2",
        dimensionSet(0, 1, -1, 1, 0, 0, 0),
        gSum(alphaEff*thetaBar)
    );
    
    Info << "Gradient 1: " << sum((theta.boundaryField()[inletID]*(flowDir & U.boundaryField()[inletID])*mesh.boundary()[inletID].magSf())) << endl;
    Info << "Gradient 2: " << sum(alphaEff.boundaryField()[inletID]*(flowDir & fvc::grad(theta)().boundaryField()[inletID])*mesh.boundary()[inletID].magSf()) << endl;

    Info << "b0: " << b0.value() << endl;
    Info << "b1: " << b1.value() << endl;
    Info << "b2: " << b2.value() << endl; 
    
    lambda = (-b1 + Foam::sqrt(pow(b1,2)-(4*b2*b0)))/(2*b2);
}

Info << "Lambda: " << lambda.value() << endl;
