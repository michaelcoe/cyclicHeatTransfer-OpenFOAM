// Calculate decay constant lambda

volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

label inletID = mesh.boundaryMesh().findPatchID("inlet");
label wallsID = mesh.boundaryMesh().findPatchID("walls");

if(runTime.value() >= delay.value())
{
    if(int(runTime.value()) % int(interval.value()) == 0)
    {
        // Calculate decay constant lambda
        dimensionedScalar lambdaNum
        (
            "lambdaNum",
            dimensionSet(0, 3, -1, 1, 0, 0, 0),
            gSum(alphaEff.boundaryField()[wallsID]*-theta.boundaryField()[wallsID].snGrad()*mesh.boundary()[wallsID].magSf())
        );

        dimensionedScalar lambdaDen
        (
            "lambdaDen",
            dimensionSet(0, 3, -1, 1, 0, 0, 0),
            gSum(thetaBar.boundaryField()[inletID]*(flowDir & U.boundaryField()[inletID])*mesh.boundary()[inletID].magSf()) + 
            gSum(alphaEff.boundaryField()[inletID]*(flowDir & fvc::grad(thetaBar)().boundaryField()[inletID])*mesh.boundary()[inletID].magSf())
        );
        
        Info << "Ux*Theta: " << gSum((thetaBar.boundaryField()[inletID]*(flowDir & U.boundaryField()[inletID])*mesh.boundary()[inletID].magSf())) << endl;
        Info << "alpha*dTheta/dx: " << gSum(alphaEff.boundaryField()[inletID]*(flowDir & fvc::grad(thetaBar)().boundaryField()[inletID])*mesh.boundary()[inletID].magSf()) << endl;

        Info << "lambda_num: " << lambdaNum.value() << endl;
        Info << "lambda_den: " << lambdaDen.value() << endl;
        Info << "lambda_num/den: " << lambdaNum.value()/lambdaDen.value() << endl;

        lambda = -(1/Lc)*Foam::log(1 - (lambdaNum/lambdaDen));
    }
}
else
{
    lambda = lambda0;
}

Info << "Lambda: " << lambda.value() << endl;
