// Solve the Reduced Temperature Equation
// label outletID = mesh.boundaryMesh().findPatchID("outlet");

forAll(theta, celli)
{
    theta[celli] = Foam::exp(-lambda.value()*(flowDir & mesh.C()[celli]))*thetaBar[celli];
}

theta = theta/gAverage(theta.boundaryField()[inletID]);

T = theta + Twall;

rhok = 1.0 - beta*(T - TRef);

dimensionedScalar TbulkNum
(
    "TbulkNum",
    dimensionSet(0, 2, -1, 1, 0, 0, 0),
    gSum(T.boundaryField()[inletID]*(flowDir & U.boundaryField()[inletID])*mesh.boundary()[inletID].magSf())
);

dimensionedScalar TbulkDen
(
    "TbulkDen",
    dimensionSet(0, 2, -1, 0, 0, 0, 0),
    gSum((flowDir & U.boundaryField()[inletID])*mesh.boundary()[inletID].magSf())
);

dimensionedScalar Tbulk
(
    TbulkNum/TbulkDen
);

Info << "Tbulk: " << Tbulk.value() << endl;

dimensionedScalar dTdn
(
    "dTdn",
    dimensionSet(0, -1, 0, 1, 0, 0, 0),
    gSum(T.boundaryField()[wallsID].snGrad()*mesh.boundary()[wallsID].magSf())/gSum(mesh.boundary()[wallsID].magSf())
);

Nu = dTdn*Lc/(Twall-Tbulk);

Info << "Nu: " << Nu.value() << endl;