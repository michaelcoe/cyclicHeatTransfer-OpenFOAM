// Solve the Reduced Temperature Equation

label outletID = mesh.boundaryMesh().findPatchID("outlet");

forAll(theta, celli)
{
    theta[celli] = Foam::exp(-lambda.value()*(flowDir & mesh.C()[celli]))*thetaBar[celli];
}

forAll(theta.boundaryField()[inletID], facei)
{
    theta.boundaryFieldRef()[inletID][facei] = Foam::exp(-lambda.value()*
    (flowDir & mesh.boundaryMesh()[inletID].faceCentres()[facei]))*thetaBar.boundaryField()[inletID][facei];
}

forAll(theta.boundaryField()[outletID], facei)
{
    theta.boundaryFieldRef()[outletID][facei] = Foam::exp(-lambda.value()*
    (flowDir & mesh.boundaryMesh()[outletID].faceCentres()[facei]))*thetaBar.boundaryField()[outletID][facei];
}

dimensionedScalar TbulkNum
(
    "TbulkNum",
    dimensionSet(0, 2, -1, 1, 0, 0, 0),
    gSum(T.boundaryField()[inletID]*U.boundaryField()[inletID].component(0)*mesh.magSf().boundaryField()[inletID])
);

dimensionedScalar TbulkDen
(
    "TbulkDen",
    dimensionSet(0, 2, -1, 0, 0, 0, 0),
    gSum(U.boundaryField()[inletID].component(0)*mesh.magSf().boundaryField()[inletID])
);

dimensionedScalar Tbulk
(
    TbulkNum/TbulkDen
);

T = theta + Twall;

T = T/Tbulk.value();
