// Solve the Reduced Temperature Equation

fvScalarMatrix thetaBarEqn
(
fvm::ddt(thetaBar)
+ fvm::div(phi, thetaBar)
==
fvm::laplacian(alphaEff, thetaBar)
+ ((alphaEff*pow(lambda,2)) + (lambda*(flowDir & U)))*thetaBar - (2*alphaEff*lambda*(flowDir & fvc::grad(thetaBar)))
// + fvOptions(T)
);

thetaBarEqn.relax();

// fvOptions.constrain(Teqn);

thetaBarEqn.solve();
