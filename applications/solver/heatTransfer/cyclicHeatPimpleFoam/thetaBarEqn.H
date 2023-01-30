 // Solve the Reduced Temperature Equation

sigma = ((alphaEff*pow(lambda,2)) + (lambda*(flowDir & U)))*thetaBar - (2*alphaEff*lambda*(flowDir & fvc::grad(thetaBar)));

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
