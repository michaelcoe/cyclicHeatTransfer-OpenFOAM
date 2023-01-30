Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// initialize the source term
// volScalarField sigma
// (
//     IOobject
//     (
//         "sigma",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("sigma",dimensionSet(0, 0, -1, 1, 0, 0, 0), 0.0)
// );

Info<< "Reading thermalProperties\n" << endl;
IOdictionary thermalProperties
(
    IOobject
    (
        "thermalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar Pr 
(
    "Pr",
    dimless,
    thermalProperties
);

dimensionedScalar Prt 
(
    "Prt",
    dimless,
    thermalProperties
);

dimensionedScalar Lc 
(
    "Lc",
    dimLength,
    thermalProperties
);

// dimensionedScalar inletH 
// (
//     "inletH",
//     dimLength,
//     thermalProperties
// );

dimensionedScalar Twall 
(
    "Twall",
    dimTemperature,
    thermalProperties
);

dimensionedScalar lambda0 
(
    "lambda0",
    dimensionSet(0, -1, 0, 0, 0, 0, 0),
    thermalProperties
);

dimensionedScalar lambda 
(
    "lambda",
    dimensionSet(0, -1, 0, 0, 0, 0, 0),
    lambda0.value()
);

dimensionedVector Ubar 
(
    "Ubar",
    dimVelocity,
    thermalProperties
);

volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    T-Twall
);

// initialize the source term
volScalarField thetaBar
(
    IOobject
    (
        "thetaBar",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    T
);

// Allocate and define variables that don't change
vector flowDir = (Ubar.value()/mag(Ubar.value()));

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
