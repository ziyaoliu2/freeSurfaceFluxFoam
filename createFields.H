Info << "Reading field C" << endl;
volScalarField C
(
  IOobject
  (
    "C",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
  ),
  mesh
);


Info << "Creating field volFluxC" << endl;
volVectorField volFluxC
( 
  IOobject
  ( 
    "volFluxC",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedVector("volFluxC", dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0))
);

Info << "Creating field fluxC" << endl;
surfaceScalarField fluxC
(
  IOobject
  (
    "fluxC",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
  ),
  mesh,
  dimensionedScalar("fluxC", dimensionSet(0,1,-1,0,0,0,0), scalar(0))
);

Info<< "Reading field p (kinematic)\n" << endl;
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

Info << "Reading transportProperties" << nl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar D("D", dimensionSet(0,2,-1,0,0,0,0), transportProperties);

dimensionedScalar nu("nu", dimensionSet(0,2,-1,0,0,0,0), transportProperties);

#include "createPhi.H"

Info << "Reading convDiff subdirectory" << nl;

dictionary convDiff = mesh.solutionDict().subDict("convDiff");

double tolerance = 0;       // convergence criteria for convection-diffuion
if( !convDiff.readIfPresent<double>("tolerance", tolerance) ){
  SeriousErrorIn("main")
      << "No tolerance parameter in convDiff subdictionary"
      << exit(FatalError);
}
int maxIter = 0;             // maximum number of iterations
if( !convDiff.readIfPresent<int>("maxIter", maxIter) ){
  SeriousErrorIn("main")
      << "No maxIter parameter in convDiff subdictionary"
      << exit(FatalError);
}


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "readGravitationalAcceleration.H"

Info<< "Creating field zeta\n" << endl;
volVectorField zeta
(
    IOobject
    (
        "zeta",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimLength, Zero)
);

Info<< "Creating field p_gh\n" << endl;
volScalarField p_gh
(
    IOobject
    (
        "p_gh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Force p_gh to be consistent with p
// Height is made relative to field 'refLevel'
p_gh =  p - (g & mesh.C());


label p_ghRefCell = 0;
scalar p_ghRefValue = 0.0;
setRefCell(p_gh, pimple.dict(), p_ghRefCell, p_ghRefValue);
mesh.setFluxRequired(p_gh.name());

#include "createMRF.H"
#include "createFvOptions.H"
