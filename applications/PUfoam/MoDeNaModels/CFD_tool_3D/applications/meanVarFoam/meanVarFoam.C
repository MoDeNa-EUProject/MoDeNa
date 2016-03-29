/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM. If not, see <http://www.gnu.org/licenses/>.

Application
    meanVarFoam

Description
    For each time: calculate the mean and variance of bubble size distribution

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "twoPhaseMixtureThermo.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class BSDProperties
{
private:
    volScalarField& M0;
    volScalarField& M1;
    volScalarField& M2;
    const volScalarField& alpha1;
    const fvMesh& mesh;

public:
    // constructors  
    BSDProperties
    (
        volScalarField &M0_, volScalarField &M1_, 
        volScalarField &M2_, volScalarField &alpha1_,
        const fvMesh &mesh_
    )
    :
    M0(M0_),
    M1(M1_),
    M2(M2_),
    alpha1(alpha1_),
    mesh(mesh_)
    {}

    // member functions
    volScalarField meanBSD()
    {
        volScalarField alphaCutOff
        (
            IOobject
            (
                "alphaCutOff",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("alphaCutOff", dimless, 0.5)
        );
        volScalarField insideFoam(pos(alphaCutOff - alpha1));
        volScalarField mean_
        (
            IOobject
            (
                "mean_",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("mean_", dimless, 0.0)
        );
        forAll(mesh.C(), celli)
        {
            M0 = Foam::max(M0, dimensionedScalar("M0min", M0.dimensions(), SMALL));
            M1 = Foam::max(M1, dimensionedScalar("M1min", M1.dimensions(), SMALL));
            M2 = Foam::max(M2, dimensionedScalar("M2min", M2.dimensions(), SMALL));
            mean_[celli] = scalar(2.0)*Foam::log(M1[celli]/M2[celli]) - 
                           scalar(0.5)*Foam::log(M2[celli]/M0[celli]);
        }
        return (insideFoam*mean_);
    }
    
    volScalarField varianceBSD()
    {
        volScalarField alphaCutOff
        (
            IOobject
            (
                "alphaCutOff",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("alphaCutOff", dimless, 0.5)
        );
        volScalarField insideFoam(pos(alphaCutOff - alpha1));
        volScalarField variance_
        (
            IOobject
            (
                "variance_",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("variance_", dimless, 0.0)
        );
        forAll(mesh.C(), celli)
        {
            M0 = Foam::max(M0, dimensionedScalar("M0min", M0.dimensions(), SMALL));
            M1 = Foam::max(M1, dimensionedScalar("M1min", M1.dimensions(), SMALL));
            M2 = Foam::max(M2, dimensionedScalar("M2min", M2.dimensions(), SMALL));
            variance_[celli] = 
                (
                    scalar(2.0)*(
                    scalar(0.5)*Foam::log(M2[celli]/M0[celli]) 
                  - Foam::log(M1[celli]/M0[celli])     
                    )
                );
            if (variance_[celli] < 0.0)
            {
                variance_[celli] = 0.0;
            }
        }
        return (insideFoam*variance_);
    }
};

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
        IOobject M0header
        (
            "M0",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
        
        IOobject M1header
        (
            "M1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        IOobject M2header
        (
            "M2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
        
        // Check M0, M1 and M2 exist
        if 
        (
            M0header.headerOk()
         && M1header.headerOk()
         && M2header.headerOk()
        )
        {
            //mesh.readUpdate();

            Info<< "    Reading M0" << endl;
            volScalarField M0(M0header, mesh);

            Info<< "    Reading M1" << endl;
            volScalarField M1(M1header, mesh);

            Info<< "    Reading M2" << endl;
            volScalarField M2(M2header, mesh);

            Info<< "    Calculating BSD properties: " << endl;

            twoPhaseMixtureThermo mixture(mesh);
            volScalarField& alpha1(mixture.alpha1());
            
            BSDProperties mu(M0, M1, M2, alpha1, mesh);
            volScalarField BSDMean
            (
                IOobject
                    (
                        "BSDMean",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mu.meanBSD()
            );
            Info<< "BSDMean = " << BSDMean.internalField() << endl;
            BSDMean.write();

            BSDProperties sigmaTwo(M0, M1, M2, alpha1, mesh);
            volScalarField BSDVariance
            (
                IOobject
                    (
                        "BSDVariance",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    sigmaTwo.varianceBSD()
            );
            Info<< "BSDVariance = " << BSDVariance.internalField() << endl;
            BSDVariance.write();
        }
        else
        {
            Info<< "    Moments are not available." << endl;
        }

        Info<< "\nEnd\n" << endl;
}

