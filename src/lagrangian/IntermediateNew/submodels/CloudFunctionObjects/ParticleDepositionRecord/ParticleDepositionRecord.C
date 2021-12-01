/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "ParticleDepositionRecord.H"

// * * * * * * * * * * * * * Protectd Member Functions * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::ParticleDepositionRecord<CloudType>::applyToPatch
(
    const label globalPatchI
) const
{
    forAll(patchIDs_, i)
    {
        if (patchIDs_[i] == globalPatchI)
        {
            return i;
        }
    }
    return -1;
}


template<class CloudType>
void Foam::ParticleDepositionRecord<CloudType>::write()
{
    if (MassDeposition_.valid())
    {
        MassDeposition_->write();
    }
    else
    {
        FatalErrorIn("void Foam::ParticleDepositionRecord<CloudType>::write()")
            << "MassDeposition not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleDepositionRecord<CloudType>::ParticleDepositionRecord
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    MassDeposition_(NULL),
    patchIDs_(),
    dmax_(this->coeffDict().template lookupOrDefault<scalar>("dmax", 1000000.0)),
    dmin_(this->coeffDict().template lookupOrDefault<scalar>("dmin", 0.0)),
    posxImax_(this->coeffDict().template lookupOrDefault<scalar>("posxImax", 1000000.0)),
    posxImin_(this->coeffDict().template lookupOrDefault<scalar>("posxImin", -1000000.0)),
    posyImax_(this->coeffDict().template lookupOrDefault<scalar>("posyImax", 1000000.0)),
    posyImin_(this->coeffDict().template lookupOrDefault<scalar>("posyImin", -1000000.0)),
    poszImax_(this->coeffDict().template lookupOrDefault<scalar>("poszImax", 1000000.0)),
    poszImin_(this->coeffDict().template lookupOrDefault<scalar>("poszImin", -1000000.0)),
    group_(this->coeffDict().template lookupOrDefault<label>("group", 1))
{
    const wordList allPatchNames = owner.mesh().boundaryMesh().names();
    wordList patchName(this->coeffDict().lookup("patches"));

    labelHashSet uniquePatchIDs;
    forAllReverse(patchName, i)
    {
        labelList patchIDs = findStrings(patchName[i], allPatchNames);

        if (patchIDs.empty())
        {
            WarningIn
            (
                "Foam::ParticleDepositionRecord<CloudType>::ParticleDepositionRecord"
                "("
                    "const dictionary&, "
                    "CloudType& "
                ")"
            )   << "Cannot find any patch names matching " << patchName[i]
                << endl;
        }

        uniquePatchIDs.insert(patchIDs);
    }

    patchIDs_ = uniquePatchIDs.toc();

    // trigger ther creation of the Q field
    preEvolve();
}


template<class CloudType>
Foam::ParticleDepositionRecord<CloudType>::ParticleDepositionRecord
(
    const ParticleDepositionRecord<CloudType>& pdr
)
:
    CloudFunctionObject<CloudType>(pdr),
    MassDeposition_(NULL),
    patchIDs_(pdr.patchIDs_),
    dmax_(pdr.dmax_),
    dmin_(pdr.dmin_),
    posxImax_(pdr.posxImax_),
    posxImin_(pdr.posxImin_),
    posyImax_(pdr.posyImax_),
    posyImin_(pdr.posyImin_),
    poszImax_(pdr.poszImax_),
    poszImin_(pdr.poszImin_),
    group_(pdr.group_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleDepositionRecord<CloudType>::~ParticleDepositionRecord()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleDepositionRecord<CloudType>::preEvolve()
{
    if (MassDeposition_.valid())
    {
        MassDeposition_->internalField() = 0.0;
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        MassDeposition_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "MassDepositionPerArea" +"_G_" + name(group_),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimVolume, 0.0)
            )
        );
    }
}


template<class CloudType>
void Foam::ParticleDepositionRecord<CloudType>::postPatch
(
    const parcelType& p,
    const polyPatch& pp,
    const scalar trackFraction,
    const tetIndices& tetIs,
    bool&
)
{
	if(((p.d()) < dmax_) && ((p.d()) >= dmin_) &&
	 ((p.PosI()).x() < posxImax_) && ((p.PosI()).x() >= posxImin_) &&
	 ((p.PosI()).y() < posyImax_) && ((p.PosI()).y() >= posyImin_) &&
	 ((p.PosI()).z() < poszImax_) && ((p.PosI()).z() >= poszImin_))
	 {
        const label patchI = pp.index();

        const label localPatchI = applyToPatch(patchI);

        if (localPatchI != -1)
        {
            vector nw;
            vector Up;

            // patch-normal direction
            this->owner().patchData(p, pp, trackFraction, tetIs, nw, Up);

            // particle velocity reletive to patch
            const vector& U = p.U() - Up;

            // quick reject if particle travelling away from the patch
            if ((nw & U) < 0)
            {
                return;
            }
        
            const label patchFaceI = pp.whichFace(p.face());
        
            const fvMesh& mesh = this->owner().mesh();
            const surfaceScalarField& magSf = mesh.magSf(); 
            scalar area = magSf.boundaryField()[patchI][patchFaceI];
        
            scalar& MassDeposition = MassDeposition_->boundaryField()[patchI][patchFaceI];
        
            MassDeposition += p.nParticle()*p.mass()/area;
        }
    }
}


// ************************************************************************* //
