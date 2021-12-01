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

Class
    Foam::ParticleDepositionRecord

Description
    Creates particle erosion field, Q

SourceFiles
    ParticleDepositionRecord.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleDepositionRecord_H
#define ParticleDepositionRecord_H

#include "CloudFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ParticleDepositionRecord Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleDepositionRecord
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


        autoPtr<volScalarField> MassDeposition_;

        //- List of patch indices to post-process
        labelList patchIDs_;

        //- dmax
        scalar dmax_;

        //- dmin
        scalar dmin_;
		
        //- posImax
        scalar posxImax_;
        
        //- posImin
        scalar posxImin_;
        
        //- posImax
        scalar posyImax_;
        
        //- posImin
        scalar posyImin_;

        //- posImax
        scalar poszImax_;
        
        //- posImin
        scalar poszImin_;
        
        label group_;

protected:

    // Protected Member Functions

        //- Returns local patchI if patch is in patchIds_ list
        label applyToPatch(const label globalPatchI) const;

        //- Write post-processing info
        virtual void write();


public:

    //- Runtime type information
    TypeName("ParticleDepositionRecord");


    // Constructors

        //- Construct from dictionary
        ParticleDepositionRecord
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ParticleDepositionRecord(const ParticleDepositionRecord<CloudType>& pe);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType> > clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType> >
            (
                new ParticleDepositionRecord<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleDepositionRecord();


    // Member Functions

        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve();

            //- Post-patch hook
            virtual void postPatch
            (
                const parcelType& p,
                const polyPatch& pp,
                const scalar trackFraction,
                const tetIndices& tetIs,
                bool& keepParticle
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ParticleDepositionRecord.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
