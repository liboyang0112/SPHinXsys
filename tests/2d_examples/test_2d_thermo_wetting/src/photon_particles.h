/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
 * @file 	photon_particles.h
 * @brief 	This is the derived class of base particles.
 * @author	Boyang Li
 */

#ifndef PHOTON_PARTICLES_H
#define PHOTON_PARTICLES_H

#include "base_particles.h"
#include "base_particles.hpp"
#include "SPHconfig.h"
#include "particle_generator_lattice.h"
namespace SPH
{

	//----------------------------------------------------------------------
	//		preclaimed classes
	//----------------------------------------------------------------------
	class Photon;
	/**
	 * @class PhotonParticles
	 * @brief A group of particles with photon body particle data.
	 */
	class PhotonParticles : public BaseParticles
	{
	public:
		PhotonParticles(SPHBody &sph_body,
					   SharedPtr<Photon> shared_photon_ptr,
					   SharedPtr<ParticleGenerator> particle_generator_ptr = makeShared<ParticleGeneratorLattice>());
		virtual ~PhotonParticles(){};

		StdLargeVec<Real> intensity_n_;
		StdLargeVec<void*> inside_body_n_;  // The pointer to the material particles containing the photon.
		StdLargeVec<Real> rho_gradiant_prev_n_; // used for reflection determination: reflect when rho_gradiant reaches max.
		StdLargeVec<Real> rho_gradiant_prev_max_n_; // used for reflection determination: reflect when rho_gradiant reaches max.

		virtual PhotonParticles *ThisObjectPtr() override { return this; };
	};
}
#endif