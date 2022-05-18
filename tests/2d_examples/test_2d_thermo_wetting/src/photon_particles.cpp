/**
 * @file photon_particles.cpp
 * @brief Definition of functions declared in photon_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 */
#include "photon_particles.h"
#include "photon.h"
namespace SPH
{
	//=============================================================================================//
	PhotonParticles::PhotonParticles(SPHBody &sph_body,
								   SharedPtr<Photon> shared_photon_ptr,
								   SharedPtr<ParticleGenerator> particle_generator_ptr)
		: BaseParticles(sph_body, shared_photon_ptr, particle_generator_ptr)
	{
		shared_photon_ptr->assignPhotonParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<indexScalar, Real>(intensity_n_, "Intensity");
		registerAVariable<indexScalar, Real>(rho_gradiant_prev_n_, "Rho_Gradiant_Prev_");
		//-----------------------------------------------------------------------------------------
		//		register sortable particle data before building up particle configuration
		//-----------------------------------------------------------------------------------------
		registerASortableVariable<indexVector, Vecd>("Position");
		registerASortableVariable<indexScalar, Real>("Intensity");
		inside_body_n_.resize(real_particles_bound_, 0);
		//sorting particle once
		//DynamicCast<RealBody>(this, body)->sortParticleWithCellLinkedList();
	}
	//=================================================================================================//
}
