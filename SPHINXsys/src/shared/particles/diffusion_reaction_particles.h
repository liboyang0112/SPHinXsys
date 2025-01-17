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
* @file 	diffusion_reaction_particles.h
* @brief 	This is the derived class of diffusion reaction particles.
* @author	Xiangyu Huand Chi Zhang
*/

#ifndef DIFFUSION_REACTION_PARTICLES_H
#define DIFFUSION_REACTION_PARTICLES_H

#include "base_particles.h"
#include "base_body.h"
#include "base_material.h"
#include "diffusion_reaction.h"
#include "xml_engine.h"
#include "solid_particles.h"

namespace SPH
{

	/**
	 * @class DiffusionReactionParticles
	 * @brief A group of particles with diffusion or/and reactions particle data.
	 */
	template <class BaseParticlesType = BaseParticles, class BaseMaterialType = BaseMaterial>
	class DiffusionReactionParticles : public BaseParticlesType
	{
	protected:
		size_t number_of_species_;			 /**< Total number of diffusion and reaction species . */
		size_t number_of_diffusion_species_; /**< Total number of diffusion species . */
		std::map<std::string, size_t> species_indexes_map_;
		StdVec<size_t> species_diffusion_source_index_;

	public:
		StdVec<StdLargeVec<Real>> species_n_;	 /**< array of diffusion/reaction scalars */
		StdVec<StdLargeVec<Real>> diffusion_dt_; /**< array of the time derivative of diffusion species */
		StdVec<StdLargeVec<Real>> diffusion_dt_prior_; /**< array of the time derivative of diffusion species */

		DiffusionReactionParticles(SPHBody &sph_body,
								   SharedPtr<DiffusionReaction<BaseParticlesType, BaseMaterialType>> diffusion_reaction_material_ptr,
								   SharedPtr<ParticleGenerator> particle_generator_ptr = makeShared<ParticleGeneratorLattice>())
			: BaseParticlesType(sph_body, diffusion_reaction_material_ptr, particle_generator_ptr)
		{
			diffusion_reaction_material_ptr->assignDiffusionReactionParticles(this);
			species_indexes_map_ = diffusion_reaction_material_ptr->SpeciesIndexMap();
			species_diffusion_source_index_ = diffusion_reaction_material_ptr->SpeciesDiffusionSourceIndex();
			number_of_species_ = diffusion_reaction_material_ptr->NumberOfSpecies();
			species_n_.resize(number_of_species_);

			std::map<std::string, size_t>::iterator itr;
			for (itr = species_indexes_map_.begin(); itr != species_indexes_map_.end(); ++itr)
			{
				//Register a species. Note that we call a template function from a template class
				this->template registerAVariable<indexScalar, Real>(species_n_[itr->second], itr->first);
				//the scalars will be sorted if particle sorting is called
				this->template registerASortableVariable<indexScalar, Real>(itr->first);
				// add species to basic output particle data
				this->template addAVariableToWrite<indexScalar, Real>(itr->first);
			}

			number_of_diffusion_species_ = diffusion_reaction_material_ptr->NumberOfSpeciesDiffusion();
			diffusion_dt_.resize(number_of_diffusion_species_);
			diffusion_dt_prior_.resize(species_diffusion_source_index_.size());
			for (size_t m = 0; m < number_of_diffusion_species_; ++m)
			{
				//----------------------------------------------------------------------
				//	register reactive change rate terms without giving variable name
				//----------------------------------------------------------------------
				std::get<indexScalar>(this->all_particle_data_).push_back(&diffusion_dt_[m]);
				diffusion_dt_[m].resize(this->total_real_particles_, Real(0));
			}
			for (size_t m = 0; m < species_diffusion_source_index_.size(); ++m)
			{
				//----------------------------------------------------------------------
				//	register reactive change rate terms without giving variable name
				//----------------------------------------------------------------------
				diffusion_dt_prior_[m].resize(this->total_real_particles_, Real(0));
			}
		};
/*
		template<class PType, class MType>
		virtual size_t initializeAParticle(DiffusionReactionParticles<PType, MType> &fromp, size_t index)
		{
			for (itr = species_indexes_map_.begin(); itr != species_indexes_map_.end(); ++itr)
			{
				if(empty_buffer)
					species_n_[i_specie][total_real_particles_] = fromp.species_n_[i_specie][index];
				else
					species_n_[i_specie].push_back(fromp.species_n_[i_specie][index]);
			}
			BaseParticlesType::initializeAParticle(fromp.pos_, fromp.vol_);
		};*/
		virtual void destroyAParticle(size_t index){

		};

		virtual ~DiffusionReactionParticles(){};

		std::map<std::string, size_t> SpeciesIndexMap() { return species_indexes_map_; };
		StdVec<size_t> SpeciesDiffusionSourceIndex() { return species_diffusion_source_index_; };

		virtual DiffusionReactionParticles<BaseParticlesType, BaseMaterialType> *
		ThisObjectPtr() override { return this; };
	};

	/**
	 * @class ElectroPhysiologyParticles
	 * @brief A group of particles with electrophysiology particle data.
	 */
	class ElectroPhysiologyParticles
		: public DiffusionReactionParticles<SolidParticles, Solid>
	{
	public:
		ElectroPhysiologyParticles(SPHBody &sph_body,
								   SharedPtr<DiffusionReaction<SolidParticles, Solid>> diffusion_reaction_material_ptr,
								   SharedPtr<ParticleGenerator> particle_generator_ptr = makeShared<ParticleGeneratorLattice>());
		virtual ~ElectroPhysiologyParticles(){};
		virtual ElectroPhysiologyParticles *ThisObjectPtr() override { return this; };
	};
	/**
	 * @class ElectroPhysiologyReducedParticles
	 * @brief A group of reduced particles with electrophysiology particle data.
	 */
	class ElectroPhysiologyReducedParticles
		: public DiffusionReactionParticles<SolidParticles, Solid>
	{
	public:
		/** Constructor. */
		ElectroPhysiologyReducedParticles(SPHBody &sph_body,
										  SharedPtr<DiffusionReaction<SolidParticles, Solid>> diffusion_reaction_material_ptr,
										  SharedPtr<ParticleGenerator> particle_generator_ptr = makeShared<ParticleGeneratorLattice>());
		/** Destructor. */
		virtual ~ElectroPhysiologyReducedParticles(){};
		virtual ElectroPhysiologyReducedParticles *ThisObjectPtr() override { return this; };

		virtual Vecd getKernelGradient(size_t particle_index_i, size_t particle_index_j, Real dW_ij, Vecd &e_ij) override
		{
			return dW_ij * e_ij;
		};
	};
}
#endif //DIFFUSION_REACTION_PARTICLES_H