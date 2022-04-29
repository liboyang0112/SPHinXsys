#include "multi_length_contact.h"
#include "/home/boyang/softwares/SPHinXsys/SPHinXsys/SPHINXsys/src/for_2D_build/meshes/cell_linked_list.hpp"
	//=================================================================================================//
	NeighborRelationContactMultiLength::
		NeighborRelationContactMultiLength(SPHBody *body, SPHBody *contact_body, StdVec<Real> &lengths) : NeighborRelationMultiLength(lengths)
	{
		Kernel *source_kernel = body->sph_adaptation_->getKernel();
		Kernel *target_kernel = contact_body->sph_adaptation_->getKernel();
		kernel_ = source_kernel->SmoothingLength() > target_kernel->SmoothingLength() ? source_kernel : target_kernel;
	}
	void NeighborRelationContactMultiLength::operator()(Neighborhood &neighborhood,
											 Vecd &displacement, size_t i_index, size_t j_index) const
	{
		Real distance = displacement.norm();
		if(neighborhood.dW_ij_n_.size() == 0){
			neighborhood.dW_ij_n_.resize(lengths_.size());
			neighborhood.W_ij_n_.resize(lengths_.size());
		}
		if (distance < kernel_->CutOffRadius(1./longest))
		{
			neighborhood.current_size_ >= neighborhood.allocated_size_
				? createRelation(neighborhood, distance, displacement, j_index)
				: initializeRelation(neighborhood, distance, displacement, j_index);
			neighborhood.current_size_++;
		}
	};

	BodyRelationContactMultiLength::BodyRelationContactMultiLength(SPHBody &sph_body, RealBodyVector contact_sph_bodies, StdVec<SPH::Real> &lengths)
		: BaseBodyRelationContact(sph_body, contact_sph_bodies), lengths_(lengths)
	{
		initialization();
	}
	//=================================================================================================//
	BodyRelationContactMultiLength::BodyRelationContactMultiLength(SPHBody &sph_body, BodyPartVector contact_body_parts, StdVec<SPH::Real> &lengths)
		: BaseBodyRelationContact(sph_body, contact_body_parts), lengths_(lengths)
	{
		initialization();
	}
	//=================================================================================================//
	void BodyRelationContactMultiLength::initialization()
	{
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			CellLinkedList *target_cell_linked_list =
				DynamicCast<CellLinkedList>(this, contact_bodies_[k]->cell_linked_list_);
			target_cell_linked_lists_.push_back(target_cell_linked_list);
			get_search_depths_multi_lengths_.push_back(
				make_unique<SearchDepthMultiLengthAndResolution>(*sph_body_, target_cell_linked_list, lengths_));
			get_contact_neighbors_multi_lengths_.push_back(
				make_unique<NeighborRelationContactMultiLength>(sph_body_, contact_bodies_[k], lengths_));
		}
	}
	//=================================================================================================//
	void BodyRelationContactMultiLength::updateConfiguration()
	{
		resetNeighborhoodCurrentSize();
		size_t total_real_particles = base_particles_->total_real_particles_;
		for (size_t k = 0; k != contact_bodies_.size(); ++k)
		{
			target_cell_linked_lists_[k]
				->searchNeighborsByParticles(total_real_particles,
											 *base_particles_, contact_configuration_[k],
											 get_particle_index_, *get_search_depths_multi_lengths_[k],
											 *get_contact_neighbors_multi_lengths_[k]);
		}
	}
