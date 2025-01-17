
#include "multi_length.h"

#include "/home/boyang/softwares/SPHinXsys/SPHinXsys/SPHINXsys/src/for_2D_build/meshes/cell_linked_list.hpp"

void NeighborRelationMultiLength::createRelation(Neighborhood &neighborhood,
									  Real &distance, Vecd &displacement, size_t j_index) const
{
	neighborhood.j_.push_back(j_index);
	for(int i=0; i<lengths_.size(); i++){
		if(distance<kernel_->CutOffRadius(1./lengths_[i])){
			neighborhood.W_ij_n_[i].push_back(kernel_->W(1./lengths_[i], distance, displacement));
			neighborhood.dW_ij_n_[i].push_back(kernel_->dW(1./lengths_[i], distance, displacement));
		}else{
			neighborhood.W_ij_n_[i].push_back(0);
			neighborhood.dW_ij_n_[i].push_back(0);
		}
	}
	if(distance < kernel_->CutOffRadius()){
		neighborhood.W_ij_.push_back(kernel_->W(distance, displacement));
		neighborhood.dW_ij_.push_back(kernel_->dW(distance, displacement));
	}else{
		neighborhood.W_ij_.push_back(0);
		neighborhood.dW_ij_.push_back(0);
	}
	neighborhood.r_ij_.push_back(distance);
	neighborhood.e_ij_.push_back(displacement / (distance + TinyReal));
	neighborhood.allocated_size_++;
}
//=================================================================================================//
void NeighborRelationMultiLength::initializeRelation(Neighborhood &neighborhood,
										  Real &distance, Vecd &displacement, size_t j_index) const
{
	size_t current_size = neighborhood.current_size_;
	neighborhood.j_[current_size] = j_index;
	for(int i=0; i<lengths_.size(); i++){
		if(distance<kernel_->CutOffRadius(1./lengths_[i])){
			neighborhood.W_ij_n_[i][current_size] = kernel_->W(1./lengths_[i], distance, displacement);
			neighborhood.dW_ij_n_[i][current_size] = kernel_->dW(1./lengths_[i], distance, displacement);
		}else{
			neighborhood.W_ij_n_[i][current_size] = 0;
			neighborhood.dW_ij_n_[i][current_size] = 0;
		}
	}
	if(distance < kernel_->CutOffRadius()){
		neighborhood.W_ij_[current_size] = kernel_->W(distance, displacement);
		neighborhood.dW_ij_[current_size] = kernel_->dW(distance, displacement);
	}else{
		neighborhood.W_ij_[current_size] = 0;
		neighborhood.dW_ij_[current_size] = 0;
	}
	neighborhood.r_ij_[current_size] = distance;
	neighborhood.e_ij_[current_size] = displacement / (distance + TinyReal);
}

void NeighborRelationInnerMultiLength::operator()(Neighborhood &neighborhood,
									   Vecd &displacement, size_t i_index, size_t j_index) const
{
	Real distance = displacement.norm();
	if(neighborhood.dW_ij_n_.size() == 0){
		neighborhood.dW_ij_n_.resize(lengths_.size());
		neighborhood.W_ij_n_.resize(lengths_.size());
	}
	if (distance < kernel_->CutOffRadius(1./longest) && i_index != j_index)
	{
		neighborhood.current_size_ >= neighborhood.allocated_size_
			? createRelation(neighborhood, distance, displacement, j_index)
			: initializeRelation(neighborhood, distance, displacement, j_index);
		neighborhood.current_size_++;
	}
};

NeighborRelationInnerMultiLength::NeighborRelationInnerMultiLength(SPHBody *body, StdVec<Real> &lengths) : NeighborRelationMultiLength(lengths)
{
	kernel_ = body->sph_adaptation_->getKernel();
}

void BodyRelationInnerMultiLength::updateConfiguration()
{
	resetNeighborhoodCurrentSize();
	cell_linked_list_
		->searchNeighborsByParticles(base_particles_->total_real_particles_,
									 *base_particles_, inner_configuration_,
									 get_particle_index_, get_longest_search_depth_,
									 get_inner_neighbor_);
}