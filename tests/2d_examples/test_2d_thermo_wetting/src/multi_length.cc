
#include "multi_length.h"

//=================================================================================================//

void NeighborRelationMultiLength::createRelation(Neighborhood &neighborhood,
									  Real &distance, Vecd &displacement, size_t j_index) const
{
	neighborhood.j_.push_back(j_index);
	for(Real &x : lengths_){
		if(distance<kernel_->CutOffRadius(1./x)){
			neighborhood.W_ij_.push_back(kernel_->W(1./x, distance, displacement));
			neighborhood.dW_ij_.push_back(kernel_->dW(1./x, distance, displacement));
		}else{
			neighborhood.W_ij_.push_back(0);
			neighborhood.dW_ij_.push_back(0);
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
	for(Real &x : lengths_){
		if(distance<kernel_->CutOffRadius(1./x)){
			neighborhood.W_ij_[current_size] = kernel_->W(1./x, distance, displacement);
			neighborhood.dW_ij_[current_size] = kernel_->dW(1./x, distance, displacement);
		}else{
			neighborhood.W_ij_[current_size] = 0;
			neighborhood.dW_ij_[current_size] = 0;
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
//=================================================================================================//
