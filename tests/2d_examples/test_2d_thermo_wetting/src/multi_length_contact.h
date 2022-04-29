#ifndef __MULTI_LENGTH_CONTACT
#define __MULTI_LENGTH_CONTACT
#include "multi_length.h"

struct SearchDepthMultiLengthAndResolution
{
	int search_depth_ = 1;
	SearchDepthMultiLengthAndResolution(SPHBody &sph_body, CellLinkedList *target_cell_linked_list, StdVec<Real> &lengths)
		: search_depth_(1)
	{
		search_depth_ = ceil(sph_body.sph_adaptation_->getKernel()->CutOffRadius() 
		/ target_cell_linked_list->GridSpacing()
		* *std::max_element(lengths.begin(),lengths.end()));
	};
	int operator()(size_t particle_index) const { return search_depth_; };
};

class NeighborRelationContactMultiLength : public NeighborRelationMultiLength
{
public:
	NeighborRelationContactMultiLength(SPHBody *body, SPHBody *contact_body, StdVec<Real> &lengths);
	virtual ~NeighborRelationContactMultiLength(){};
	void operator()(Neighborhood &neighborhood,
					Vecd &displacement, size_t i_index, size_t j_index) const;
};

class BodyRelationContactMultiLength : public BaseBodyRelationContact
{
protected:
	StdVec<Real> &lengths_;
	SPHBodyParticlesIndex get_particle_index_;
	StdVec<unique_ptr<SearchDepthMultiLengthAndResolution>> get_search_depths_multi_lengths_;
	StdVec<unique_ptr<NeighborRelationContactMultiLength>> get_contact_neighbors_multi_lengths_;
	void initialization();
public:
	BodyRelationContactMultiLength(SPHBody &sph_body, RealBodyVector contact_bodies, StdVec<SPH::Real> &lengths);
	BodyRelationContactMultiLength(SPHBody &sph_body, BodyPartVector contact_body_parts, StdVec<SPH::Real> &lengths);
	virtual ~BodyRelationContactMultiLength(){};
	virtual void updateConfiguration() override;
};
#endif
