
#include "sphinxsys.h"

using namespace SPH;

struct SearchDepthMultiLength
{
	int search_depth_ = 1;
	SearchDepthMultiLength(StdVec<Real> &lengths){
		if(lengths.size())
			search_depth_ = ceil(*std::max_element(lengths.begin(),lengths.end()));
	}
	int operator()(size_t particle_index) const { 
		return search_depth_; 
	};
};


struct SearchDepthMultiResolution
{
	int search_depth_;
	Real longest = 1;
	SearchDepthMultiResolution(SPHBody &sph_body, CellLinkedList *target_cell_linked_list)
		: search_depth_(1)
	{
		Real inv_grid_spacing_ = 1.0 / target_cell_linked_list->GridSpacing();
		Kernel *kernel_ = sph_body.sph_adaptation_->getKernel();
		search_depth_ = ceil(kernel_->CutOffRadius() * inv_grid_spacing_);
	};
	int operator()(size_t particle_index) const { return search_depth_; };
};

class NeighborRelationMultiLength: public NeighborRelation
{
protected:
	StdVec<Real> &lengths_;
	Real longest;
	//----------------------------------------------------------------------
	//	Below are for constant smoothing length.
	//----------------------------------------------------------------------
	void createRelation(Neighborhood &neighborhood, Real &distance,
						Vecd &displacement, size_t j_index) const override;
	void initializeRelation(Neighborhood &neighborhood, Real &distance,
							Vecd &displacement, size_t j_index) const override;
public:
	NeighborRelationMultiLength(StdVec<Real> &lengths) : NeighborRelation(), lengths_(lengths){
		longest = *std::max_element(lengths.begin(),lengths.end());
	};
	virtual ~NeighborRelationMultiLength(){};
};


class NeighborRelationInnerMultiLength : public NeighborRelationMultiLength
{
public:
	NeighborRelationInnerMultiLength(SPHBody *body, StdVec<Real> &lengths);
	void operator()(Neighborhood &neighborhood,
					Vecd &displacement, size_t i_index, size_t j_index) const;
};

class BodyRelationInnerMultiLength : public BaseBodyRelationInner
{
protected:
	StdVec<Real> &lengths_;
	SPHBodyParticlesIndex get_particle_index_;
	SearchDepthMultiLength get_longest_search_depth_;
	NeighborRelationInnerMultiLength get_inner_neighbor_;
	CellLinkedList *cell_linked_list_;
public:
	explicit BodyRelationInnerMultiLength(RealBody &real_body, StdVec<Real> &lengths)
	:BaseBodyRelationInner(real_body), lengths_(lengths), get_inner_neighbor_(&real_body, lengths),
	cell_linked_list_(DynamicCast<CellLinkedList>(this, real_body.cell_linked_list_)),
	get_longest_search_depth_(lengths){};
	virtual ~BodyRelationInnerMultiLength(){};
	virtual void updateConfiguration() override;
};
/*
class BodyRelationContactMultiLength : public BaseBodyRelationContact
{
protected:
	StdVec<Real> &lengths_;
	UniquePtrVectorKeeper<SearchDepthMultiResolution> search_depth_multi_resolution_ptr_vector_keeper_;
	UniquePtrVectorKeeper<NeighborRelationContact> neighbor_relation_contact_ptr_vector_keeper_;
protected:
	StdVec<CellLinkedList *> target_cell_linked_lists_;
	StdVec<SearchDepthMultiResolution *> get_search_depths_;
	StdVec<NeighborRelationContact *> get_contact_neighbors_;
	virtual void resetNeighborhoodCurrentSize();
public:
	RealBodyVector contact_bodies_;
	ContatcParticleConfiguration contact_configuration_; /**< Configurations for particle interaction between bodies. */
/*	BaseBodyRelationContact(SPHBody &sph_body, RealBodyVector contact_bodies);
	BaseBodyRelationContact(SPHBody &sph_body, BodyPartVector contact_body_parts);
	virtual ~BaseBodyRelationContact(){};
	virtual void updateConfigurationMemories() override;
};*/