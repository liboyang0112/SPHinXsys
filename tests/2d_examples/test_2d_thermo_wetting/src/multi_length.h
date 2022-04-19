
#include "sphinxsys.h"

using namespace SPH;

struct SearchDepthMultiLength
{
	Real longest = 1;
	int operator()(size_t particle_index) const { 
		if(longest>1) return ceil(longest); 
		return 1;
	};
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
	NeighborRelationInnerMultiLength(StdVec<Real> &lengths) : NeighborRelationMultiLength(lengths){}
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
	explicit BodyRelationInnerMultiLength(RealBody &real_body, StdVec<Real> lengths)
	:BaseBodyRelationInner(real_body), lengths_(lengths), get_inner_neighbor_(lengths){
		get_longest_search_depth_.longest = *std::max_element(lengths.begin(),lengths.end());
	};
	virtual ~BodyRelationInnerMultiLength(){};
	virtual void updateConfiguration() override
	{
		resetNeighborhoodCurrentSize();
		cell_linked_list_
			->searchNeighborsByParticles(base_particles_->total_real_particles_,
										 *base_particles_, inner_configuration_,
										 get_particle_index_, get_longest_search_depth_,
										 get_inner_neighbor_);
	}
};