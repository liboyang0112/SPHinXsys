#ifndef __MULTI_LENGTH
#define __MULTI_LENGTH
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

#endif