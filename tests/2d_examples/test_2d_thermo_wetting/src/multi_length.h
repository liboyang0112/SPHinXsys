
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

	class BodyRelationContactMultiLength : public BaseBodyRelationContact
	{
	protected:
		SPHBodyParticlesIndex get_particle_index_;

		void initialization();

	public:
		BodyRelationContact(SPHBody &sph_body, RealBodyVector contact_bodies);
		BodyRelationContact(SPHBody &sph_body, BodyPartVector contact_body_parts);
		virtual ~BodyRelationContact(){};
		virtual void updateConfiguration() override;
	};


	class NeighborRelationContact : public NeighborRelationMultiLength
	{
	public:
		NeighborRelationContact(SPHBody *body, SPHBody *contact_body);
		virtual ~NeighborRelationContact(){};
		void operator()(Neighborhood &neighborhood,
						Vecd &displacement, size_t i_index, size_t j_index) const;
	};


class BodyRelationContactMultiLength : public BaseBodyRelationContact
{
protected:
	StdVec<Real> &lengths_;
	UniquePtrVectorKeeper<SearchDepthMultiResolution> search_depth_multi_resolution_ptr_vector_keeper_;
	UniquePtrVectorKeeper<NeighborRelationContact> neighbor_relation_contact_ptr_vector_keeper_;
protected:
	StdVec<CellLinkedList *> target_cell_linked_lists_;
	StdVec<SearchDepthMultiResolution *> get_search_depths_;
	StdVec<NeighborRelationContactMultiLength *> get_contact_neighbors_;
	virtual void resetNeighborhoodCurrentSize();
public:
	RealBodyVector contact_bodies_;
	ContatcParticleConfiguration contact_configuration_; /**< Configurations for particle interaction between bodies. */
	BaseBodyRelationContact(SPHBody &sph_body, RealBodyVector contact_bodies);
	BaseBodyRelationContact(SPHBody &sph_body, BodyPartVector contact_body_parts);
	virtual ~BaseBodyRelationContact(){};
	virtual void updateConfigurationMemories() override;
};


	class ComplexBodyRelation : public SPHBodyRelation
	{
	private:
		UniquePtrKeeper<BaseBodyRelationInner> base_body_relation_inner_ptr_keeper_;
		UniquePtrKeeper<BaseBodyRelationContact> base_body_relation_contact_ptr_keeper_;

	public:
		BaseBodyRelationInner &inner_relation_;
		BaseBodyRelationContact &contact_relation_;
		RealBodyVector contact_bodies_;
		ParticleConfiguration &inner_configuration_;
		ContatcParticleConfiguration &contact_configuration_;

		ComplexBodyRelation(BaseBodyRelationInner &inner_relation, BaseBodyRelationContact &contact_relation);
		ComplexBodyRelation(RealBody &real_body, RealBodyVector contact_bodies);
		ComplexBodyRelation(BaseBodyRelationInner &inner_relation, RealBodyVector contact_bodies);
		ComplexBodyRelation(RealBody &real_body, BodyPartVector contact_body_parts);
		virtual ~ComplexBodyRelation(){};

		virtual void updateConfigurationMemories() override;
		virtual void updateConfiguration() override;
	};