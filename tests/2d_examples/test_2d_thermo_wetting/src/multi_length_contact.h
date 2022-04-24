#include "multi_length.h"

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