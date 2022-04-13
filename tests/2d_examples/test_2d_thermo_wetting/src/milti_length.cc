
struct SearchDepthMultiResolution
{
	Real longest = 1;
	int operator()(size_t particle_index) const { return ceil(longest); };
};

class BodyRelationInnerMultiLength : public BodyRelationInner
{
public:
	BodyRelationInnerMultiLength(RealBody &real_body, Real ratio):BodyRelationInner(real_body){};
};

class BodyRelationInnerMultiLength : public BaseBodyRelationInner
{
protected:
	Real longest_;
	SPHBodyParticlesIndex get_particle_index_;
	SearchDepthMultiResolution get_longest_search_depth_;
	NeighborRelationInner MultiLength get_inner_neighbor_;
	CellLinkedList *cell_linked_list_;
public:
	explicit BodyRelationInnerMultiLength(RealBody &real_body, Real longest)
	:longest_(longest){
		get_longest_search_depth_.longest = longest;
	};
	virtual ~BodyRelationInnerMultiLength(){};
	virtual void updateConfiguration() override;
};
void BodyRelationInnerMultiLength::updateConfiguration()
{
	resetNeighborhoodCurrentSize();
	cell_linked_list_
		->searchNeighborsByParticles(base_particles_->total_real_particles_,
									 *base_particles_, inner_configuration_,
									 get_particle_index_, get_longest_search_depth_,
									 get_inner_neighbor_);
}

	class NeighborhoodMultiLength
	{
	public:
		size_t current_size_;	/**< the current number of neighors */
		size_t allocated_size_; /**< the limit of neighors does not require memory allocation  */
		StdLargeVec<size_t> j_;	  /**< index of the neighbor particle. */
		StdLargeVec<Real> W_ij_;  /**< kernel value or particle volume contribution */
		StdLargeVec<Real> dW_ij_; /**< derivative of kernel function or inter-particle surface contribution */
		StdLargeVec<Real> r_ij_;  /**< distance between j and i. */
		StdLargeVec<Vecd> e_ij_;  /**< unit vector pointing from j to i or inter-particle surface direction */
		StdVec<Real> lengths_;
		StdVec<StdLargeVec<Real>> W_ij_n_;  /**< kernel value or particle volume contribution */
		StdVec<StdLargeVec<Real>> dW_ij_n_; /**< derivative of kernel function or inter-particle surface contribution */

		NeighborhoodMultiLength() : current_size_(0), allocated_size_(0){};
		~NeighborhoodMultiLength(){};

		void removeANeighbor(size_t neighbor_n);
	};
	class NeighborRelationInnerMultiLength : public NeighborRelationMultiLength
	{
	public:
		explicit NeighborRelationInnerMultiLength(SPHBody *body);
		void operator()(NeighborhoodMultiLength &neighborhood,
						Vecd &displacement, size_t i_index, size_t j_index) const;
	};
	NeighborRelationInnerMultiLength::NeighborRelationInnerMultiLength(SPHBody *body) : NeighborRelation()
	{
		kernel_ = body->sph_adaptation_->getKernel();
	}
	//=================================================================================================//
	void NeighborRelationInnerMultiLength::operator()(NeighborhoodMultiLength &neighborhood,
										   Vecd &displacement, size_t i_index, size_t j_index) const
	{
		Real distance = displacement.norm();
		if (distance < kernel_->CutOffRadius() && i_index != j_index)
		{
			neighborhood.current_size_ >= neighborhood.allocated_size_
				? createRelation(neighborhood, distance, displacement, j_index)
				: initializeRelation(neighborhood, distance, displacement, j_index);
			neighborhood.current_size_++;
		}
	};
	void NeighborRelationMultiLength::createRelation(NeighborhoodMultiLength &neighborhood,
										  Real &distance, Vecd &displacement, size_t j_index) const
	{
		neighborhood.j_.push_back(j_index);
		neighborhood.W_ij_.push_back(kernel_->W(distance, displacement));
		neighborhood.dW_ij_.push_back(kernel_->dW(distance, displacement));
		neighborhood.r_ij_.push_back(distance);
		neighborhood.e_ij_.push_back(displacement / (distance + TinyReal));
		neighborhood.allocated_size_++;
	}
	//=================================================================================================//
	void NeighborRelationMultiLength::initializeRelation(NeighborhoodMultiLength &neighborhood,
											  Real &distance, Vecd &displacement, size_t j_index) const
	{
		size_t current_size = neighborhood.current_size_;
		neighborhood.j_[current_size] = j_index;
		neighborhood.W_ij_[current_size] = kernel_->W(distance, displacement);
		neighborhood.dW_ij_[current_size] = kernel_->dW(distance, displacement);
		neighborhood.r_ij_[current_size] = distance;
		neighborhood.e_ij_[current_size] = displacement / (distance + TinyReal);
	}
	//=================================================================================================//
	void NeighborRelationMultiLength::createRelation(NeighborhoodMultiLength &neighborhood, Real &distance,
										  Vecd &displacement, size_t j_index, Real i_h_ratio, Real h_ratio_min) const
	{
		neighborhood.j_.push_back(j_index);
		Real weight = distance < kernel_->CutOffRadius(i_h_ratio) ? kernel_->W(i_h_ratio, distance, displacement) : 0.0;
		neighborhood.W_ij_.push_back(weight);
		neighborhood.dW_ij_.push_back(kernel_->dW(h_ratio_min, distance, displacement));
		neighborhood.r_ij_.push_back(distance);
		neighborhood.e_ij_.push_back(displacement / (distance + TinyReal));
		neighborhood.allocated_size_++;
	}
	//=================================================================================================//
	void NeighborRelationMultiLength::
		initializeRelation(NeighborhoodMultiLength &neighborhood, Real &distance,
						   Vecd &displacement, size_t j_index, Real i_h_ratio, Real h_ratio_min) const
	{
		size_t current_size = neighborhood.current_size_;
		neighborhood.j_[current_size] = j_index;
		neighborhood.W_ij_[current_size] = distance < kernel_->CutOffRadius(i_h_ratio)
											   ? kernel_->W(i_h_ratio, distance, displacement)
											   : 0.0;
		neighborhood.dW_ij_[current_size] = kernel_->dW(h_ratio_min, distance, displacement);
		neighborhood.r_ij_[current_size] = distance;
		neighborhood.e_ij_[current_size] = displacement / (distance + TinyReal);
	}
