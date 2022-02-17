/**
 * @file 	case.h
 * @brief 	Numerical parameters and body defination for 2D two-phase dambreak flow.
 * @author 	Chi Zhang and Xiangyu Hu
 */
 /**
  * @brief 	SPHinXsys Library.
  */
#include "sphinxsys.h"
  /**
 * @brief Namespace cite here.
 */
using namespace SPH;
/**
 * @brief Basic geometry parameters and numerical setup.
 */
const Real k_B = 3.285e4;
Real DL = 2.0; 							/**< Tank length. */
Real DH = 1.0; 							/**< Tank height. */
Real particle_spacing_ref = DL / 60.0; 	/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; 	/**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));


/** Material properties. */
Real diffusion_coff = 1.0e-3;
Real bias_coff = 0.0;
Real alpha = Pi / 6.0;
Vec2d bias_direction(cos(alpha), sin(alpha));
/**
*@brief Temperatures.
*/
Real phi_upper_wall = 00.0;
Real phi_lower_wall = 40.0;
Real phi_side_wall = 0;
Real phi_fluid_initial = 0.6;
Real phi_gas_initial = 0.6;

/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.4;					/**< Reference density of water. */
Real rho0_a = 1.0e-3;				/**< Reference density of air. */
Real gravity_g = 9.8/2.81/1e9;				/**< Gravity force of fluid. */
Real U_max = 1.0;					/**< Characteristic velocity. */
Real c_f = 10.0 * U_max;			/**< Reference sound speed. */
Real mu_f = 5.0e-2;					/**< Water viscosity. */
Real mu_a = 5.0e-5;					/**< Air visocsity. */
Real contact_angle = (150.0 / 180.0) * 3.1415926; 	/**< Contact angle with Wall. */
Real tension_force = 0.008;
/** vdW properties. */
Real gamma_vdw=7./5;  //(dim+2)/dim
Real rho0_water=1./0.5824;  // molecule size
Real alpha_water=7.438e4;  // attraction force
Real molmass_water=18;
Real molmass_air=27;
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(0.375 * DL, 0.0));
	water_block_shape.push_back(Vecd(0.375 * DL, 0.35 * DH));
	water_block_shape.push_back(Vecd(0.625 * DL, 0.35 * DH));
	water_block_shape.push_back(Vecd(0.625 *DL, 0.0));
	water_block_shape.push_back(Vecd(0.375 * DL, 0.0));
	return water_block_shape;
}
/** create outer wall shape */
std::vector<Vecd> createOuterWallShape()
{
	std::vector<Vecd> outer_wall_shape;
	outer_wall_shape.push_back(Vecd(-BW, -BW));
	outer_wall_shape.push_back(Vecd(-BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, -BW));
	outer_wall_shape.push_back(Vecd(-BW, -BW));

	return outer_wall_shape;
}
/**
* @brief create inner wall shape
*/
std::vector<Vecd> createInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(0.0, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, DH));
	inner_wall_shape.push_back(Vecd(DL, DH));
	inner_wall_shape.push_back(Vecd(DL, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, 0.0));

	return inner_wall_shape;
}
/**
 * @brief 	Case dependent material properties definition.
 */



/**
 * application dependent solid body initial condition
 */
class ThermosolidBodyInitialCondition
	: public  DiffusionReactionInitialCondition<SolidBody, SolidParticles, Solid>
{
protected:
	size_t phi_;

	void Update(size_t index_i, Real dt) override
	{
		
		if (-BW <= pos_n_[index_i][1] && pos_n_[index_i][1] <= 0.0)
		{
			species_n_[phi_][index_i] = phi_lower_wall;
		}
		else if (DH <= pos_n_[index_i][1] && pos_n_[index_i][1] <= DH+BW)
		{
			species_n_[phi_][index_i] = phi_upper_wall;
		}else{
			species_n_[phi_][index_i] = phi_side_wall;
		}
		
	};
public: 
	ThermosolidBodyInitialCondition(SolidBody &diffusion_solid_body)
		: DiffusionReactionInitialCondition<SolidBody, SolidParticles, Solid>(diffusion_solid_body) {
		phi_ = material_->SpeciesIndexMap()["Temperature"];
	};
};

class vdWFluid : public CompressibleFluid
{
protected:
	Real alpha_, molmass_, rho_max_; // alpha is the attraction force
public:
	explicit vdWFluid(Real rho0, Real rho_max, Real gamma, Real alpha, Real molmass, Real mu = 0.0) 
	//rho0 is the density maximum
		: CompressibleFluid(rho0, gamma, mu), alpha_(alpha), molmass_(molmass), rho_max_(rho_max)
	{
		material_type_ = "vdWFluid";
	};
	virtual ~vdWFluid(){};
	Real getAlpha(){ return alpha_; }
	Real HeatCapacityRatio() { return gamma_; };
	virtual Real DensityFromPT(Real p, Real T) {};
	virtual Real getPressure(Real rho, Real T) override;
	virtual Real getSoundSpeed(Real p, Real rho) override;
	virtual vdWFluid *ThisObjectPtr() override { return this; };
};

Real vdWFluid::getPressure(Real rho, Real T)
{
	Real ratio = 1./(1-rho/rho_max_);
	return rho*k_B*T*ratio-alpha_*rho*rho;
}

Real vdWFluid::getSoundSpeed(Real p, Real rho)
{
	Real ratio = 1./(1-rho/rho_max_);
	return sqrt(gamma_ * ratio * (alpha_ * rho + p / rho) - 2*alpha_*rho);
}



/**
 * application dependent fluid body initial condition
 */
class ThermofluidBodyInitialCondition
	: public  DiffusionReactionInitialCondition< FluidBody, CompressibleFluidParticles, vdWFluid>
{
protected:
	size_t phi_;

	void Update(size_t index_i, Real dt) override
	{

		if (0 <= pos_n_[index_i][1] && pos_n_[index_i][1] <= DH)
		{
			species_n_[phi_][index_i] = phi_fluid_initial;
		}

	};
public:
	ThermofluidBodyInitialCondition(FluidBody &diffusion_fluid_body)
		: DiffusionReactionInitialCondition<FluidBody, CompressibleFluidParticles, vdWFluid >(diffusion_fluid_body) {
		phi_ = material_->SpeciesIndexMap()["Temperature"];
	};
};

/**
 * application dependent gas body initial condition
 */
class ThermogasBodyInitialCondition
	: public  DiffusionReactionInitialCondition< FluidBody, CompressibleFluidParticles, vdWFluid>
{
protected:
	size_t phi_;

	void Update(size_t index_i, Real dt) override
	{

		if (0 <= pos_n_[index_i][1] && pos_n_[index_i][1] <= DH)
		{
			species_n_[phi_][index_i] = phi_gas_initial;
		}

	};
public:
	ThermogasBodyInitialCondition(FluidBody &diffusion_fluid_body)
		: DiffusionReactionInitialCondition<FluidBody, CompressibleFluidParticles, vdWFluid >(diffusion_fluid_body) {
		phi_ = material_->SpeciesIndexMap()["Temperature"];
	};
};

//----------------------------------------------------------------------
//	An observer particle generator
//----------------------------------------------------------------------
class ObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	ObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		/** A measuring point at the center of the channel */
		Vec2d point_coordinate(0.0, DH * 0.5);
		positions_volumes_.push_back(std::make_pair(point_coordinate, 0.0));
	}
};

/**
 *Set thermal relaxation between different bodies 
 */
class ThermalRelaxationComplex
	: public RelaxationOfAllDiffusionSpeciesRK2<
		  FluidBody, CompressibleFluidParticles, vdWFluid,
		  RelaxationOfAllDiffussionSpeciesComplex<
			  FluidBody, CompressibleFluidParticles, vdWFluid, SolidBody, SolidParticles, Solid>,
		  ComplexBodyRelation>
{
public:
	explicit ThermalRelaxationComplex(ComplexBodyRelation &body_complex_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_complex_relation){};
	virtual ~ThermalRelaxationComplex(){};
};


class ThermalRelaxationComplexWA
	: public RelaxationOfAllDiffusionSpeciesRK2<
		  FluidBody, FluidParticles, vdWFluid,
		  RelaxationOfAllDiffussionSpeciesComplex<
			  FluidBody, FluidParticles, vdWFluid, FluidBody, FluidParticles, vdWFluid>,
		  ComplexBodyRelation>
{
public:
	explicit ThermalRelaxationComplexWA(ComplexBodyRelation &body_complex_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_complex_relation){};
	virtual ~ThermalRelaxationComplexWA(){};
};

/**
*@brief 	Water body definition.
*/
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &sph_system, const string &body_name)
		: FluidBody(sph_system, body_name)//, makeShared<SPHAdaptation>(1.3, 1))
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
/**
*@brief 	Air body definition.
*/
class AirBlock : public FluidBody
{
public:
	AirBlock(SPHSystem &sph_system, const std::string &body_name)
		: FluidBody(sph_system, body_name)//, makeShared<SPHAdaptation>(1.3, 1))
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createInnerWallShape(), ShapeBooleanOps::add);
		multi_polygon.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};

/**
 * @brief 	Wall boundary body definition.
 */

class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem& sph_system, string body_name)
		: SolidBody(sph_system, body_name)//, makeShared<SPHAdaptation>(1.3, 1))
	{
		/** Geomtry definition. */
		std::vector<Vecd> outer_shape = createOuterWallShape();
		std::vector<Vecd> inner_shape = createInnerWallShape();
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(outer_shape, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(inner_shape, ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};

/**
 * @brief 	Case dependent material properties definition.
 */
/**
 * Setup heat conduction material properties for diffusion fluid body 
 */
/*
class WaterMaterial
	: public DiffusionReaction<FluidParticles, WeaklyCompressibleFluid>
{
public:
	WaterMaterial()
		: DiffusionReaction<FluidParticles, WeaklyCompressibleFluid>({"Phi"}, rho0_f, c_f, mu_f)
	{
		initializeAnDiffusion<DirectionalDiffusion>("Phi", "Phi", diffusion_coff, bias_coff, bias_direction);
	};
};
 */
class WaterMaterial
	: public DiffusionReaction<CompressibleFluidParticles, vdWFluid>
{
public:
	WaterMaterial()
		: DiffusionReaction<CompressibleFluidParticles, vdWFluid>({"Temperature"},
			rho0_f, rho0_water, gamma_vdw, alpha_water, molmass_water, mu_f)
	{
		initializeAnDiffusion<DirectionalDiffusion>("Temperature", "Temperature",
			diffusion_coff, bias_coff, bias_direction);
		initializeASource("Temperature");
	};
};

class AirMaterial
	: public DiffusionReaction<CompressibleFluidParticles, vdWFluid>
{
public:
	AirMaterial()
		: DiffusionReaction<CompressibleFluidParticles, vdWFluid>({"Temperature"},
			1e-3, 1.4, gamma_vdw, 0, molmass_air, mu_a)
	{
		initializeAnDiffusion<DirectionalDiffusion>("Temperature", "Temperature",
			diffusion_coff, bias_coff, bias_direction);
		initializeASource("Temperature");
	};
};
/*
class AirMaterial
	: public DiffusionReaction<FluidParticles, WeaklyCompressibleFluid>
{
public:
	AirMaterial()
		: DiffusionReaction<FluidParticles, WeaklyCompressibleFluid>({"Phi"}, rho0_a, c_f, mu_a)
	{
		initializeAnDiffusion<DirectionalDiffusion>("Phi", "Phi", diffusion_coff, bias_coff, bias_direction);
		initializeASource("Phi");
	};
};
*/
/**
 * Setup heat conduction material properties for diffusion solid body
 */

class WallMaterial
	: public DiffusionReaction<SolidParticles, Solid>
{
public:
	WallMaterial()
		: DiffusionReaction<SolidParticles, Solid>({"Temperature"})
	{
		initializeAnDiffusion<DirectionalDiffusion>("Temperature", "Temperature", diffusion_coff, bias_coff, bias_direction);
		initializeASource("Temperature");
	};
};


class HeatSource // External heat source, laser for example
{
public:
	HeatSource(){};
	~HeatSource(){};
	Real InducedHeating(Vecd& position){return 1;};
};

template <class BodyType, class BaseParticlesType, class BaseMaterialType>
class DiffusionSourceInitialization // Similar to TimeStepInitialization
	: public ParticleDynamicsSimple,
	  public DiffusionReactionSimpleData<BodyType, BaseParticlesType, BaseMaterialType>
{
private:
	UniquePtrKeeper<HeatSource> heat_source_ptr_keeper_;
	size_t source_index_;
public:
	DiffusionSourceInitialization(SPHBody &sph_body, HeatSource &heat_source, size_t source_index);
	virtual ~DiffusionSourceInitialization(){};

protected:
	StdLargeVec<Vecd> &pos_n_;
	StdVec<StdLargeVec<Real>>  &diffusion_dt_prior_;
	HeatSource *heat_source_;
	virtual void setupDynamics(Real dt = 0.0) override;
	virtual void Update(size_t index_i, Real dt = 0.0) override;
};

template <class BodyType, class BaseParticlesType, class BaseMaterialType>
DiffusionSourceInitialization<BodyType, BaseParticlesType, BaseMaterialType>::DiffusionSourceInitialization(SPHBody &sph_body, HeatSource &heat_source, size_t source_index)
	: ParticleDynamicsSimple(sph_body), DiffusionReactionSimpleData<BodyType, BaseParticlesType, BaseMaterialType>(sph_body),
	  pos_n_(this->particles_->pos_n_), diffusion_dt_prior_(this->particles_->diffusion_dt_prior_),
	  heat_source_(&heat_source), source_index_(source_index) {}
//=================================================================================================//

template <class BodyType, class BaseParticlesType, class BaseMaterialType>
void DiffusionSourceInitialization<BodyType, BaseParticlesType, BaseMaterialType>::setupDynamics(Real dt)
{
	this->particles_->total_ghost_particles_ = 0;
}
//=================================================================================================//

template <class BodyType, class BaseParticlesType, class BaseMaterialType>
void DiffusionSourceInitialization<BodyType, BaseParticlesType, BaseMaterialType>::Update(size_t index_i, Real dt)
{
	diffusion_dt_prior_[source_index_][index_i] = heat_source_->InducedHeating(pos_n_[index_i]);
}


template <class BodyType, class BaseParticlesType, class BaseMaterialType>
class StressTensorHeatSource  // The Heat from viscous force and pressure, only inner forces are included.
	: public InteractionDynamics,
	  public DiffusionReactionInnerData<BodyType, BaseParticlesType, BaseMaterialType>
{
private:
	size_t source_index_;
public:
	StressTensorHeatSource(BaseBodyRelationInner &inner_relation, size_t source_index_);
	virtual ~StressTensorHeatSource(){};
	NoRiemannSolver riemann_solver_;
protected:
	Real mu_;
	Real smoothing_length_;
	StdLargeVec<Real> &Vol_, &rho_n_, &p_;
	StdVec<StdLargeVec<Real>> &diffusion_dt_prior_;
	StdLargeVec<Vecd> &vel_n_;
	virtual void Interaction(size_t index_i, Real dt = 0.0) override;
};

template <class BodyType, class BaseParticlesType, class BaseMaterialType>
StressTensorHeatSource<BodyType, BaseParticlesType, BaseMaterialType>::StressTensorHeatSource(BaseBodyRelationInner &inner_relation, size_t source_index)
	: InteractionDynamics(*inner_relation.sph_body_),
	  DiffusionReactionInnerData<BodyType, BaseParticlesType, BaseMaterialType>(inner_relation),
	  Vol_(this->particles_->Vol_), rho_n_(this->particles_->rho_n_), p_(this->particles_->p_),
	  vel_n_(this->particles_->vel_n_),
	  mu_(this->material_->ReferenceViscosity()),
	  smoothing_length_(sph_adaptation_->ReferenceSmoothingLength()),
	  diffusion_dt_prior_(this->particles_->diffusion_dt_prior_),
	  source_index_(source_index),
	  riemann_solver_(*(this->material_), *(this->material_)) {}
//=================================================================================================//

template <class BodyType, class BaseParticlesType, class BaseMaterialType>
void StressTensorHeatSource<BodyType, BaseParticlesType, BaseMaterialType>::Interaction(size_t index_i, Real dt)
{
	Real rho_i = rho_n_[index_i];
	Real internalEnergyIncrease(0);
	const Vecd &vel_i = vel_n_[index_i];
	FluidState state_i(rho_n_[index_i], vel_n_[index_i], p_[index_i]);

	Vecd vel_derivative(0);
	Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
	for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
	{
		size_t index_j = inner_neighborhood.j_[n];
		Real dW_ij = inner_neighborhood.dW_ij_[n];
		Vecd& e_ij = inner_neighborhood.e_ij_[n];
		Vecd vij = vel_i - vel_n_[index_j];
		FluidState state_j(rho_n_[index_j], vel_n_[index_j], p_[index_j]);
		Real p_star = riemann_solver_.getPStar(state_i, state_j, e_ij);
		//pressure
		internalEnergyIncrease += dot(2.0 * p_star * Vol_[index_j] * dW_ij * e_ij / rho_i,
			vij);
		//viscous
		vel_derivative = vij / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
		internalEnergyIncrease += dot(2.0 * mu_ * vel_derivative * Vol_[index_j] * dW_ij / rho_i,
			vij);
	}
	diffusion_dt_prior_[source_index_][index_i] += internalEnergyIncrease*2/3/k_B;
}


template <class BodyType, class BaseParticlesType, class BaseMaterialType>
class vdWAttractionHeatSource // The heat from vdW fluid "a" term
	: public ParticleDynamicsSimple,
	  public DiffusionReactionSimpleData<BodyType, BaseParticlesType, BaseMaterialType>
{
private:
	size_t source_index_;
public:
	vdWAttractionHeatSource(SPHBody &sph_body, size_t source_index);
	virtual ~vdWAttractionHeatSource(){};

protected:
	StdLargeVec<Real> &drho_dt_;
	StdVec<StdLargeVec<Real>>  &diffusion_dt_prior_;
	virtual void Update(size_t index_i, Real dt = 0.0) override;
};

template <class BodyType, class BaseParticlesType, class BaseMaterialType>
vdWAttractionHeatSource<BodyType, BaseParticlesType, BaseMaterialType>::vdWAttractionHeatSource(SPHBody &sph_body, size_t source_index)
	: ParticleDynamicsSimple(sph_body), DiffusionReactionSimpleData<BodyType, BaseParticlesType, BaseMaterialType>(sph_body),
	  drho_dt_(this->particles_->drho_dt_), diffusion_dt_prior_(this->particles_->diffusion_dt_prior_),
	  source_index_(source_index) {}
//=================================================================================================//

template <class BodyType, class BaseParticlesType, class BaseMaterialType>
void vdWAttractionHeatSource<BodyType, BaseParticlesType, BaseMaterialType>::Update(size_t index_i, Real dt)
{
	diffusion_dt_prior_[source_index_][index_i] -= this->material_->getAlpha()*drho_dt_[index_i]*2/3/k_B;
}

namespace SPH::fluid_dynamics{

template <class RiemannSolverType>
class vdWPressureRelaxationInner :
public BasePressureRelaxationInner<RiemannSolverType>
{
protected:
	StdLargeVec<Real> &temperature_;
public:
	explicit vdWPressureRelaxationInner(BaseBodyRelationInner &inner_relation) : 
			BasePressureRelaxationInner<RiemannSolverType>(inner_relation), 
			temperature_(*(std::get<indexScalar>(this->particles_->all_particle_data_)[this->particles_->all_variable_maps_[indexScalar]["Temperature"]])){};
	virtual ~vdWPressureRelaxationInner(){};
protected:
	virtual void Initialization(size_t index_i, Real dt = 0.0) override;
};

template <class RiemannSolverType>
void vdWPressureRelaxationInner<RiemannSolverType>::Initialization(size_t index_i, Real dt)
{
	this->rho_n_[index_i] += this->drho_dt_[index_i] * dt * 0.5;
	this->Vol_[index_i] = this->mass_[index_i] / this->rho_n_[index_i];
	this->p_[index_i] = this->material_->getPressure(this->rho_n_[index_i],temperature_[index_i]);
	this->pos_n_[index_i] += this->vel_n_[index_i] * dt * 0.5;
}

template <class BaseDensityRelaxationType>
class vdWDensityRelaxation : public BaseDensityRelaxationWithWall<BaseDensityRelaxationType>
{
	
protected:
	StdLargeVec<Real> &temperature_;
	virtual void Interaction(size_t index_i, Real dt = 0.0) override;
public:
	vdWDensityRelaxation(ComplexBodyRelation &fluid_wall_relation)
	: BaseDensityRelaxationWithWall<BaseDensityRelaxationType>(fluid_wall_relation),
	temperature_(*(std::get<indexScalar>(this->particles_->all_particle_data_)[this->particles_->all_variable_maps_[indexScalar]["Temperature"]])) 
	{}
	vdWDensityRelaxation(BaseBodyRelationInner &fluid_inner_relation,
									BaseBodyRelationContact &wall_contact_relation)
	:BaseDensityRelaxationWithWall<BaseDensityRelaxationType>(fluid_inner_relation,wall_contact_relation),
	temperature_(*(std::get<indexScalar>(this->particles_->all_particle_data_)[this->particles_->all_variable_maps_[indexScalar]["Temperature"]])) 
	{}
	vdWDensityRelaxation(ComplexBodyRelation &fluid_complex_relation,
									BaseBodyRelationContact &wall_contact_relation)
	:BaseDensityRelaxationWithWall<BaseDensityRelaxationType>(fluid_complex_relation,wall_contact_relation),
	temperature_(*(std::get<indexScalar>(this->particles_->all_particle_data_)[this->particles_->all_variable_maps_[indexScalar]["Temperature"]])) 
	{}
	virtual ~vdWDensityRelaxation(){};
};
//=================================================================================================//

template <class BaseDensityRelaxationType>
void vdWDensityRelaxation<BaseDensityRelaxationType>::Interaction(size_t index_i, Real dt)
{
	BaseDensityRelaxationType::Interaction(index_i, dt);

	FluidState state_i(this->rho_n_[index_i], this->vel_n_[index_i], this->p_[index_i]);
	Real density_change_rate = 0.0;
	for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
	{
		Vecd &dvel_dt_prior_i = this->dvel_dt_prior_[index_i];

		StdLargeVec<Real> &Vol_k = *(this->wall_Vol_[k]);
		StdLargeVec<Vecd> &vel_ave_k = *(this->wall_vel_ave_[k]);
		StdLargeVec<Vecd> &dvel_dt_ave_k = *(this->wall_dvel_dt_ave_[k]);
		StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
		Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
		for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
		{
			size_t index_j = wall_neighborhood.j_[n];
			Vecd &e_ij = wall_neighborhood.e_ij_[n];
			Real r_ij = wall_neighborhood.r_ij_[n];
			Real dW_ij = wall_neighborhood.dW_ij_[n];

			Real face_wall_external_acceleration = dot((dvel_dt_prior_i - dvel_dt_ave_k[index_j]), -e_ij);
			Vecd vel_in_wall = 2.0 * vel_ave_k[index_j] - state_i.vel_;
			Real dp = state_i.rho_ * r_ij * SMAX(0.0, face_wall_external_acceleration);
			Real p_in_wall = state_i.p_ + dp;
			Real rho_in_wall = state_i.rho_ + dp/this->material_->getSoundSpeed(state_i.p_ + 0.5*dp, temperature_[index_i]); //sound speed = dp/drho
			FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);
			Vecd vel_star = this->riemann_solver_.getVStar(state_i, state_j, n_k[index_j]);
			density_change_rate += 2.0 * state_i.rho_ * Vol_k[index_j] * dot(state_i.vel_ - vel_star, e_ij) * dW_ij;
		}
	}
	this->drho_dt_[index_i] += density_change_rate;
}
	using vdWPressureRelaxationRiemannWithWall = BasePressureRelaxationWithWall<PressureRelaxation<vdWPressureRelaxationInner<NoRiemannSolver>>>;
	using vdWDensityRelaxationRiemannWithWall = vdWDensityRelaxation<DensityRelaxationRiemannInner>;
	using vdWExtendMultiPhasePressureRelaxationRiemannWithWall = ExtendPressureRelaxationWithWall<ExtendPressureRelaxation<BasePressureRelaxationMultiPhase<vdWPressureRelaxationInner<NoRiemannSolver>>>>;
	using vdWMultiPhaseDensityRelaxationRiemannWithWall = vdWDensityRelaxation<MultiPhaseDensityRelaxationRiemann>;

}