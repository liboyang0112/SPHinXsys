#ifndef __CASE
#define __CASE
#include <cassert>
/**
 * @file 	case.h
 * @brief 	Numerical parameters and body defination for 2D two-phase dambreak flow.
 * @author 	Chi Zhang and Xiangyu Hu
 */
/**
 * @brief 	SPHinXsys Library.
 */
#include "multi_length.h"
#include "multi_length_contact.h"
/**
 * @brief Namespace cite here.
 */
using namespace SPH;
/**
 * @brief Basic geometry parameters and numerical setup.
 */
auto getParticleTemperature(BaseParticles *base_particles_){
	return std::get<indexScalar>(base_particles_->all_particle_data_)[base_particles_->all_variable_maps_[indexScalar]["Temperature"]];
}
int dim=2;
int dof = 2;
const Real k_B = 1;
#include "vdW_dynamics_multi_phase.hpp"
Real DL = 50;						   /**< Tank length. */
Real DH = 50;						   /**< Tank height. */
Real particle_spacing_ref = DL / 50.0; /**< Initial reference particle spacing. */

/** Material properties. */
Real bias_coff = 0.0;
Real alpha = Pi / 6.0;
Vec2d bias_direction(cos(alpha), sin(alpha));
/**
 *@brief Temperatures.
 */
//Real phi_upper_wall = 273.15;
//Real phi_lower_wall = 273.15;
//Real phi_side_wall = 273.15;
//Real phi_fluid_initial = 273.15;
//Real phi_gas_initial = 273.15;

Real phi_upper_wall = 0.2;
Real phi_lower_wall = 0.2;
Real phi_side_wall = 0.2;
Real phi_gas_initial = 0.2;

/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.5;							  // 1000;						  /**< Reference density of water. */
Real rho0_a = 1.100859077806;					  /**< Reference density of air. */
Real gravity_g = 0.01;								  /**< Gravity force of fluid. */
Real U_max = 1.0;								  /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;						  /**< Reference sound speed. */
Real mu_f = 2e3;								  /**< Water viscosity. */
Real mu_a = 5.0e-5;								  /**< Air visocsity. */
Real contact_angle = (150.0 / 180.0) * 3.1415926; /**< Contact angle with Wall. */
Real tension_force = 0.008;
/** vdW properties. */
Real gamma_vdw = (dof+2) / dof;	//(dim+2)/dim
Real rho0_water = 2;		// max density kg/m3
Real alpha_water = 2;		// attraction force
Real molmass_water = 18e-3; // kg/mol
Real molmass_air = 27e-3;
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
	// geometry
	std::vector<Vecd> water_block_shape;
	//solid bulk
	water_block_shape.push_back(Vecd(0.02 * DL, 0.02 * DH));
	water_block_shape.push_back(Vecd(0.02 * DL, 0.2 * DH));
	water_block_shape.push_back(Vecd(0.98 * DL, 0.2 * DH));
	water_block_shape.push_back(Vecd(0.98 * DL, 0.02 * DH));
	water_block_shape.push_back(Vecd(0.02 * DL, 0.02 * DH));
	//floating water
	//water_block_shape.push_back(Vecd(0.2 * DL, 0.2 * DH));
	//water_block_shape.push_back(Vecd(0.2 * DL, 0.5 * DH));
	//water_block_shape.push_back(Vecd(0.8 * DL, 0.5 * DH));
	//water_block_shape.push_back(Vecd(0.8 * DL, 0.2 * DH));
	//water_block_shape.push_back(Vecd(0.2 * DL, 0.2 * DH));
	return water_block_shape;
}
std::vector<Vecd> createAirBlockShape()
{
	// geometry
	std::vector<Vecd> water_block_shape;
	//solid bulk
	water_block_shape.push_back(Vecd(0.1 * DL, 0.3 * DH));
	water_block_shape.push_back(Vecd(0.1 * DL, 0.8 * DH));
	water_block_shape.push_back(Vecd(0.9 * DL, 0.8 * DH));
	water_block_shape.push_back(Vecd(0.9 * DL, 0.3 * DH));
	water_block_shape.push_back(Vecd(0.1 * DL, 0.3 * DH));
	//floating water
	//water_block_shape.push_back(Vecd(0.2 * DL, 0.2 * DH));
	//water_block_shape.push_back(Vecd(0.2 * DL, 0.5 * DH));
	//water_block_shape.push_back(Vecd(0.8 * DL, 0.5 * DH));
	//water_block_shape.push_back(Vecd(0.8 * DL, 0.2 * DH));
	//water_block_shape.push_back(Vecd(0.2 * DL, 0.2 * DH));
	return water_block_shape;
}
std::vector<Vecd> createPhotonBlockShape()
{
	// geometry
	std::vector<Vecd> photon_block_shape;
	photon_block_shape.push_back(Vecd(0.35 * DL, 0.56 * DH));
	photon_block_shape.push_back(Vecd(0.35 * DL, 0.6 * DH));
	photon_block_shape.push_back(Vecd(0.45 * DL, 0.6 * DH));
	photon_block_shape.push_back(Vecd(0.45 * DL, 0.56 * DH));
	photon_block_shape.push_back(Vecd(0.35 * DL, 0.56 * DH));
	return photon_block_shape;
}
/** create outer wall shape */
std::vector<Vecd> createOuterWallShape(Real BW)
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

class vdWFluid : public Fluid
{
public:
	explicit vdWFluid(Real rho0, Real rho_max, Real gamma, Real alpha, Real molmass, Real mu = 0.0)
		// rho0 is the density maximum
		: Fluid(rho0, mu), alpha_(alpha), molmass_(molmass), rho_max_(rho_max),  gamma_(gamma)
	{
		material_type_ = "vdWFluid";
	};
	Real alpha_, molmass_, rho_max_, gamma_; // alpha is the attraction force
	virtual ~vdWFluid(){};
	Real getAlpha() { return alpha_; }
	Real HeatCapacityRatio() { return gamma_; };
	virtual Real DensityFromPT(Real p, Real T){};
	virtual Real getPressure(Real rho) override {return 0;};
	virtual Real DensityFromPressure(Real p) override {return 0;};
	virtual Real getPressure(Real rho, Real T) override;
	virtual Real getViscosity(Real rho, Real T) override;
	virtual Real getSoundSpeed(Real p, Real rho) override;
	virtual vdWFluid *ThisObjectPtr() override { return this; };
};

Real vdWFluid::getViscosity(Real rho, Real T){
	return mu_*rho/rho0_;  //constant kinematic viscosity
}
Real vdWFluid::getPressure(Real rho, Real T)
{
	Real ratio = 1. / (1 - rho / rho_max_);
	if(ratio<0 || ratio!=ratio) ratio=100;
	//Real ret = rho * k_B * T * ratio;
	Real ret = rho * k_B * T * ratio - alpha_ * rho * rho;
	assert (ret == ret);
	return ret;
}

Real vdWFluid::getSoundSpeed(Real p, Real rho)
{
	Real ratio = 1. / (1 - rho / rho_max_);
	Real ret = sqrt(fabs(gamma_ * ratio * (alpha_ * rho + p / rho) - 2 * alpha_ * rho));
	if (ret != ret)
	{
		printf("WARNING: Sound Speed is NAN");
		exit(0);
	}
	return ret;
}
using vdWParticles = DiffusionReactionParticles<FluidParticles, vdWFluid>;
using WallParticles = DiffusionReactionParticles<SolidParticles, Solid>;
/**
 * application dependent solid body initial condition
 */
class ThermosolidBodyInitialCondition
	: public DiffusionReactionInitialCondition<SolidBody, SolidParticles, Solid>
{
protected:
	size_t phi_;
	Real T0_;
	Real BW;
	void Update(size_t index_i, Real dt) override
	{

		if (-BW <= pos_n_[index_i][1] && pos_n_[index_i][1] <= 0.0)
		{
			species_n_[phi_][index_i] = T0_;
		}
		else if (DH <= pos_n_[index_i][1] && pos_n_[index_i][1] <= DH + BW)
		{
			species_n_[phi_][index_i] = phi_upper_wall;
		}
		else
		{
			species_n_[phi_][index_i] = phi_side_wall;
		}
	};

public:
	ThermosolidBodyInitialCondition(SolidBody &diffusion_solid_body, Real T0, Real BW_)
		: DiffusionReactionInitialCondition<SolidBody, SolidParticles, Solid>(diffusion_solid_body),
		T0_(T0)
	{
		phi_ = material_->SpeciesIndexMap()["Temperature"];
		BW = BW_;
	};
};

/**
 * application dependent fluid body initial condition
 */
class ThermofluidBodyInitialCondition
	: public DiffusionReactionInitialCondition<FluidBody, FluidParticles, vdWFluid>
{
protected:
	size_t phi_;
	Real T0_;
	void Update(size_t index_i, Real dt) override
	{

		if (0 <= pos_n_[index_i][1] && pos_n_[index_i][1] <= DH)
		{
			species_n_[phi_][index_i] = T0_;
		}
	};

public:
	ThermofluidBodyInitialCondition(FluidBody &diffusion_fluid_body, Real T0)
		: DiffusionReactionInitialCondition<FluidBody, FluidParticles, vdWFluid>(diffusion_fluid_body),
		T0_(T0)
	{
		phi_ = material_->SpeciesIndexMap()["Temperature"];
	};
};

/**
 * application dependent gas body initial condition
 */
class ThermogasBodyInitialCondition
	: public DiffusionReactionInitialCondition<FluidBody, FluidParticles, vdWFluid>
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
		: DiffusionReactionInitialCondition<FluidBody, FluidParticles, vdWFluid>(diffusion_fluid_body)
	{
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
		  FluidBody, FluidParticles, vdWFluid,
		  RelaxationOfAllDiffussionSpeciesComplex<
			  FluidBody, FluidParticles, vdWFluid, SolidBody, SolidParticles, Solid>,
		  ComplexBodyRelation>
{
public:
	explicit ThermalRelaxationComplex(ComplexBodyRelation &body_complex_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_complex_relation){};
	virtual ~ThermalRelaxationComplex(){};
};
class ThermalRelaxationComplexWall
	: public RelaxationOfAllDiffusionSpeciesRK2<
		  SolidBody, SolidParticles, Solid,
		  RelaxationOfAllDiffussionSpeciesComplex<
			  SolidBody, SolidParticles, Solid, FluidBody, FluidParticles, vdWFluid>,
		  ComplexBodyRelation>
{
public:
	explicit ThermalRelaxationComplexWall(ComplexBodyRelation &body_complex_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_complex_relation){};
	virtual ~ThermalRelaxationComplexWall(){};
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
	WaterBlock(SPHSystem &sph_system, const string &body_name, shared_ptr<SPHAdaptation> adp)
		: FluidBody(sph_system, body_name, adp)
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
class PhotonBlock : public FluidBody
{
public:
	PhotonBlock(SPHSystem &sph_system, const string &body_name, shared_ptr<SPHAdaptation> adp)
		: FluidBody(sph_system, body_name, adp)
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createPhotonBlockShape(), ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
/**
 *@brief 	Air body definition.
 */
class AirBlock : public FluidBody
{
public:
        AirBlock(SPHSystem &sph_system, const string &body_name, shared_ptr<SPHAdaptation> adp)
                : FluidBody(sph_system, body_name, adp)
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createAirBlockShape(), ShapeBooleanOps::add);
		//multi_polygon.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};

/**
 * @brief 	Wall boundary body definition.
 */

class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &sph_system, string body_name, Real BW, shared_ptr<SPHAdaptation> adp)
		: SolidBody(sph_system, body_name, adp) //, makeShared<SPHAdaptation>(1.3, 1))
	{
		/** Geomtry definition. */
		std::vector<Vecd> outer_shape = createOuterWallShape(BW);
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
class FluidMaterial
	: public DiffusionReaction<FluidParticles, vdWFluid>
{
public:
	FluidMaterial(Real diffusion_coff)
		: DiffusionReaction<FluidParticles, vdWFluid>({"Temperature"},
																  rho0_f, rho0_water, gamma_vdw, alpha_water, molmass_water, mu_f)
	{
		initializeAnDiffusion<DirectionalDiffusion>("Temperature", "Temperature",
													diffusion_coff, bias_coff, bias_direction);
		initializeASource("Temperature");
	};
	FluidMaterial(Real rho0, Real rho0_m, Real gamma, Real alpha, Real molmass, Real mu, Real diffusion_coff)
		: DiffusionReaction<FluidParticles, vdWFluid>({"Temperature"},
																  rho0, rho0_m, gamma, alpha, molmass, mu)
	{
		initializeAnDiffusion<DirectionalDiffusion>("Temperature", "Temperature",
													diffusion_coff, bias_coff, bias_direction);
		initializeASource("Temperature");
	};
};

class WallMaterial
	: public DiffusionReaction<SolidParticles, Solid>
{
public:
	WallMaterial(Real diffusion_coff, Real rho0, Real stiffness)
		: DiffusionReaction<SolidParticles, Solid>({"Temperature"}, rho0, stiffness)
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
	Real InducedHeating(Vecd &position) {
		return 0;
		return position[0]<20?2:0; 
	}
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
	StdVec<StdLargeVec<Real>> &diffusion_dt_prior_;
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
class StressTensorHeatSource // The Heat from viscous force and pressure, only inner forces are included.
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
	Real smoothing_length_;
	StdLargeVec<Real> &Vol_, &rho_n_, &p_;
	StdVec<StdLargeVec<Real>> &diffusion_dt_prior_;
	StdLargeVec<Vecd> &vel_n_;
	StdLargeVec<Real> &temperature_;
	virtual void Interaction(size_t index_i, Real dt = 0.0) override;
};

template <class BodyType, class BaseParticlesType, class BaseMaterialType>
StressTensorHeatSource<BodyType, BaseParticlesType, BaseMaterialType>::StressTensorHeatSource(BaseBodyRelationInner &inner_relation, size_t source_index)
	: InteractionDynamics(*inner_relation.sph_body_),
	  DiffusionReactionInnerData<BodyType, BaseParticlesType, BaseMaterialType>(inner_relation),
	  Vol_(this->particles_->Vol_), rho_n_(this->particles_->rho_n_), p_(this->particles_->p_),
	  vel_n_(this->particles_->vel_n_),
	  smoothing_length_(sph_adaptation_->ReferenceSmoothingLength()),
	  diffusion_dt_prior_(this->particles_->diffusion_dt_prior_),
	  source_index_(source_index),
	  riemann_solver_(*(this->material_), *(this->material_)),
	  temperature_(*getParticleTemperature(this->particles_)){}

//=================================================================================================//

template <class BodyType, class BaseParticlesType, class BaseMaterialType>
void StressTensorHeatSource<BodyType, BaseParticlesType, BaseMaterialType>::Interaction(size_t index_i, Real dt)
{
	Real rho_i = rho_n_[index_i];
	Real internalEnergyIncrease(0);
	const Vecd &vel_i = vel_n_[index_i];
	Real pi_att=-((vdWFluid*)this->material_)->getAlpha()*rho_n_[index_i]*rho_n_[index_i];
	Real pi = p_[index_i]-pi_att;
	FluidState state_i(rho_n_[index_i], vel_n_[index_i], pi);
//	FluidState state_i_att(rho_n_[index_i], vel_n_[index_i], pi_att);
	Vecd vel_derivative(0);
	Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
	for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
	{
		size_t index_j = inner_neighborhood.j_[n];
		Real dW_ij = inner_neighborhood.dW_ij_[n];
		Real dW_ij_n = inner_neighborhood.dW_ij_n_[0][n];
		Vecd &e_ij = inner_neighborhood.e_ij_[n];
		Vecd vij = vel_i - vel_n_[index_j];
		Real pj_att=-((vdWFluid*)this->material_)->getAlpha()*rho_n_[index_j]*rho_n_[index_j];
		Real pj = p_[index_j]-pj_att;
		FluidState state_j(rho_n_[index_j], vel_n_[index_j], pj);
//		FluidState state_j_att(rho_n_[index_j], vel_n_[index_j], pj_att);
		Real p_star = riemann_solver_.getPStar(state_i, state_j, e_ij);
//		Real p_star_att = riemann_solver_.getPStar(state_i_att, state_j_att, e_ij);
		// pressure
		internalEnergyIncrease += dot((p_star * dW_ij) * Vol_[index_j] * e_ij / rho_i, vij);
//		internalEnergyIncrease -= 2.0 * ((vdWFluid*)this->material_)->getAlpha() * Vol_k[index_j] / state_j.rho_ * dot(state_i.vel_ - state_j.vel_,e_ij) * dW_ij_att;
		// viscous
		vel_derivative = vij / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
		Real visheat = -dot(this->material_->getViscosity(this->rho_n_[index_i], temperature_[index_i]) * vel_derivative * Vol_[index_j] * dW_ij / rho_i, vij);
		internalEnergyIncrease += visheat;
	}
	diffusion_dt_prior_[source_index_][index_i] += internalEnergyIncrease * 2 / dof / k_B;
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
	StdVec<StdLargeVec<Real>> &diffusion_dt_prior_;
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
	Real ret = -this->material_->getAlpha() * drho_dt_[index_i] * 2 / dof / k_B;
	diffusion_dt_prior_[source_index_][index_i] += ret;
}

namespace SPH::fluid_dynamics
{

	template <class RiemannSolverType>
	class vdWPressureRelaxationInner : public BasePressureRelaxationInner<RiemannSolverType>
	{
	protected:
		StdLargeVec<Real> &temperature_;
		size_t source_index_;
		Real smoothing_length_;
		StdVec<StdLargeVec<Real>> &diffusion_dt_prior_;
	public:
		explicit vdWPressureRelaxationInner(BaseBodyRelationInner &inner_relation) : 
		BasePressureRelaxationInner<RiemannSolverType>(inner_relation),
	  	smoothing_length_(this->sph_adaptation_->ReferenceSmoothingLength()),
		diffusion_dt_prior_(((vdWParticles*)(this->particles_))->diffusion_dt_prior_),
		source_index_(0),
		temperature_(*getParticleTemperature(this->particles_)){};
		virtual ~vdWPressureRelaxationInner(){};

	protected:
		virtual void Initialization(size_t index_i, Real dt = 0.0) override;
		virtual void Interaction(size_t index_i, Real dt = 0.0) override;
	};
	template<class RiemannSolverType>
    void vdWPressureRelaxationInner<RiemannSolverType>::Interaction(size_t index_i, Real dt)
	{
		Real pi_att=-((vdWFluid*)this->material_)->getAlpha()*this->rho_n_[index_i]*this->rho_n_[index_i];
		Real pi_perf = this->rho_n_[index_i] * k_B * temperature_[index_i];
		Real pi_rep = this->p_[index_i]-pi_att-pi_perf;
		FluidState state_i(this->rho_n_[index_i], this->vel_n_[index_i], pi_perf);
		FluidState state_i_r(this->rho_n_[index_i], this->vel_n_[index_i], pi_rep);
		Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
		this->dvel_dt_[index_i] = this->dvel_dt_prior_[index_i];
		Real internalEnergyIncrease(0);
		Vecd acceleration_i;
		Vecd vel_derivative(0);
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Real dW_ij = inner_neighborhood.dW_ij_[n];
			Real dW_ij_n = inner_neighborhood.dW_ij_n_[0][n];
			Vecd& e_ij = inner_neighborhood.e_ij_[n];
			Real r_ij = inner_neighborhood.r_ij_[n];
			Real pj_att = -((vdWFluid*)this->material_)->getAlpha()*this->rho_n_[index_j]*this->rho_n_[index_j];
			Real pj_perf = this->rho_n_[index_j] * k_B * temperature_[index_j];
			Real pj_rep = this->p_[index_j] - pj_att-pj_perf;
			Vecd vij = this->vel_n_[index_i] - this->vel_n_[index_j];
			//repulsion H heat
			FluidState state_r(this->rho_n_[index_j], this->vel_n_[index_j], pj_rep);
			Real p_star = this->riemann_solver_.getPStar(state_i_r, state_r, e_ij);
			acceleration_i = -2.0 * (p_star * dW_ij) * this->Vol_[index_j] * e_ij / state_i_r.rho_;
			internalEnergyIncrease -= dot(acceleration_i, vij)*state_i_r.p_/(state_i_r.p_+pj_rep);
			this->dvel_dt_[index_i] += acceleration_i;
			//perfect gas 2H heat
			FluidState state_j(this->rho_n_[index_j], this->vel_n_[index_j], pj_perf);
			p_star = this->riemann_solver_.getPStar(state_i, state_j, e_ij);
			acceleration_i = -2.0 * (p_star * dW_ij_n) * this->Vol_[index_j] * e_ij / state_i.rho_;
			internalEnergyIncrease -= dot(acceleration_i, vij)*state_i.p_/(state_i.p_+pj_perf);
			//att 2H no heat
			this->dvel_dt_[index_i] += acceleration_i + 2.0 * ((vdWFluid*)this->material_)->getAlpha() * this->Vol_[index_j] * state_j.rho_ * e_ij * dW_ij_n;
			//visc
			vel_derivative = vij / (r_ij + 0.01 * smoothing_length_);
			acceleration_i = 2.0 * this->material_->getViscosity(this->rho_n_[index_i], temperature_[index_i]) * vel_derivative * this->Vol_[index_j] * dW_ij_n / state_i.rho_;
			internalEnergyIncrease -= 0.5*dot(acceleration_i, vij);
			this->dvel_dt_[index_i] += acceleration_i;
			//Core force:
			//10.1016/j.ijheatmasstransfer.2018.10.119
			//https://www.docin.com/p-1406155771.html
			//W.G. Hoover, Smooth Particle Applied Mechanics: The State of the Art (Advanced Series in Nonlinear Dynamics), World Scientific Publishing Company, River Edge, NJ, 2006.
			Real sigma2 = pow(this->sph_adaptation_->ReferenceSpacing(),2);
			Real r_ij2 = r_ij*r_ij;
			if(r_ij2 < sigma2)
			this->dvel_dt_[index_i] += this->material_->ReferenceDensity()*e_ij*(3*r_ij/sigma2)*pow(1-r_ij2/sigma2,3) / state_i.rho_ / this->Vol_[index_i];
		}
		temperature_[index_i] += internalEnergyIncrease * 2 / dof / k_B * dt;
	}    

	template <class RiemannSolverType>
	void vdWPressureRelaxationInner<RiemannSolverType>::Initialization(size_t index_i, Real dt)
	{
		//this->rho_n_[index_i] += this->drho_dt_[index_i] * dt * 0.5;
		this->Vol_[index_i] = this->mass_[index_i] / this->rho_n_[index_i];
		this->p_[index_i] = this->material_->getPressure(this->rho_n_[index_i], temperature_[index_i]);
		this->pos_n_[index_i] += this->vel_n_[index_i] * dt * 0.5;
	}

	template <class BasePressureRelaxationType>
	class vdWPressureRelaxation : public RelaxationWithWall<BasePressureRelaxationType>
	{
	public:
		// template for different combination of constructing body relations
		template <class BaseBodyRelationType>
		vdWPressureRelaxation(BaseBodyRelationType &base_body_relation,
						   BaseBodyRelationContact &wall_contact_relation);
		virtual ~vdWPressureRelaxation(){};
	protected:
		virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		virtual Vecd computeNonConservativeAcceleration(size_t index_i) override;
	};
	template <class BasePressureRelaxationType>
	template <class BaseBodyRelationType>
	vdWPressureRelaxation<BasePressureRelaxationType>::
		vdWPressureRelaxation(BaseBodyRelationType &base_body_relation,
						   BaseBodyRelationContact &wall_contact_relation)
		: RelaxationWithWall<BasePressureRelaxationType>(base_body_relation, wall_contact_relation) {}
	//=================================================================================================//
	template <class BasePressureRelaxationType>
	void vdWPressureRelaxation<BasePressureRelaxationType>::Interaction(size_t index_i, Real dt)
	{
		if(index_i == 57){
			if(this->temperature_[index_i]>0.5){
				int halt = 1;
			}
		}
		BasePressureRelaxationType::Interaction(index_i, dt);
		Real pi_att = -((vdWFluid*)this->material_)->getAlpha()*this->rho_n_[index_i]*this->rho_n_[index_i];
		Real pi = this->p_[index_i]-pi_att;
		FluidState state_i(this->rho_n_[index_i], this->vel_n_[index_i], pi);
		FluidState state_i_att(this->rho_n_[index_i], this->vel_n_[index_i], pi_att);
		Vecd dvel_dt_prior_i = computeNonConservativeAcceleration(index_i);
		Vecd acceleration(0.0);
		Vecd acceleration_surface(0);
		Vecd acceleration_i;
		Real internalEnergyIncrease = 0;
		Vecd vel_derivative(0);
		for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
		{
			StdLargeVec<Real> &Vol_k = *(this->wall_Vol_[k]);
			StdLargeVec<Vecd> &vel_ave_k = *(this->wall_vel_ave_[k]);
			StdLargeVec<Vecd> &dvel_dt_ave_k = *(this->wall_dvel_dt_ave_[k]);
			StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
			Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
			//Real rho_in_wall = 1.6;// + dp / this->material_->getSoundSpeed(state_i.p_ + 0.5 * dp, state_i.rho_); // sound speed = dp/drho
			//Real pj_att = ((vdWFluid*)this->material_)->getAlpha()*rho_in_wall*rho_in_wall;
			Real p_in_wall = 1.6;
			Vecd rhograd(0);
			Real densum = 0;
			for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
			{
				size_t index_j = wall_neighborhood.j_[n];
				Vecd &e_ij = wall_neighborhood.e_ij_[n];
				Real dW_ij = wall_neighborhood.dW_ij_[n];
				//Real dW_ij_att = 0;
				Real dW_ij_att = wall_neighborhood.dW_ij_n_[0][n];
				Real r_ij = wall_neighborhood.r_ij_[n];
				densum += wall_neighborhood.W_ij_[n];
				rhograd += dW_ij*e_ij;
				Real face_wall_external_acceleration = dot((dvel_dt_prior_i - dvel_dt_ave_k[index_j]), -e_ij);
				Vecd vel_in_wall = 2.0 * vel_ave_k[index_j] - state_i.vel_;
				Vecd vij = state_i.vel_ - vel_in_wall;
				Real dp = state_i.rho_ * r_ij * SMAX(0.0, face_wall_external_acceleration);
				Real rho_in_wall = state_i.rho_;// + dp / this->material_->getSoundSpeed(state_i.p_ + 0.5 * dp, state_i.rho_); // sound speed = dp/drho
				//Real rho_in_wall = 1.6;// + dp / this->material_->getSoundSpeed(state_i.p_ + 0.5 * dp, state_i.rho_); // sound speed = dp/drho
				Real pj_att = -((vdWFluid*)this->material_)->getAlpha()*rho_in_wall*rho_in_wall;
				//Real p_in_wall = state_i.p_-pj_att;
				FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);
				//FluidState state_j_att(this->rho_n_[index_j], this->vel_n_[index_j], pj_att);
				//Real p_star = this->riemann_solver_.getPStar(state_i, state_j, n_k[index_j]);
				//Real p_star_att = this->riemann_solver_.getPStar(state_i_att, state_j_att, e_ij);
				vel_derivative = vij / (r_ij + 0.01 * this->smoothing_length_);
				//acceleration_i = -2.0 * (p_star * dW_ij) * e_ij * Vol_k[index_j] / state_i.rho_
				acceleration_i = //-5.0 * dW_ij * e_ij+ 
				2.0 * this->material_->getViscosity(this->rho_n_[index_i], this->temperature_[index_i]) * vel_derivative * Vol_k[index_j] * dW_ij / state_i.rho_;
				internalEnergyIncrease -= 0.5*dot(acceleration_i, vij);
				acceleration += acceleration_i;
				//acceleration_surface += e_ij * dW_ij;
				Real sigma2 = pow(this->sph_adaptation_->ReferenceSpacing(),2);
				Real r_ij2 = r_ij*r_ij;
				if(r_ij2 < sigma2)
				this->dvel_dt_[index_i] += e_ij*(3*r_ij/sigma2)*pow(1-r_ij2/sigma2,3) / state_i.rho_ / this->Vol_[index_i];
			}
			if(densum >= 0.2*FluidWallData::contact_material_[k]->ReferenceDensity()){
				if(dot(rhograd,state_i.vel_)>0){
					state_i.vel_ -= rhograd*2.*dot(rhograd,state_i.vel_)/pow(rhograd.norm(),2);
					acceleration -= rhograd*2.*dot(rhograd,acceleration)/pow(rhograd.norm(),2);
					acceleration_surface = Vecd(0);
					this->dvel_dt_[index_i] = 0;
				}
			}
		}
		this->temperature_[index_i] += internalEnergyIncrease * 2 / dof / k_B * dt;
		this->dvel_dt_[index_i] += acceleration + acceleration_surface;
	}
	//=================================================================================================//
	template <class BasePressureRelaxationType>
	Vecd vdWPressureRelaxation<BasePressureRelaxationType>::computeNonConservativeAcceleration(size_t index_i)
	{
		return this->dvel_dt_prior_[index_i];
	}

	template <class BaseDensityRelaxationType>
	class vdWDensityRelaxation : public BaseDensityRelaxationWithWall<BaseDensityRelaxationType>
	{

	protected:
		StdLargeVec<Real> &temperature_;
		virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		//	virtual void Initialization(size_t index_i, Real dt = 0.0) override;
	public:
		vdWDensityRelaxation(ComplexBodyRelation &fluid_wall_relation)
			: BaseDensityRelaxationWithWall<BaseDensityRelaxationType>(fluid_wall_relation),
			  temperature_(*getParticleTemperature(this->particles_))
		{
		}
		vdWDensityRelaxation(BaseBodyRelationInner &fluid_inner_relation,
							 BaseBodyRelationContact &wall_contact_relation)
			: BaseDensityRelaxationWithWall<BaseDensityRelaxationType>(fluid_inner_relation, wall_contact_relation),
			  temperature_(*getParticleTemperature(this->particles_))
		{
		}
		vdWDensityRelaxation(ComplexBodyRelation &fluid_complex_relation,
							 BaseBodyRelationContact &wall_contact_relation)
			: BaseDensityRelaxationWithWall<BaseDensityRelaxationType>(fluid_complex_relation, wall_contact_relation),
			  temperature_(*getParticleTemperature(this->particles_))
		{
		}
		virtual ~vdWDensityRelaxation(){};
	};
	//=================================================================================================//
	template <class RiemannSolverType>
	class vdWBaseDensityRelaxationInner : public BaseDensityRelaxationInner<RiemannSolverType>
	{
		public:
		vdWBaseDensityRelaxationInner(BaseBodyRelationInner &inner_relation) : BaseDensityRelaxationInner<RiemannSolverType> (inner_relation){};
		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
	};
 	template<class RiemannSolverType>
    void vdWBaseDensityRelaxationInner<RiemannSolverType>::Interaction(size_t index_i, Real dt)
	{
		FluidState state_i(this->rho_n_[index_i], this->vel_n_[index_i], this->p_[index_i]);
		Real density_change_rate = 0.0;
		Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			Vecd& e_ij = inner_neighborhood.e_ij_[n];
			Real dW_ij = inner_neighborhood.dW_ij_n_[0][n];
			FluidState state_j(this->rho_n_[index_j], this->vel_n_[index_j], this->p_[index_j]);
			Vecd vel_star = this->riemann_solver_.getVStar(state_i, state_j, e_ij);
			density_change_rate += 2.0 * state_i.rho_ * this->Vol_[index_j] * dot(state_i.vel_ - vel_star, e_ij) * dW_ij;
		}
		this->drho_dt_[index_i] = density_change_rate;
	};
	// template <class BaseDensityRelaxationType>
	// void vdWDensityRelaxation<BaseDensityRelaxationType>::Initialization(size_t index_i, Real dt)
	//{
	// }
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
				Real dW_ij = wall_neighborhood.dW_ij_n_[0][n];

				Real face_wall_external_acceleration = dot((dvel_dt_prior_i - dvel_dt_ave_k[index_j]), -e_ij);
				Vecd vel_in_wall = 2.0 * vel_ave_k[index_j] - state_i.vel_;
				Real dp = state_i.rho_ * r_ij * SMAX(0.0, face_wall_external_acceleration);
				Real p_in_wall = state_i.p_ + dp;
				Real rho_in_wall = state_i.rho_ + dp / this->material_->getSoundSpeed(state_i.p_ + 0.5 * dp, state_i.rho_); // sound speed = dp/drho
				FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);
				Vecd vel_star = this->riemann_solver_.getVStar(state_i, state_j, n_k[index_j]);
				density_change_rate += 2.0 * state_i.rho_ * Vol_k[index_j] * dot(state_i.vel_ - vel_star, e_ij) * dW_ij;
				assert (density_change_rate == density_change_rate);
			}
		}
		this->drho_dt_[index_i] += density_change_rate;
	}


	using vdWPressureRelaxationRiemannWithWall = BasePressureRelaxationWithWall<vdWPressureRelaxation<vdWPressureRelaxationMultiPhase<vdWPressureRelaxationInner<NoRiemannSolver>>>>;
	using vdWDensityRelaxationRiemannWithWall = vdWDensityRelaxation<vdWBaseDensityRelaxationInner<NoRiemannSolver>>;
	using vdWExtendMultiPhasePressureRelaxationRiemannWithWall = ExtendPressureRelaxationWithWall<ExtendPressureRelaxation<BasePressureRelaxationMultiPhase<vdWPressureRelaxationInner<NoRiemannSolver>>>>;
	using vdWMultiPhaseDensityRelaxationRiemannWithWall = vdWDensityRelaxation<BaseDensityRelaxationMultiPhase<BaseDensityRelaxationInner<NoRiemannSolver>>>;

	class DensitySummationFreeSurfaceComplexWithoutUpdate : public DensitySummationFreeSurfaceComplex
	{
	public:
		DensitySummationFreeSurfaceComplexWithoutUpdate(BaseBodyRelationInner &inner_relation, BaseBodyRelationContact &contact_relation) : DensitySummationFreeSurfaceComplex(inner_relation, contact_relation){};
		DensitySummationFreeSurfaceComplexWithoutUpdate(ComplexBodyRelation &complex_relation, BaseBodyRelationContact &extra_contact_relation) : DensitySummationFreeSurfaceComplex(complex_relation, extra_contact_relation){};
		virtual ~DensitySummationFreeSurfaceComplexWithoutUpdate(){};

	protected:
		// virtual void Update(size_t index_i, Real dt = 0.0) override;
		virtual Real ReinitializedDensity(Real rho_sum, Real rho_0, Real rho_n) override { return rho_sum; };
		virtual void Interaction(size_t index_i, Real dt) override;
	};

	// void DensitySummationFreeSurfaceComplexWithoutUpdate::Update(size_t index_i, Real dt){}

	void DensitySummationFreeSurfaceComplexWithoutUpdate::Interaction(size_t index_i, Real dt)
	{
		Real sigma = W0_ * mass_[index_i];
		const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			sigma += inner_neighborhood.W_ij_n_[0][n] * mass_[inner_neighborhood.j_[n]];
		/** Contact interaction. */
		//Real inv_Vol_0_i = rho0_ / mass_[index_i];
		for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
		{
			StdLargeVec<Real> &contact_mass_k = *(this->contact_mass_[k]);
			Real contact_inv_rho0_k = contact_inv_rho0_[k];
			Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
			//	sigma += contact_neighborhood.W_ij_n_[0][n] * inv_Vol_0_i * contact_inv_rho0_k * contact_mass_k[contact_neighborhood.j_[n]];
				sigma += contact_neighborhood.W_ij_n_[0][n] * mass_[index_i];
			}
		}
		//Real rho_new = sigma * rho0_ * inv_sigma0_;
		Real rho_new = sigma;
		this->particles_->drho_dt_[index_i] = (rho_new-this->rho_sum_[index_i])/dt;
		this->rho_sum_[index_i] = rho_new;
	}
	
	template <int DataTypeIndex, typename VariableType>
	struct averageAParticleDataValue
	{
		void operator()(ParticleData &particle_data, size_t this_index, Real this_weight, size_t another_index, Real another_weight) const {
        	for (size_t i = 0; i != std::get<DataTypeIndex>(particle_data).size(); ++i){
        	    (*std::get<DataTypeIndex>(particle_data)[i])[this_index] =
					(*std::get<DataTypeIndex>(particle_data)[i])[this_index]*this_weight+
        	    	(*std::get<DataTypeIndex>(particle_data)[i])[another_index]*another_weight;
			}
		};
	};	
	template <>
	struct averageAParticleDataValue<indexPointer, void*>
	{
		void operator()(ParticleData &particle_data, size_t this_index, Real this_weight, size_t another_index, Real another_weight) const {
        	for (size_t i = 0; i != std::get<indexPointer>(particle_data).size(); ++i){
				if(this_weight < another_weight)
        	    	(*std::get<indexPointer>(particle_data)[i])[this_index] =
        	    		(*std::get<indexPointer>(particle_data)[i])[another_index];
			}
		};
	};

	class DensitySummationWithMergingAndSplitting : public DensitySummationFreeSurfaceComplex
	{
		ParticleDataOperation<averageAParticleDataValue> average_a_particle_value_;
		size_t n_particles_to_merge_;
		std::atomic_size_t n_particles_to_split_;
		ParticleData &all_particle_data_;
		StdLargeVec<size_t> &merge_chain_; // V[index_i] = index_j, meaning merge i to j
		StdLargeVec<size_t> &merge_list_; // list of indices of particles to merge
		StdLargeVec<size_t> &split_list_; // list of indices of particles to split
		StdLargeVec<size_t> &split_to_;  // list of indices of split result particles
		int splitN;
		Real mergeThreshold = 0.2;
		Real splitThreshold = 0.8;
	public:
		DensitySummationWithMergingAndSplitting(BaseBodyRelationInner &inner_relation,
		BaseBodyRelationContact &contact_relation,
		StdLargeVec<size_t> &cache1, StdLargeVec<size_t> &cache2,
		StdLargeVec<size_t> &cache3, StdLargeVec<size_t> &cache4) : 
		DensitySummationFreeSurfaceComplex(inner_relation, contact_relation),
		merge_chain_(cache1), merge_list_(cache2), split_list_(cache3), split_to_(cache4),
		all_particle_data_(base_particles_->all_particle_data_){
			splitN = pow(2,Vecd(0).size());
		};

		DensitySummationWithMergingAndSplitting(ComplexBodyRelation &complex_relation, BaseBodyRelationContact &extra_contact_relation,
		StdLargeVec<size_t> &cache1, StdLargeVec<size_t> &cache2,
		StdLargeVec<size_t> &cache3, StdLargeVec<size_t> &cache4) : 
		DensitySummationFreeSurfaceComplex(complex_relation, extra_contact_relation),
		merge_chain_(cache1), merge_list_(cache2), split_list_(cache3), split_to_(cache4),
		all_particle_data_(base_particles_->all_particle_data_){
			splitN = pow(2,Vecd(0).size());
		};
		virtual ~DensitySummationWithMergingAndSplitting(){};

	protected:
		// virtual void Update(size_t index_i, Real dt = 0.0) override;
		virtual Real ReinitializedDensity(Real rho_sum, Real rho_0, Real rho_n) override { return rho_sum; };
		virtual void Interaction(size_t index_i, Real dt) override;
		void createMergeList(){
			n_particles_to_merge_ = 0;
			for(size_t idx = base_particles_->total_real_particles_-1 ; idx >= 0 ; idx--){
				size_t idxj = merge_chain_[idx];
				if(idxj != idx){
					merge_chain_[idxj] = idxj;
					merge_list_[n_particles_to_merge_++] = idx;
				}
			}
		}
		void mergeParticles(){
			parallel_for(
				blocked_range<size_t>(0, n_particles_to_merge_),
				[&](const blocked_range<size_t> &r)
				{
					for (size_t i = r.begin(); i != r.end(); ++i)
					{
						size_t &to_merge = merge_list_[i];
						size_t &target = merge_chain_[to_merge];
						average_a_particle_value_(all_particle_data_,target, mass_[target], to_merge, mass_[to_merge]);
					}
				},
				ap
			);
		}
		void splitParticles(){
			size_t reservedSpace = splitN * n_particles_to_split_;
			for(size_t i = 0; i < reservedSpace; i++){
				if(i < n_particles_to_merge_){
					split_to_[i] = merge_list_[i];
				}else{
					split_to_[i] = i-n_particles_to_merge_+base_particles_->total_real_particles_;
				}
			}
			parallel_for(
				blocked_range<size_t>(0, n_particles_to_split_),
				[&](const blocked_range<size_t> &r)
				{
					for (size_t i = r.begin(); i != r.end(); ++i)
					{
						int spliti = i*splitN;
						mass_[split_list_[i]] /= splitN;
						for(size_t dim = 0; dim < Vecd(0).size(); dim++){
							for(size_t i = -1; i < 2 ; i+=2){
								base_particles_->copyFromAnotherParticle(split_to_[spliti++],split_list_[i]);
								base_particles_->pos_n_[split_to_[spliti]][dim] += i*sph_adaptation_->ReferenceSmoothingLength();
							}
						}
					}
				},
				ap
			);
		}
		void reset(){

		}
		void exec(Real dt) override{
			//Calculate density, create merge_chain_, split_list_
			DensitySummationFreeSurfaceComplex::exec();
			//Merging is more complex in parallel executing,
			//Need to decide which particle to merge first.
			//Translate merge_chain_ to merge_list_
			createMergeList();
			//executing particle merging,
			//Modify target particle, but do not kill merged particle.
			mergeParticles();
			//executing particle splitting, result particles take id's of merged particles.
			splitParticles();
			//kill the rest merged particles, we do it here 
			//because particle splitting can take the id of merged particles,
			//to reduce particle copying which takes relative more time.
			killParticles(*this->base_particles_, 
			n_particles_to_merge_ - splitN * n_particles_to_split_,
			merge_list_, split_to_, splitN * n_particles_to_split_);  
			reset();
		};
	};

	// void DensitySummationWithMergingAndSplitting::Update(size_t index_i, Real dt){}

	void DensitySummationWithMergingAndSplitting::Interaction(size_t index_i, Real dt)
	{
		Real sigma = W0_ * mass_[index_i];
		Real ratio = sigma;
		const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
		size_t closest_idx;
		Real closest_W = 0;
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n){
			sigma += inner_neighborhood.W_ij_n_[0][n] * mass_[inner_neighborhood.j_[n]];
			if(closest_W < inner_neighborhood.W_ij_n_[0][n]){
				closest_W = inner_neighborhood.W_ij_n_[0][n];
				closest_idx = inner_neighborhood.j_[n];
			}
		}
		ratio /= sigma;
		if(ratio < mergeThreshold){
			merge_chain_[index_i] = closest_idx;
		}else{
			merge_chain_[index_i] = index_i;
			if(ratio > splitThreshold){
				split_list_[n_particles_to_split_++] = index_i;
			}
		}
		/** Contact interaction. */
		//Real inv_Vol_0_i = rho0_ / mass_[index_i];
		for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
		{
			StdLargeVec<Real> &contact_mass_k = *(this->contact_mass_[k]);
			Real contact_inv_rho0_k = contact_inv_rho0_[k];
			Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
			//	sigma += contact_neighborhood.W_ij_n_[0][n] * inv_Vol_0_i * contact_inv_rho0_k * contact_mass_k[contact_neighborhood.j_[n]];
				sigma += contact_neighborhood.W_ij_n_[0][n] * mass_[index_i];
			}
		}
		//Real rho_new = sigma * rho0_ * inv_sigma0_;
		Real rho_new = sigma;
		this->particles_->drho_dt_[index_i] = (rho_new-this->rho_sum_[index_i])/dt;
		this->rho_sum_[index_i] = rho_new;
	}
	

	class vdWViscousAccelerationInner : public ViscousAccelerationInner
	{
	protected:
		StdLargeVec<Real> &temperature_;
		virtual void Interaction(size_t index_i, Real dt = 0.0) override;
	public:
		explicit vdWViscousAccelerationInner(BaseBodyRelationInner &inner_relation):ViscousAccelerationInner(inner_relation),
		temperature_(*getParticleTemperature(this->particles_)){};
		virtual ~vdWViscousAccelerationInner(){};
	};

	void vdWViscousAccelerationInner::Interaction(size_t index_i, Real dt)
	{
		Real rho_i = rho_n_[index_i];
		const Vecd &vel_i = vel_n_[index_i];
		Vecd acceleration(0), vel_derivative(0);
		const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			//viscous force
			vel_derivative = (vel_i - vel_n_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
			acceleration += 2.0 * this->material_->getViscosity(this->rho_n_[index_i], temperature_[index_i]) * vel_derivative * Vol_[index_j] * inner_neighborhood.dW_ij_[n] / rho_i;
		}
		dvel_dt_prior_[index_i] += acceleration;
	}


	// add particle data
	class KortewegTermCalc : public InteractionDynamics, public FluidDataInner
	{
	public:
		explicit KortewegTermCalc(BaseBodyRelationInner &inner_relation, Real weberNumber);
		virtual ~KortewegTermCalc(){};

	protected:
		Real weberNumber_;
		StdLargeVec<Real> &Vol_, &rho_n_, &p_;
		StdLargeVec<Vecd> &vel_n_, &dvel_dt_prior_;
		StdLargeVec<Matd> korteweg_tensor_;
		virtual void Interaction(size_t index_i, Real dt = 0.0) override;
	};
	KortewegTermCalc::KortewegTermCalc(BaseBodyRelationInner &inner_relation, Real weberNumber)
		: InteractionDynamics(*inner_relation.sph_body_),
		  FluidDataInner(inner_relation),
		  Vol_(particles_->Vol_), rho_n_(particles_->rho_n_), p_(particles_->p_),
		  vel_n_(particles_->vel_n_),
		  dvel_dt_prior_(particles_->dvel_dt_prior_),
		  weberNumber_(weberNumber)
	{
		// register Tensor to particle data
		particles_->registerAVariable<indexMatrix, Matd>(korteweg_tensor_, "Korteweg_tensor");
	}
	//=================================================================================================//
	void KortewegTermCalc::Interaction(size_t index_i, Real dt)
	{
		Real rho_i = rho_n_[index_i];
		const Vecd &vel_i = vel_n_[index_i];
		Vecd acceleration(0), vel_derivative(0);
		const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
		Vecd grad_rho(0, 0);
		Real laplace_rho(0);
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			grad_rho += rho_n_[index_j] * Vol_[index_j] * (rho_n_[index_j] / rho_i - 1) * inner_neighborhood.dW_ij_[n];
			laplace_rho += (4 - dim) * Vol_[index_j] * (rho_n_[index_j] - rho_i) / inner_neighborhood.r_ij_[n] * inner_neighborhood.dW_ij_[n];
		}
		Matd ret = 0.5 * (rho_i * laplace_rho + 0.5 * dot(grad_rho, grad_rho)) * Matd(1.0) - SimTK::outer(grad_rho, grad_rho);
		korteweg_tensor_[index_i] = weberNumber_ * ret;
	}

	class KortewegAccelerationInner : public InteractionDynamics, public FluidDataInner
	{
	public:
		explicit KortewegAccelerationInner(BaseBodyRelationInner &inner_relation);
		virtual ~KortewegAccelerationInner(){};

	protected:
		Real weberNumber_;
		StdLargeVec<Real> &Vol_, &rho_n_, &p_;
		StdLargeVec<Vecd> &vel_n_, &dvel_dt_prior_;
		StdLargeVec<Matd> &korteweg_tensor_;
		virtual void Interaction(size_t index_i, Real dt = 0.0) override;
	};
	KortewegAccelerationInner::KortewegAccelerationInner(BaseBodyRelationInner &inner_relation)
	: InteractionDynamics(*inner_relation.sph_body_),
	  FluidDataInner(inner_relation),
	  Vol_(particles_->Vol_), rho_n_(particles_->rho_n_), p_(particles_->p_),
	  vel_n_(particles_->vel_n_),
	  dvel_dt_prior_(particles_->dvel_dt_prior_),
	  korteweg_tensor_(*(std::get<indexMatrix>(this->particles_->all_particle_data_)[this->particles_->all_variable_maps_[indexMatrix]["Korteweg_tensor"]])) {}
//=================================================================================================//
	void KortewegAccelerationInner::Interaction(size_t index_i, Real dt)
	{
		Real rho_i = rho_n_[index_i];
		const Vecd &vel_i = vel_n_[index_i];
		Vecd acceleration(0), vel_derivative(0);
		const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
		for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
		{
			size_t index_j = inner_neighborhood.j_[n];
			// Kortewegs force
			acceleration += rho_n_[index_j] * Vol_[index_j] *
							(korteweg_tensor_[index_i] / rho_i / rho_i - korteweg_tensor_[index_j] / rho_n_[index_j] / rho_n_[index_j]) * inner_neighborhood.e_ij_[n] * inner_neighborhood.dW_ij_[n];
		}
		dvel_dt_prior_[index_i] += acceleration;
	}
}
#endif
