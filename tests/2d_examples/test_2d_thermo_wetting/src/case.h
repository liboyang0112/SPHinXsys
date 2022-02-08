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
Real phi_fluid_initial = 40.0;
Real phi_gas_initial = 40.0;

/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;					/**< Reference density of water. */
Real rho0_a = 1.0e-3;				/**< Reference density of air. */
Real gravity_g = 0.0;				/**< Gravity force of fluid. */
Real U_max = 1.0;					/**< Characteristic velocity. */
Real c_f = 10.0 * U_max;			/**< Reference sound speed. */
Real mu_f = 5.0e-2;					/**< Water viscosity. */
Real mu_a = 5.0e-5;					/**< Air visocsity. */
Real contact_angle = (150.0 / 180.0) * 3.1415926; 	/**< Contact angle with Wall. */
Real tension_force = 0.008;
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
		phi_ = material_->SpeciesIndexMap()["Phi"];
	};
};

/**
 * application dependent fluid body initial condition
 */
class ThermofluidBodyInitialCondition
	: public  DiffusionReactionInitialCondition< FluidBody, FluidParticles, WeaklyCompressibleFluid>
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
		: DiffusionReactionInitialCondition<FluidBody, FluidParticles, WeaklyCompressibleFluid >(diffusion_fluid_body) {
		phi_ = material_->SpeciesIndexMap()["Phi"];
	};
};

/**
 * application dependent gas body initial condition
 */
class ThermogasBodyInitialCondition
	: public  DiffusionReactionInitialCondition< FluidBody, FluidParticles, WeaklyCompressibleFluid>
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
		: DiffusionReactionInitialCondition<FluidBody, FluidParticles, WeaklyCompressibleFluid >(diffusion_fluid_body) {
		phi_ = material_->SpeciesIndexMap()["Phi"];
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
		  FluidBody, FluidParticles, WeaklyCompressibleFluid,
		  RelaxationOfAllDiffussionSpeciesComplex<
			  FluidBody, FluidParticles, WeaklyCompressibleFluid, SolidBody, SolidParticles, Solid>,
		  ComplexBodyRelation>
{
public:
	explicit ThermalRelaxationComplex(ComplexBodyRelation &body_complex_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_complex_relation){};
	virtual ~ThermalRelaxationComplex(){};
};


class ThermalRelaxationComplexWA
	: public RelaxationOfAllDiffusionSpeciesRK2<
		  FluidBody, FluidParticles, WeaklyCompressibleFluid,
		  RelaxationOfAllDiffussionSpeciesComplex<
			  FluidBody, FluidParticles, WeaklyCompressibleFluid, FluidBody, FluidParticles, WeaklyCompressibleFluid>,
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


class AirMaterial
	: public DiffusionReaction<FluidParticles, WeaklyCompressibleFluid>
{
public:
	AirMaterial()
		: DiffusionReaction<FluidParticles, WeaklyCompressibleFluid>({"Phi"}, rho0_a, c_f, mu_a)
	{
		initializeAnDiffusion<DirectionalDiffusion>("Phi", "Phi", diffusion_coff, bias_coff, bias_direction);
	};
};

/**
 * Setup heat conduction material properties for diffusion solid body
 */
class WallMaterial
	: public DiffusionReaction<SolidParticles, Solid>
{
public:
	WallMaterial()
		: DiffusionReaction<SolidParticles, Solid>({"Phi"})
	{
		initializeAnDiffusion<DirectionalDiffusion>("Phi", "Phi", diffusion_coff, bias_coff, bias_direction);
	};
};