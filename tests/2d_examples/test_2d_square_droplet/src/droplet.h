/**
 * @file 	droplet.h
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
Real DL = 2.0;						   /**< Tank length. */
Real DH = 2.0;						   /**< Tank height. */
Real LL = 1.0;						   /**< Liquid colume length. */
Real LH = 1.0;						   /**< Liquid colume height. */
Real particle_spacing_ref = DL / 40.0; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 2;	   /**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1.0;		  /**< Reference density of water. */
Real rho0_a = 0.001;	  /**< Reference density of air. */
Real gravity_g = 0.0;	  /**< Gravity force of fluid. */
Real U_max = 1.0;		  /**< Characteristic velocity. */
Real c_f = 10.0 * U_max; /**< Reference sound speed. */
Real mu_f = 0.2;		  /**< Water viscosity. */
Real mu_a = 0.0002;		  /**< Air viscosity. */
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(0.5, 0.5));
	water_block_shape.push_back(Vecd(0.5, LH + 0.5));
	water_block_shape.push_back(Vecd(LL + 0.5, LH + 0.5));
	water_block_shape.push_back(Vecd(LL + 0.5, 0.5));
	water_block_shape.push_back(Vecd(0.5, 0.5));
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
*@brief 	Water body definition.
*/
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &sph_system, const std::string &body_name)
		: FluidBody(sph_system, body_name, makeShared<SPHAdaptation>(1.3, 1))
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
		: FluidBody(sph_system, body_name, makeShared<SPHAdaptation>(1.3, 1.0))
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
	WallBoundary(SPHSystem &sph_system, const std::string &body_name)
		: SolidBody(sph_system, body_name, makeShared<SPHAdaptation>(1.3, 1))
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
