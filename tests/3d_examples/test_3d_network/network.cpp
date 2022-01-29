/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
 * @file 	network.cpp
 * @brief 	This is the example of generating a neural network on a sphere 
 * @author 	Chi Zhang and Xiangyu Hu
 */
/** header file and namespace. */
#include "sphinxsys.h"
using namespace SPH;

/** Set the file path to the stl file. */
std::string full_path_to_stl = "./input/sphere.stl";
Vec3d domain_lower_bound(-1.0, -1.0, -1.0);
Vec3d domain_upper_bound(1.0, 1.0, 1.0);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 100.0;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);

/** Define the my heart body. */
class MyPolygonBody : public SolidBody
{
public:
	MyPolygonBody(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name, makeShared<SPHAdaptation>(1.15, 1.0))
	{
		Vecd translation(-1.0, -1.0, -1.0);
		TriangleMeshShapeSTL triangle_mesh_shape(full_path_to_stl, translation, 0.025);
		body_shape_.add<LevelSetShape>(this, triangle_mesh_shape);
	}
};
/**
 *  The main program
 */
int main()
{
	/** Setup the system. */
	SPHSystem system(system_domain_bounds, dp_0);
	/** Output */
	In_Output in_output(system);
	/** Creat a sphere, corresponding material and particles. */
	MyPolygonBody polygon_body(system, "Polygon");
	SolidParticles body_particles(
		polygon_body, makeShared<ParticleGeneratorNetwork>(Vecd(-1.0, 0.0, 0.0), Vecd(-0.964, 0.0, 0.266), 15, 5.0));
	/** Write particle data. */
	BodyStatesRecordingToVtp write_states(in_output, system.real_bodies_);
	write_states.writeToFile(0);
	return 0;
}