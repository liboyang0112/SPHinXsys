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
#ifndef SMALL_VECTORS_H
#define SMALL_VECTORS_H

#include "base_data_type.h"

namespace SPH {

	Vec2d FirstAxisVector(const Vec2d& zero_vector);
	Vec3d FirstAxisVector(const Vec3d& zero_vector);
	Real getMaxAbsoluteElement(const Vec2d& input);
	Real getMaxAbsoluteElement(const Vec3d& input);
	Vec3d upgradeToVector3D(const Real& input);
	Vec3d upgradeToVector3D(const Vec2d& input);
	Vec3d upgradeToVector3D(const Vec3d& input);
	Mat3d upgradeToMatrix3D(const Mat2d& input);
	Mat3d upgradeToMatrix3D(const Mat3d& input);

	template<typename OutVectorType>
	OutVectorType upgradeVector(const Real& input)
	{
		OutVectorType out_vector(0);
		out_vector[0] = input;
		return out_vector;
	};
	template<typename OutVectorType>
	OutVectorType upgradeVector(const Vec2d& input)
	{
		OutVectorType out_vector(0);
		out_vector[0] = input[0];
		out_vector[1] = input[1];
		return out_vector;
	};
	template<typename OutVectorType>
	OutVectorType upgradeVector(const Vec3d& input)
	{
		OutVectorType out_vector(0);
		out_vector[0] = input[0];
		out_vector[1] = input[1];
		out_vector[2] = input[2];
		return out_vector;
	};

	Mat2d getInverse(const Mat2d& A);
	Mat3d getInverse(const Mat3d& A);
	Mat2d getAverageValue(const Mat2d &A, const Mat2d& B);
	Mat3d getAverageValue(const Mat3d &A, const Mat3d& B);
	Mat2d inverseCholeskyDecomposition(const Mat2d& A);
	Mat3d inverseCholeskyDecomposition(const Mat3d& A);
	Mat2d getDiagonal(const Mat2d& A);
	Mat3d getDiagonal(const Mat3d& A);

	/** get transformation matrix. */
	Mat2d getTransformationMatrix(const Vec2d& direction_of_y);
	Mat3d getTransformationMatrix(const Vec3d& direction_of_z);

	/** get angle between two vectors. */
	Real getCosineOfAngleBetweenTwoVectors (Vec2d vector_1, Vec2d vector_2);
	Real getCosineOfAngleBetweenTwoVectors (Vec3d vector_1, Vec3d vector_2);

	/** get orthogonal projection of a vactor. */
	Vec2d getVectorProjectionOfVector (Vec2d vector_1, Vec2d vector_2);
	Vec3d getVectorProjectionOfVector (Vec3d vector_1, Vec3d vector_2);
}

#endif //SMALL_VECTORS_H
