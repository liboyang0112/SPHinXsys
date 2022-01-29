#include "triangle_mesh_shape.h"

namespace SPH
{
	//=================================================================================================//
	SimTK::ContactGeometry::TriangleMesh *TriangleMeshShape::
		generateTriangleMesh(SimTK::PolygonalMesh &ploy_mesh)
	{
		SimTK::ContactGeometry::TriangleMesh *triangle_mesh;
		triangle_mesh = triangle_mesh_ptr_keeper_.createPtr<SimTK::ContactGeometry::TriangleMesh>(ploy_mesh);
		if (!SimTK::ContactGeometry::TriangleMesh::isInstance(*triangle_mesh))
		{
			std::cout << "\n Error the triangle mesh is not valid" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			throw;
		}
		std::cout << "num of faces:" << triangle_mesh->getNumFaces() << std::endl;

		return triangle_mesh;
	}
	//=================================================================================================//
	bool TriangleMeshShape::checkContain(const Vec3d &pnt, bool BOUNDARY_INCLUDED)
	{

		SimTK::Vec2 uv_coordinate;
		bool inside = false;
		int face_id;
		Vec3d closest_pnt = triangle_mesh_->findNearestPoint(pnt, inside, face_id, uv_coordinate);

		StdVec<int> neighbor_face(4);
		neighbor_face[0] = face_id;
		/** go throught the neighbor faces. */
		for (int i = 1; i < 4; i++)
		{
			int edge = triangle_mesh_->getFaceEdge(face_id, i - 1);
			int face = triangle_mesh_->getEdgeFace(edge, 0);
			neighbor_face[i] = face != face_id ? face : triangle_mesh_->getEdgeFace(edge, 1);
		}

		Vec3d from_face_to_pnt = pnt - closest_pnt;
		Real sum_weights = 0.0;
		Real weighted_dot_product = 0.0;
		for (int i = 0; i < 4; i++)
		{
			SimTK::UnitVec3 normal_direction = triangle_mesh_->getFaceNormal(neighbor_face[i]);
			Real dot_product = dot(normal_direction, from_face_to_pnt);
			Real weight = dot_product * dot_product;
			weighted_dot_product += weight * dot_product;
			sum_weights += weight;
		}

		weighted_dot_product /= sum_weights;

		bool weighted_inside = false;
		if (weighted_dot_product < 0.0)
			weighted_inside = true;

		if (face_id < 0 && face_id > triangle_mesh_->getNumFaces())
		{
			std::cout << "\n Error the nearest point is not valid" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			throw;
		}

		return weighted_inside;
	}
	//=================================================================================================//
	Vec3d TriangleMeshShape::findClosestPoint(const Vec3d &input_pnt)
	{
		bool inside = false;
		int face_id;
		SimTK::Vec2 normal;
		Vec3d closest_pnt;
		closest_pnt = triangle_mesh_->findNearestPoint(input_pnt, inside, face_id, normal);
		if (face_id < 0 && face_id > triangle_mesh_->getNumFaces())
		{
			std::cout << "\n Error the nearest point is not valid" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			throw;
		}
		return closest_pnt;
	}
	//=================================================================================================//
	BoundingBox TriangleMeshShape::findBounds()
	{
		int number_of_vertices = triangle_mesh_->getNumVertices();
		//initial reference values
		Vec3d lower_bound = Vec3d(Infinity);
		Vec3d upper_bound = Vec3d(-Infinity);

		for (int i = 0; i != number_of_vertices; ++i)
		{
			Vec3d vertex_position = triangle_mesh_->getVertexPosition(i);
			for (int j = 0; j != 3; ++j)
			{
				lower_bound[j] = SMIN(lower_bound[j], vertex_position[j]);
				upper_bound[j] = SMAX(upper_bound[j], vertex_position[j]);
			}
		}
		return BoundingBox(lower_bound, upper_bound);
	}
	//=================================================================================================//
	TriangleMeshShapeSTL::
		TriangleMeshShapeSTL(const std::string &filepathname, Vec3d translation, Real scale_factor,
							 const std::string &shape_name)
		: TriangleMeshShape(shape_name)
	{
		if (!fs::exists(filepathname))
		{
			std::cout << "\n Error: the input file:" << filepathname << " is not exists" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			throw;
		}
		SimTK::PolygonalMesh polymesh;
		polymesh.loadStlFile(filepathname);
		polymesh.scaleMesh(scale_factor);
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(translation));
	}
	//=================================================================================================//
	TriangleMeshShapeBrick::
		TriangleMeshShapeBrick(Vec3d halfsize, int resolution, Vec3d translation,
							   const std::string &shape_name)
		: TriangleMeshShape(shape_name)
	{
		SimTK::PolygonalMesh polymesh = SimTK::PolygonalMesh::createBrickMesh(halfsize, resolution);
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(translation));
	}
	//=================================================================================================//
	TriangleMeshShapeBrick::
		TriangleMeshShapeBrick(const TriangleMeshShapeBrick::ShapeParameters &shape_parameters,
							   const std::string &shape_name)
		: TriangleMeshShapeBrick(shape_parameters.halfsize_, shape_parameters.resolution_,
								 shape_parameters.translation_, shape_name) {}
	//=================================================================================================//
	TriangleMeshShapeShere::TriangleMeshShapeShere(Real radius, int resolution, Vec3d translation,
												   const std::string &shape_name)
		: TriangleMeshShape(shape_name)
	{
		SimTK::PolygonalMesh polymesh = SimTK::PolygonalMesh::createSphereMesh(radius, resolution);
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(translation));
	}
	//=================================================================================================//
	TriangleMeshShapeCylinder::
		TriangleMeshShapeCylinder(SimTK::UnitVec3 axis, Real radius, Real halflength, int resolution,
								  Vec3d translation, const std::string &shape_name)
		: TriangleMeshShape(shape_name)
	{
		SimTK::PolygonalMesh polymesh =
			SimTK::PolygonalMesh::createCylinderMesh(axis, radius, halflength, resolution);
		triangle_mesh_ = generateTriangleMesh(polymesh.transformMesh(translation));
	}
	//=================================================================================================//
}