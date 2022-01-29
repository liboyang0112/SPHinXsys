/**
 * @file 	load_image.cpp
 * @brief 	This is the test of using distance map to generate body fitted particles (3D).
 * @details We use this case to test the particle generation and relaxation for a complex geometry. 
 *			Before particle generation, we clean the sharp corners of the model. 
 * @author 	Yijin Mao, Yongchuan Yu and Xiangyu Hu
 */

#include "sphinxsys.h"

/** case file to setup the test case */
#include "load_image.h"

using namespace SPH;

int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up -- a SPHSystem
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, dp_0);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = true;
//handle command line arguments
#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(ac, av);
#endif
	/** output environment. */
	In_Output in_output(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	ImportedModel imported_model(system, "ImportedModel");
	SolidParticles imported_model_particles(imported_model, makeShared<ParticleGeneratorMultiResolution>());
	imported_model_particles.addAVariableToWrite<indexScalar, Real>("SmoothingLengthRatio");
	//----------------------------------------------------------------------
	//	Define simple file input and outputs functions.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_imported_model_to_vtp(in_output, {imported_model});
	MeshRecordingToPlt cell_linked_list_recording(in_output, imported_model, imported_model.cell_linked_list_);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInnerVariableSmoothingLength imported_model_inner(imported_model);
	//BaseBodyRelationInner* imported_model_inner(imported_model);
	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	RandomizePartilePosition random_imported_model_particles(imported_model);
	/** A  Physics relaxation step. */
	//relax_dynamics::RelaxationStepInner relaxation_step_inner(imported_model_inner.get(), true);
	relax_dynamics::RelaxationStepInner relaxation_step_inner(imported_model_inner, true);
	relax_dynamics::UpdateSmoothingLengthRatioByBodyShape update_smoothing_length_ratio(imported_model);
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	random_imported_model_particles.parallel_exec(0.25);
	relaxation_step_inner.surface_bounding_.parallel_exec();
	update_smoothing_length_ratio.parallel_exec();
	write_imported_model_to_vtp.writeToFile();
	imported_model.updateCellLinkedList();
	cell_linked_list_recording.writeToFile(0);
	//----------------------------------------------------------------------
	//	Particle relaxation time stepping start here.
	//----------------------------------------------------------------------
	int ite_p = 0;
	while (ite_p < 1000)
	{
		//update_smoothing_length_ratio.parallel_exec();
		relaxation_step_inner.parallel_exec();
		ite_p += 1;
		if (ite_p % 100 == 0)
		{
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the imported model N = " << ite_p << "\n";
			write_imported_model_to_vtp.writeToFile(ite_p);
		}
	}
	std::cout << "The physics relaxation process of imported model finish !" << std::endl;

	return 0;
}
