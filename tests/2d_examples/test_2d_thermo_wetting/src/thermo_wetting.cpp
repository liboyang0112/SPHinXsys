/**
 * @file 	two_phase_dambreak.cpp
 * @brief 	2D two-phase dambreak flow.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for multi-phase simulation.
 * @author 	Chi Zhang, Xiangyu Hu and Boyang Li
 */
#include "case.h"
#include "sphinxsys.h"
#include "SPHconfig.h"
using namespace SPH;
/**
 * @brief 	Main program starts here.
 */

struct vdWFluidPar
{
	Real T0;
	Real rho0;
	Real alpha;
	Real rho_m;
	Real mu;
	Real molmass;
	int dof;
	Real gamma;
};
struct SolidPar
{
	Real T0;
};

int main(int argc, char* argv[])
{
	/**
	 * @brief Build up -- a SPHSystem --
	 */
	SPHconfig cfg;
	vdWFluidPar water;
	vdWFluidPar air;
	SolidPar wall;
	(*cfg.vdWFluids)[0].lookupValue("viscousity",air.mu);
	(*cfg.vdWFluids)[0].lookupValue("rho_0",air.rho0);
	(*cfg.vdWFluids)[0].lookupValue("a",air.alpha);
	(*cfg.vdWFluids)[0].lookupValue("rho_m",air.rho_m);
	(*cfg.vdWFluids)[0].lookupValue("gamma",air.gamma);
	(*cfg.vdWFluids)[1].lookupValue("temperature",air.T0);
	(*cfg.vdWFluids)[1].lookupValue("viscousity",water.mu);
	(*cfg.vdWFluids)[1].lookupValue("rho_0",water.rho0);
	(*cfg.vdWFluids)[1].lookupValue("a",water.alpha);
	(*cfg.vdWFluids)[1].lookupValue("rho_m",water.rho_m);
	(*cfg.vdWFluids)[1].lookupValue("gamma",water.gamma);
	(*cfg.vdWFluids)[1].lookupValue("temperature",water.T0);
	(*cfg.Solids)[0].lookupValue("temperature",wall.T0);
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation from restart files. 0: not from restart files. */
	sph_system.restart_step_ = 0;
	/** I/O environment. */
	In_Output in_output(sph_system);
	/**
	 * @brief Material property, partilces and body creation of water.
	 */
	printf("%f, %f, %f, %f, %f, %f\n", water.rho0, water.rho_m, water.gamma, water.alpha, water.molmass, water.mu);

	WaterBlock water_block(sph_system, "WaterBody");
	DiffusionReactionParticles<CompressibleFluidParticles, vdWFluid>
		diffusion_fluid_body_particles(water_block, makeShared<FluidMaterial>(water.rho0, water.rho_m, water.gamma, water.alpha, water.molmass, water.mu));
	diffusion_fluid_body_particles.addAVariableToWrite<indexScalar, Real>("Density");
	WallBoundary wall_boundary(sph_system, "Wall");
	DiffusionReactionParticles<SolidParticles, Solid>
		wall_particles(wall_boundary, makeShared<WallMaterial>());

	ObserverBody temperature_observer(sph_system, "FluidObserver");
	ObserverParticles temperature_observer_particles(temperature_observer, makeShared<ObserverParticleGenerator>());
	/**
	 * @brief Material property, partilces and body creation of air.
	 */
	/** topology */
	//BodyRelationInnerMultiLength fluid_body_inner(water_block); //TBD
	//BodyRelationInnerMultiLength solid_body_inner(wall_boundary); //TBD
	BodyRelationInner fluid_body_inner(water_block);
	BodyRelationInner solid_body_inner(wall_boundary);
	ComplexBodyRelation water_wall_complex(fluid_body_inner, { &wall_boundary });
	BaseBodyRelationContact &water_wall_contact = water_wall_complex.contact_relation_;
	BodyRelationContact fluid_observer_contact(temperature_observer, {&water_block});
	/**
	 * @brief 	Define all numerical methods which are used in this case.
	 */
	 /** Define external force. */
	Gravity		gravity(Vecd(0.0, -gravity_g));
	HeatSource heat_source;
	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	ThermosolidBodyInitialCondition thermosolid_condition(wall_boundary, wall.T0);
	ThermofluidBodyInitialCondition thermofluid_initial_condition(water_block, water.T0);

	TimeStepInitialization		initialize_a_water_step(water_block, gravity);

	DiffusionSourceInitialization<FluidBody, CompressibleFluidParticles, vdWFluid> initialize_a_water_step_thermo(water_block, heat_source, 0);
	StressTensorHeatSource<FluidBody, CompressibleFluidParticles, vdWFluid>  stress_tensor_heat_water(fluid_body_inner, 0);
	vdWAttractionHeatSource<FluidBody, CompressibleFluidParticles, vdWFluid>  vdW_attr_heat_water(water_block, 0);
	/** Corrected strong configuration for diffusion solid body. */
	solid_dynamics::CorrectConfiguration 			correct_configuration(solid_body_inner);
	/** Time step size calculation. */
	GetDiffusionTimeStepSize<FluidBody, CompressibleFluidParticles, vdWFluid> get_thermal_time_step(water_block);
	/** Diffusion process between three diffusion bodies. */
	//ThermalRelaxationComplexWA 	thermal_relaxation_complex_wa(water_air_complex);
	ThermalRelaxationComplex 	thermal_relaxation_complex_ww(water_wall_complex);
	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	 /** Evaluation of density by summation approach. */
	fluid_dynamics::DensitySummationFreeSurfaceComplexWithoutUpdate 	
		update_water_density_by_summation(fluid_body_inner, water_wall_contact);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize 	get_water_advection_time_step_size(water_block, U_max);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_water_time_step_size(water_block);
	/** Pressure relaxation for water by using position verlet time stepping. */
	fluid_dynamics::vdWPressureRelaxationRiemannWithWall 
		water_pressure_relaxation(fluid_body_inner, water_wall_contact);
	fluid_dynamics::vdWDensityRelaxationRiemannWithWall 
		water_density_relaxation(fluid_body_inner, water_wall_contact);
	//fluid_dynamics::ViscousAccelerationMultiPhase
	fluid_dynamics::ViscousAccelerationInner
		water_viscou_acceleration(fluid_body_inner);
	/** Suface tension and wetting effects. */
	fluid_dynamics::FreeSurfaceIndicationComplex
		surface_detection(fluid_body_inner, water_wall_contact);
	fluid_dynamics::ColorFunctionGradientComplex
		color_gradient(fluid_body_inner, water_wall_contact);
	fluid_dynamics::ColorFunctionGradientInterplationInner
		color_gradient_interpolation(fluid_body_inner);
	fluid_dynamics::SurfaceTensionAccelerationInner
		surface_tension_acceleration(fluid_body_inner, tension_force);
	/** Wetting effects. */
	fluid_dynamics::SurfaceNormWithWall
		wetting_norm(water_wall_contact, contact_angle);
	/** Computing vorticity in the flow. */
	fluid_dynamics::VorticityInner 	compute_vorticity(fluid_body_inner);
	/** Output the body states. */
	BodyStatesRecordingToVtp 		body_states_recording(in_output, sph_system.real_bodies_);
	RegressionTestEnsembleAveraged<ObservedQuantityRecording<indexScalar, Real>>
		write_fluid_phi("Temperature", in_output, fluid_observer_contact);
	ObservedQuantityRecording<indexVector, Vecd>
		write_fluid_velocity("Velocity", in_output, fluid_observer_contact);
	/** Output the body states for restart simulation. */
	RestartIO		restart_io(in_output, sph_system.real_bodies_);
	/** Pre-simulation*/
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromBodyShape();
	/**
	 * @brief The time stepping starts here.
	 */
	correct_configuration.parallel_exec();
	thermosolid_condition.parallel_exec();
	thermofluid_initial_condition.parallel_exec();
	Real dt_thermal = get_thermal_time_step.parallel_exec();
	 /** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (sph_system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.restart_step_);
		water_block.updateCellLinkedList();
		fluid_body_inner.updateConfiguration();
		water_wall_contact.updateConfiguration();
		//air_water_complex.updateConfiguration();
		//air_wall_contact.updateConfiguration();
	}
	/** Output the start states of bodies. */
	body_states_recording.writeToFile(0);
	/**
	 * @brief 	Basic parameters.
	 */
	size_t number_of_iterations = sph_system.restart_step_;
	int screen_output_interval = 30;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 1000; 	/**< End time. */
	Real D_Time = 1;		/**< Time stamps for output of body states. */
	Real Dt = 0.0;			/**< Default advection time step sizes. */
	Real dt = 0.0; 			/**< Default acoustic time step sizes. */
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	tick_count::interval_t interval_computing_time_step;
	tick_count::interval_t interval_computing_pressure_relaxation;
	tick_count::interval_t interval_updating_configuration;
	tick_count time_instance;
	bool firstiter = 1;
	bool denseWriting = 0;
	int iframe = 1;
	/**
	 * @brief 	Main loop starts here.
	 */

	update_water_density_by_summation.parallel_exec(1);
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** Acceleration due to viscous force and gravity. */
			time_instance = tick_count::now();
			//initialize_a_air_step.parallel_exec();

			Real Dt_f = get_water_advection_time_step_size.parallel_exec();
			//Real Dt_a = get_air_advection_time_step_size.parallel_exec();
			Dt = Dt_f;//SMIN(Dt_f, Dt_a);

			//update_air_density_by_summation.parallel_exec();
			//air_transport_correction.parallel_exec(Dt);

			//air_viscou_acceleration.parallel_exec();

			//surface_detection.parallel_exec();
			//color_gradient.parallel_exec();
			//color_gradient_interpolation.parallel_exec();
			//wetting_norm.parallel_exec();
			//surface_tension_acceleration.parallel_exec();

			interval_computing_time_step += tick_count::now() - time_instance;

			//initialize_a_air_step_thermo.parallel_exec();
			/** Dynamics including pressure relaxation. */
			time_instance = tick_count::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt_f = get_water_time_step_size.parallel_exec();
				//Real dt_a = get_air_time_step_size.parallel_exec();
				dt = SMIN(SMIN(dt_f, Dt),dt_thermal);
				if(iframe==89){
					int halt = 1;
				}
				initialize_a_water_step_thermo.parallel_exec();
				initialize_a_water_step.parallel_exec();
				water_viscou_acceleration.parallel_exec();
				water_pressure_relaxation.parallel_exec(dt);
				//air_pressure_relaxation.parallel_exec(dt);

				update_water_density_by_summation.parallel_exec(dt);
				water_density_relaxation.parallel_exec(dt);
				//air_density_relaxation.parallel_exec(dt);

				//thermal_relaxation_complex_wa.parallel_exec(dt);
				stress_tensor_heat_water.parallel_exec();
				vdW_attr_heat_water.parallel_exec();
				thermal_relaxation_complex_ww.parallel_exec(dt);
				//thermal_relaxation_complex_aw.parallel_exec(dt);
				//stress_tensor_heat_air.parallel_exec();
				//vdW_attr_heat_air.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				//if(firstiter) {
				//	body_states_recording.writeToFile();
				//	firstiter =0;
				//}
				water_block.updateCellLinkedList();
				water_wall_complex.updateConfiguration();
				if(denseWriting) {
					body_states_recording.writeToFile();
					iframe++;
				}
			}
			interval_computing_pressure_relaxation += tick_count::now() - time_instance;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";

				if(!denseWriting) body_states_recording.writeToFile();
				if (number_of_iterations % restart_output_interval == 0)
					restart_io.writeToFile(number_of_iterations);
			}
			number_of_iterations++;
			time_instance = tick_count::now();
			

			/** Update cell linked list and configuration. */

			interval_updating_configuration += tick_count::now() - time_instance;
		}


		tick_count t2 = tick_count::now();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;

	}

	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
		<< " seconds." << std::endl;
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
		<< interval_computing_time_step.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
		<< interval_computing_pressure_relaxation.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
		<< interval_updating_configuration.seconds() << "\n";

	return 0;
}
