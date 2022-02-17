/**
 * @file 	two_phase_dambreak.cpp
 * @brief 	2D two-phase dambreak flow.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for multi-phase simulation.
 * @author 	Chi Zhang, Xiangyu Hu and Boyang Li
 */
#include "case.h"
#include "sphinxsys.h"
using namespace SPH;
/**
 * @brief 	Main program starts here.
 */


int main()
{
	/**
	 * @brief Build up -- a SPHSystem --
	 */
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

	AirBlock air_block(sph_system, "AirBody");
	DiffusionReactionParticles<CompressibleFluidParticles, vdWFluid>
		diffusion_gas_body_particles(air_block, makeShared<AirMaterial>());

	WaterBlock water_block(sph_system, "WaterBody");
	DiffusionReactionParticles<CompressibleFluidParticles, vdWFluid>
		diffusion_fluid_body_particles(water_block, makeShared<WaterMaterial>());

	WallBoundary wall_boundary(sph_system, "Wall");
	DiffusionReactionParticles<SolidParticles, Solid>
		wall_particles(wall_boundary, makeShared<WallMaterial>());

	ObserverBody temperature_observer(sph_system, "FluidObserver");
	ObserverParticles temperature_observer_particles(temperature_observer, makeShared<ObserverParticleGenerator>());
	/**
	 * @brief Material property, partilces and body creation of air.
	 */
	/** topology */
	BodyRelationInner fluid_body_inner(water_block);
	BodyRelationInner gas_body_inner(air_block);
	BodyRelationInner solid_body_inner(wall_boundary);
	ComplexBodyRelation water_air_complex(water_block, {&air_block});
	ComplexBodyRelation water_wall_complex(fluid_body_inner, { &wall_boundary });
	BaseBodyRelationContact &water_wall_contact = water_wall_complex.contact_relation_;
	ComplexBodyRelation air_water_complex(air_block, { &water_block });
	ComplexBodyRelation air_wall_complex(gas_body_inner, { &wall_boundary });
	BaseBodyRelationContact &air_wall_contact = air_wall_complex.contact_relation_;
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
	ThermosolidBodyInitialCondition thermosolid_condition(wall_boundary);
	ThermogasBodyInitialCondition thermogas_initial_condition(air_block);
	ThermofluidBodyInitialCondition thermofluid_initial_condition(water_block);

	TimeStepInitialization		initialize_a_water_step(water_block, gravity);
	TimeStepInitialization		initialize_a_air_step(air_block, gravity);

	DiffusionSourceInitialization<FluidBody, CompressibleFluidParticles, vdWFluid> initialize_a_air_step_thermo(air_block, heat_source, 0);
	DiffusionSourceInitialization<FluidBody, CompressibleFluidParticles, vdWFluid> initialize_a_water_step_thermo(water_block, heat_source, 0);
	StressTensorHeatSource<FluidBody, CompressibleFluidParticles, vdWFluid>  stress_tensor_heat_air(gas_body_inner, 0);
	StressTensorHeatSource<FluidBody, CompressibleFluidParticles, vdWFluid>  stress_tensor_heat_water(fluid_body_inner, 0);
	vdWAttractionHeatSource<FluidBody, CompressibleFluidParticles, vdWFluid>  vdW_attr_heat_air(air_block, 0);
	vdWAttractionHeatSource<FluidBody, CompressibleFluidParticles, vdWFluid>  vdW_attr_heat_water(water_block, 0);
	/** Corrected strong configuration for diffusion solid body. */
	solid_dynamics::CorrectConfiguration 			correct_configuration(solid_body_inner);
	/** Time step size calculation. */
	GetDiffusionTimeStepSize<FluidBody, CompressibleFluidParticles, vdWFluid> get_thermal_time_step(water_block);
	GetDiffusionTimeStepSize<FluidBody, CompressibleFluidParticles, vdWFluid> get_thermal_time_step_air(air_block);
	/** Diffusion process between three diffusion bodies. */
	//ThermalRelaxationComplexWA 	thermal_relaxation_complex_wa(water_air_complex);
	ThermalRelaxationComplex 	thermal_relaxation_complex_ww(water_wall_complex);
	ThermalRelaxationComplex 	thermal_relaxation_complex_aw(air_wall_complex);
	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	 /** Evaluation of density by summation approach. */
	fluid_dynamics::DensitySummationFreeSurfaceComplex 	
		update_water_density_by_summation(water_air_complex.inner_relation_, water_wall_contact);
	fluid_dynamics::DensitySummationComplex
		update_air_density_by_summation(air_water_complex, air_wall_contact);
	fluid_dynamics::TransportVelocityCorrectionComplex
		air_transport_correction(air_water_complex, air_wall_contact);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize 	get_water_advection_time_step_size(water_block, U_max);
	fluid_dynamics::AdvectionTimeStepSize 	get_air_advection_time_step_size(air_block, U_max);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_water_time_step_size(water_block);
	fluid_dynamics::AcousticTimeStepSize get_air_time_step_size(air_block);
	/** Pressure relaxation for water by using position verlet time stepping. */
	fluid_dynamics::vdWPressureRelaxationRiemannWithWall 
		water_pressure_relaxation(water_air_complex.inner_relation_, water_wall_contact);
	fluid_dynamics::vdWDensityRelaxationRiemannWithWall 
		water_density_relaxation(water_air_complex.inner_relation_, water_wall_contact);
	/** Extend Pressure relaxation is used for air. */
	fluid_dynamics::vdWExtendMultiPhasePressureRelaxationRiemannWithWall
		air_pressure_relaxation(air_water_complex, air_wall_contact, 2.0);
	fluid_dynamics::vdWMultiPhaseDensityRelaxationRiemannWithWall
		air_density_relaxation(air_water_complex, air_wall_contact);
	/** Viscous acceleration. */
	fluid_dynamics::ViscousAccelerationMultiPhase
		air_viscou_acceleration(air_water_complex);
	fluid_dynamics::ViscousAccelerationMultiPhase
		water_viscou_acceleration(water_air_complex);
	/** Suface tension and wetting effects. */
	fluid_dynamics::FreeSurfaceIndicationComplex
		surface_detection(water_air_complex.inner_relation_, water_wall_contact);
	fluid_dynamics::ColorFunctionGradientComplex
		color_gradient(water_air_complex.inner_relation_, water_wall_contact);
	fluid_dynamics::ColorFunctionGradientInterplationInner
		color_gradient_interpolation(water_air_complex.inner_relation_);
	fluid_dynamics::SurfaceTensionAccelerationInner
		surface_tension_acceleration(water_air_complex.inner_relation_, tension_force);
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
	correct_configuration.exec();
	thermosolid_condition.exec();
	thermofluid_initial_condition.exec();
	thermogas_initial_condition.exec();
	Real dt_thermal = SMIN(get_thermal_time_step.exec(),get_thermal_time_step_air.exec());
	 /** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (sph_system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.restart_step_);
		water_block.updateCellLinkedList();
		air_block.updateCellLinkedList();
		water_air_complex.updateConfiguration();
		water_wall_contact.updateConfiguration();
		air_water_complex.updateConfiguration();
		air_wall_contact.updateConfiguration();
	}
	/** Output the start states of bodies. */
	body_states_recording.writeToFile(0);
	/**
	 * @brief 	Basic parameters.
	 */
	size_t number_of_iterations = sph_system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 50; 	/**< End time. */
	Real D_Time = End_Time / 1000;		/**< Time stamps for output of body states. */
	Real Dt = 0.0;			/**< Default advection time step sizes. */
	Real dt = 0.0; 			/**< Default acoustic time step sizes. */
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	tick_count::interval_t interval_computing_time_step;
	tick_count::interval_t interval_computing_pressure_relaxation;
	tick_count::interval_t interval_updating_configuration;
	tick_count time_instance;

	/**
	 * @brief 	Main loop starts here.
	 */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** Acceleration due to viscous force and gravity. */
			time_instance = tick_count::now();
			initialize_a_water_step.exec();
			initialize_a_air_step.exec();

			Real Dt_f = get_water_advection_time_step_size.exec();
			Real Dt_a = get_air_advection_time_step_size.exec();
			Dt = SMIN(Dt_f, Dt_a);

			update_water_density_by_summation.exec();
			update_air_density_by_summation.exec();
			air_transport_correction.exec(Dt);

			air_viscou_acceleration.exec();
			water_viscou_acceleration.exec();

			surface_detection.exec();
			color_gradient.exec();
			color_gradient_interpolation.exec();
			//wetting_norm.exec();
			//surface_tension_acceleration.exec();

			interval_computing_time_step += tick_count::now() - time_instance;

			initialize_a_air_step_thermo.exec();
			initialize_a_water_step_thermo.exec();
			stress_tensor_heat_air.exec();
			stress_tensor_heat_water.exec();
			vdW_attr_heat_air.exec();
			vdW_attr_heat_water.exec();
			/** Dynamics including pressure relaxation. */
			time_instance = tick_count::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt_f = get_water_time_step_size.exec();
				Real dt_a = get_air_time_step_size.exec();
				dt = SMIN(SMIN(SMIN(dt_f, dt_a), Dt),dt_thermal);
				water_pressure_relaxation.exec(dt);
				air_pressure_relaxation.exec(dt);

				water_density_relaxation.exec(dt);
				air_density_relaxation.exec(dt);

				//thermal_relaxation_complex_wa.exec(dt);
				thermal_relaxation_complex_ww.exec(dt);
				thermal_relaxation_complex_aw.exec(dt);
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
		body_states_recording.writeToFile();
			}
			interval_computing_pressure_relaxation += tick_count::now() - time_instance;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0)
					restart_io.writeToFile(number_of_iterations);
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			time_instance = tick_count::now();
			
			water_block.updateCellLinkedList();
			water_air_complex.updateConfiguration();
			water_wall_complex.updateConfiguration();

			air_block.updateCellLinkedList();
			air_water_complex.updateConfiguration();
			air_wall_complex.updateConfiguration();

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
