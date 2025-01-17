/**
 * @file 	two_phase_dambreak.cpp
 * @brief 	2D two-phase dambreak flow.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for multi-phase simulation.
 * @author 	Chi Zhang, Xiangyu Hu and Boyang Li
 */
#include "SPHconfig.h"
#include "phase_transition_solid.h"
#include "case.h"
#include "photon.h"
#include "photon_particles.h"
#include "photon_dynamics.h"
using vdWSolidPhaseTransitionDynamics=
	phaseTransitionDynamics<	
							WallParticles,
							vdWParticles
	>;
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
	Real diffusion_coff;
};
struct SolidPar
{
	Real T0;
	Real diffusion_coff;
	Real rho0;
	Real stiffness;
};

struct execClass{
	bool doParallel;
	template<class execType>
	auto doexec(execType &target, Real dt = 0){
		if(doParallel) return target.parallel_exec(dt);
		return target.exec(dt);
	}
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
	execClass exec;
	Real adaptation = 1.3;
	bool writePhoton = 0;
	(*cfg.vdWFluids)[0].lookupValue("viscosity",air.mu);
	(*cfg.vdWFluids)[0].lookupValue("rho_0",air.rho0);
	(*cfg.vdWFluids)[0].lookupValue("a",air.alpha);
	(*cfg.vdWFluids)[0].lookupValue("rho_m",air.rho_m);
	(*cfg.vdWFluids)[0].lookupValue("gamma",air.gamma);
	(*cfg.vdWFluids)[0].lookupValue("temperature",air.T0);
	(*cfg.vdWFluids)[0].lookupValue("diffusion_coff",air.diffusion_coff);
	(*cfg.vdWFluids)[1].lookupValue("viscosity",water.mu);
	(*cfg.vdWFluids)[1].lookupValue("rho_0",water.rho0);
	(*cfg.vdWFluids)[1].lookupValue("a",water.alpha);
	(*cfg.vdWFluids)[1].lookupValue("rho_m",water.rho_m);
	(*cfg.vdWFluids)[1].lookupValue("diffusion_coff",water.diffusion_coff);
	(*cfg.vdWFluids)[1].lookupValue("gamma",water.gamma);
	(*cfg.vdWFluids)[1].lookupValue("temperature",water.T0);
	(*cfg.Solids)[0].lookupValue("temperature",wall.T0);
	(*cfg.Solids)[0].lookupValue("diffusion_coff",wall.diffusion_coff);
	(*cfg.Solids)[0].lookupValue("rho_0",wall.rho0);
	(*cfg.Solids)[0].lookupValue("stiffness",wall.stiffness);
	cfg.Job->lookupValue("particle_spacing_ref",particle_spacing_ref);
	cfg.Job->lookupValue("doParallel",exec.doParallel);
	cfg.Job->lookupValue("adaptation",adaptation);
	cfg.Job->lookupValue("writePhoton",writePhoton);
	Real BW = particle_spacing_ref * 4 * adaptation;
	BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
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
	printf("%f, %f, %f, %f, %f, %f\n", water.rho0, water.rho_m, water.gamma, water.alpha, water.molmass, water.mu, water.diffusion_coff);

	auto adp = makeShared<SPHAdaptation>(adaptation, 1);
	WaterBlock water_block(sph_system, "WaterBody", adp);
	vdWParticles diffusion_fluid_body_particles(water_block, makeShared<FluidMaterial>(water.rho0, water.rho_m, water.gamma, water.alpha, water.molmass, water.mu, water.diffusion_coff));
	diffusion_fluid_body_particles.addAVariableToWrite<indexScalar, Real>("Density");

	AirBlock air_block(sph_system, "AirBody", adp);
	vdWParticles diffusion_air_body_particles(air_block, makeShared<FluidMaterial>(air.rho0, air.rho_m, air.gamma, air.alpha, air.molmass, air.mu, air.diffusion_coff));
	diffusion_air_body_particles.addAVariableToWrite<indexScalar, Real>("Density");

	WallBoundary wall_boundary(sph_system, "Wall", BW, adp);
	WallParticles
		wall_particles(wall_boundary, makeShared<WallMaterial>(wall.diffusion_coff, wall.rho0, wall.stiffness));

	auto photon_adp = makeShared<SPHAdaptation>(1.5, 1);
	PhotonBlock photon_block(sph_system, "PhotonBody", adp);
	PhotonParticles photon_particles(photon_block, makeShared<Photon>(&(*cfg.Photons)[0]));
	photon_particles.addAVariableToWrite<indexScalar, Real>("Intensity");

	size_t cacheSize = diffusion_fluid_body_particles.total_real_particles_
	+ wall_particles.total_real_particles_ + diffusion_air_body_particles.total_real_particles_;
	StdLargeVec<size_t> caches1;
	caches1.resize(cacheSize);
	StdLargeVec<size_t> caches2;
	caches2.resize(cacheSize);
	StdLargeVec<size_t> caches3;
	caches3.resize(cacheSize);
	StdLargeVec<size_t> caches4;
	caches4.resize(cacheSize);
	StdLargeVec<size_t> caches5;
	caches5.resize(cacheSize);


	Real laserIntensity = 0;
	(*cfg.Photons)[0].lookupValue("Intensity",laserIntensity);
	PlainWave laser(Vecd(0,-1),laserIntensity);
	PhotonInitialization photon_init(photon_block, laser, 0);
	Real photonExist = 0;
	//exec.doexec(photon_init,0);
	Real lightspeed = 0;
	(*cfg.Photons)[0].lookupValue("c0",lightspeed);
	Real lightdt = 0.1 * adp->ReferenceSmoothingLength() / lightspeed;
	
	BodyRelationContact photon_contact(photon_block,{&wall_boundary, &water_block});
	PhotonPropagation photon_propagation(photon_contact);

	vdWSolidPhaseTransitionDynamics phaseTransition(wall_particles, diffusion_fluid_body_particles, cfg.PhaseTransition, caches1, caches2);
	ObserverBody temperature_observer(sph_system, "FluidObserver");
	ObserverParticles temperature_observer_particles(temperature_observer, makeShared<ObserverParticleGenerator>());
	/**
	 * @brief Material property, partilces and body creation of air.
	 */
	/** topology */
	StdVec<Real> lengths = {1};
	(*cfg.vdWFluids)[1].lookupValue("attr_length", lengths[0]);
	BodyRelationInnerMultiLength air_body_inner(air_block,lengths);
	BodyRelationInnerMultiLength fluid_body_inner(water_block,lengths);
	BodyRelationInner solid_body_inner(wall_boundary);
	BodyRelationContact water_air_contact(water_block,{&air_block});
	BodyRelationContact air_water_contact(air_block,{&water_block});
	BodyRelationContactMultiLength water_wall_contact(water_block,{&wall_boundary},lengths);
	BodyRelationContactMultiLength air_wall_contact(air_block,{&wall_boundary},lengths);
	BodyRelationContact wall_contact(wall_boundary,{&water_block,&air_block});
	//ComplexBodyRelation water_wall_complex(fluid_body_inner, { &wall_boundary });
	ComplexBodyRelation water_air_complex(fluid_body_inner, water_air_contact);
	ComplexBodyRelation water_wall_complex(fluid_body_inner, water_wall_contact);
	ComplexBodyRelation air_water_complex(air_body_inner, air_water_contact);
	ComplexBodyRelation air_wall_complex(air_body_inner, air_wall_contact);
	ComplexBodyRelation wall_complex(solid_body_inner, wall_contact);
	//BaseBodyRelationContact &water_wall_contact = water_wall_complex.contact_relation_;
	BodyRelationContact fluid_observer_contact(temperature_observer, {&water_block});

	/**
	 * @brief 	Define all numerical methods which are used in this case.
	 */
	 /** Define external force. */

	cfg.ExternalForce->lookupValue("gravity",gravity_g);
	Gravity		gravity(Vecd(0.0, -gravity_g));
	HeatSource heat_source;
	/**
	 * @brief 	Methods used for time stepping.
	 */
	 /** Initialize particle acceleration. */
	ThermosolidBodyInitialCondition thermosolid_condition(wall_boundary, wall.T0, BW);
	ThermofluidBodyInitialCondition thermofluid_initial_condition(water_block, water.T0);
	ThermofluidBodyInitialCondition thermoair_initial_condition(air_block, air.T0);

	TimeStepInitialization		initialize_a_water_step(water_block, gravity);
	TimeStepInitialization          initialize_a_air_step(air_block, gravity);

	DiffusionSourceInitialization<FluidBody, FluidParticles, vdWFluid> initialize_a_water_step_thermo(water_block, heat_source, 0);
	DiffusionSourceInitialization<FluidBody, FluidParticles, vdWFluid> initialize_a_air_step_thermo(air_block, heat_source, 0);
	//StressTensorHeatSource<FluidBody, FluidParticles, vdWFluid>  stress_tensor_heat_water(fluid_body_inner, 0);
	//vdWAttractionHeatSource<FluidBody, FluidParticles, vdWFluid>  vdW_attr_heat_water(water_block, 0);
	/** Corrected strong configuration for diffusion solid body. */
	solid_dynamics::CorrectConfiguration 			correct_configuration(solid_body_inner);
	/** Time step size calculation. */
	GetDiffusionTimeStepSize<FluidBody, FluidParticles, vdWFluid> get_thermal_time_step(water_block);
	/** Diffusion process between three diffusion bodies. */
	//ThermalRelaxationComplexWA 	thermal_relaxation_complex_wa(water_air_complex);
	ThermalRelaxationComplex        thermal_relaxation_complex_water(water_wall_complex);//, water_wall_contact);
	ThermalRelaxationComplex        thermal_relaxation_complex_air(air_wall_complex);//, air_wall_contact);
	ThermalRelaxationComplexWall    thermal_relaxation_complex_wall(wall_complex);
	/**
	 * @brief 	Algorithms of fluid dynamics.
	 */
	 /** Evaluation of density by summation approach. */
	//fluid_dynamics::DensitySummationWithMergingAndSplitting 	
	//	update_water_density_by_summation(fluid_body_inner, water_wall_contact, caches1, caches2, caches3, caches4, caches5);
	fluid_dynamics::DensitySummationFreeSurfaceComplexWithoutUpdate 	
		update_water_density_by_summation(fluid_body_inner, water_wall_contact);
	fluid_dynamics::DensitySummationFreeSurfaceComplexWithoutUpdate 	
		update_air_density_by_summation(air_body_inner, air_wall_contact);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize 	get_water_advection_time_step_size(water_block, U_max);
	fluid_dynamics::AdvectionTimeStepSize 	get_air_advection_time_step_size(air_block, U_max);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_water_time_step_size(water_block);
	fluid_dynamics::AcousticTimeStepSize get_air_time_step_size(air_block);
	/** Pressure relaxation for water by using position verlet time stepping. */
	fluid_dynamics::vdWPressureRelaxationRiemannWithWall 
		water_pressure_relaxation(water_air_complex, water_wall_contact);
	fluid_dynamics::vdWPressureRelaxationRiemannWithWall 
		air_pressure_relaxation(air_water_complex, air_wall_contact);
	fluid_dynamics::vdWDensityRelaxationRiemannWithWall 
		water_density_relaxation(fluid_body_inner, water_wall_contact);
	fluid_dynamics::vdWDensityRelaxationRiemannWithWall 
		air_density_relaxation(air_body_inner, air_wall_contact);
	//fluid_dynamics::ViscousAccelerationMultiPhase
	//fluid_dynamics::vdWViscousAccelerationInner
	//	water_viscou_acceleration(fluid_body_inner);
	/** Suface tension and wetting effects. */
//	fluid_dynamics::FreeSurfaceIndicationComplex
//		surface_detection(fluid_body_inner, water_wall_contact);
//	fluid_dynamics::ColorFunctionGradientComplex
//		color_gradient(fluid_body_inner, water_wall_contact);
//	fluid_dynamics::ColorFunctionGradientInterplationInner
//		color_gradient_interpolation(fluid_body_inner);
	//fluid_dynamics::SurfaceTensionAccelerationInner
	//	surface_tension_acceleration(fluid_body_inner, tension_force);
	/** Wetting effects. */
	//fluid_dynamics::SurfaceNormWithWall
	//	wetting_norm(water_wall_contact, contact_angle);
	/** Computing vorticity in the flow. */
	//fluid_dynamics::VorticityInner 	compute_vorticity(fluid_body_inner);
	/** Output the body states. */
	BodyStatesRecordingToVtp 		body_states_recording(in_output, sph_system.real_bodies_);
	RegressionTestEnsembleAveraged<ObservedQuantityRecording<indexScalar, Real>>
		write_fluid_phi("Temperature", in_output, fluid_observer_contact);
	RegressionTestEnsembleAveraged<ObservedQuantityRecording<indexScalar, Real>>
		write_fluid_p("Pressure", in_output, fluid_observer_contact);
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
	exec.doexec(correct_configuration);
	exec.doexec(thermosolid_condition);
	exec.doexec(thermofluid_initial_condition);
	exec.doexec(thermoair_initial_condition);
	Real dt_thermal = exec.doexec(get_thermal_time_step);
	 /** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (sph_system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.restart_step_);

		//air_water_complex.updateConfiguration();
		//air_wall_contact.updateConfiguration();
	}		
	exec.doexec(phaseTransition);
	water_block.updateCellLinkedList();
	air_block.updateCellLinkedList();
	wall_boundary.updateCellLinkedList();
	water_air_complex.updateConfiguration();
	water_wall_contact.updateConfiguration();
	air_wall_complex.updateConfiguration();
	air_water_contact.updateConfiguration();
	wall_complex.updateConfiguration();
	/** Output the start states of bodies. */
	body_states_recording.writeToFile(0);
	/**
	 * @brief 	Basic parameters.
	 */
	int denseWriting = 0;
	cfg.Job->lookupValue("denseWriting",denseWriting);
	size_t number_of_iterations = sph_system.restart_step_;
	int screen_output_interval = 30;
	cfg.Job->lookupValue("screenOut",screen_output_interval);
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 1000; 	/**< End time. */
	Real D_Time = 1;		/**< Time stamps for output of body states. */
	Real Dt = 0.0;			/**< Default advection time step sizes. */
	Real dt = 0.0; 			/**< Default acoustic time step sizes. */
	cfg.Job->lookupValue("EndTime",End_Time);
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	tick_count::interval_t interval_computing_time_step;
	tick_count::interval_t interval_computing_pressure_relaxation;
	tick_count::interval_t interval_updating_configuration;
	tick_count time_instance;
	bool firstiter = 1;
	int debugframe = 99999;
	Real gent = 20;
	cfg.Job->lookupValue("debugFrame",debugframe);
	/**
	 * @brief 	Main loop starts here.
	 */

	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (1)
		{
			time_instance = tick_count::now();
			exec.doexec(update_water_density_by_summation,Dt);
			exec.doexec(update_air_density_by_summation,Dt);
			exec.doexec(phaseTransition);
			water_block.updateCellLinkedList();
			air_block.updateCellLinkedList();
			wall_boundary.updateCellLinkedList();
			water_air_complex.updateConfiguration();
			water_wall_contact.updateConfiguration();
			air_wall_complex.updateConfiguration();
			air_water_contact.updateConfiguration();
			wall_complex.updateConfiguration();

			interval_computing_pressure_relaxation += tick_count::now() - time_instance;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";

				if(!denseWriting) {
					photon_block.setNewlyUpdated();
					body_states_recording.writeToFile();
					debugframe--;
				}
				if (number_of_iterations % restart_output_interval == 0)
					restart_io.writeToFile(number_of_iterations);
			}
			number_of_iterations++;
			if(integration_time >= D_Time) break;
			/** Acceleration due to viscous force and gravity. */
			time_instance = tick_count::now();
			//exec.doexec(initialize_a_air_step);

			Real Dt_f = exec.doexec(get_water_advection_time_step_size);
			Real Dt_a = exec.doexec(get_air_advection_time_step_size);
			Dt = SMIN(Dt_f, Dt_a);


			//exec.doexec(air_transport_correction,Dt);

			//exec.doexec(air_viscou_acceleration);

			//exec.doexec(surface_detection);
			//exec.doexec(color_gradient);
			//exec.doexec(color_gradient_interpolation);
			//exec.doexec(wetting_norm);
			//exec.doexec(surface_tension_acceleration);

			interval_computing_time_step += tick_count::now() - time_instance;

			//exec.doexec(initialize_a_air_step_thermo);
			/** Dynamics including pressure relaxation. */
			time_instance = tick_count::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				if(debugframe==0) {
					int halt = 1;
				}
				Real dt_f = exec.doexec(get_water_time_step_size);
				Real dt_a = exec.doexec(get_air_time_step_size);
				dt = SMIN(SMIN(dt_f, Dt),dt_a)/2;
				exec.doexec(initialize_a_water_step_thermo);
				exec.doexec(initialize_a_air_step_thermo);
				exec.doexec(initialize_a_water_step);
				//exec.doexec(water_viscou_acceleration);
				exec.doexec(air_pressure_relaxation,dt);
				exec.doexec(water_pressure_relaxation,dt);
				//exec.doexec(air_pressure_relaxation,dt);
				//exec.doexec(stress_tensor_heat_water);
				exec.doexec(water_density_relaxation,dt);
				exec.doexec(air_density_relaxation,dt);
				//exec.doexec(air_density_relaxation,dt);
				int steps = 1+dt/lightdt;
				//printf("light steps: %d\n",steps);
				if(photonExist){
					auto t = GlobalStaticVariables::physical_time_;
					for(int i=0; i < steps;i++){
						if(!photon_particles.total_real_particles_) {
							photonExist = 0;
							break;
						}
						photon_block.updateCellLinkedList();
						photon_contact.updateConfiguration();
						exec.doexec(photon_propagation,dt/steps);
						if(writePhoton) {
							GlobalStaticVariables::physical_time_ += dt/steps;
							water_block.setNewlyUpdated();
							wall_boundary.setNewlyUpdated();
							air_block.setNewlyUpdated();
							body_states_recording.writeToFile();
						}
					}
					if(writePhoton) {
						GlobalStaticVariables::physical_time_ = t;
					}
				}else{
					if(GlobalStaticVariables::physical_time_-gent>0){
						exec.doexec(photon_init);
						photonExist = 1;
						gent+=20;
					}
				}
				//exec.doexec(thermal_relaxation_complex_wa,dt);
				//exec.doexec(vdW_attr_heat_water);
				if(dt_thermal < dt){
					int steps = 1+dt/dt_thermal;
					for(int i = 0; i < steps;i++){
						exec.doexec(thermal_relaxation_complex_air,dt/steps);
						exec.doexec(thermal_relaxation_complex_wall,dt/steps);
						exec.doexec(thermal_relaxation_complex_water,dt/steps);
					}
				}else{
					exec.doexec(thermal_relaxation_complex_water,dt);
					exec.doexec(thermal_relaxation_complex_air,dt);
					exec.doexec(thermal_relaxation_complex_wall,dt);
				}
				//exec.doexec(thermal_relaxation_complex_aw,dt);
				//exec.doexec(stress_tensor_heat_air);
				//exec.doexec(vdW_attr_heat_air);
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				//if(firstiter) {
				//	body_states_recording.writeToFile();
				//	firstiter =0;
				//}
				if(denseWriting) {
					photon_block.setNewlyUpdated();
					body_states_recording.writeToFile();
					debugframe--;
				}
			}

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
