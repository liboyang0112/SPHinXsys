/**
 * @file 	excitation-contraction.cpp
 * @brief 	This is the case studying the electromechanics on a biventricular heart model in 3D.
 * @author 	Chi Zhang and Xiangyu Hu
 * 			Unit :
 *			time t = ms = 12.9 [-]
 * 			length l = mm
 * 			mass m = g
 *			density rho = g * (mm)^(-3)
 *			Pressure pa = g * (mm)^(-1) * (ms)^(-2)
 *			diffusion d = (mm)^(2) * (ms)^(-2) 
 */
/** 
 * SPHinXsys Library. 
 */
#include "sphinxsys.h"
/** Namespace cite here. */
using namespace SPH;
/** Geometry parameter. */
/** Set the file path to the stl file. */
std::string full_path_to_stl_file = "./input/heart-new.stl";
Real length_scale = 1.0;
Real time_scale = 1.0 / 12.9;
Real stress_scale = 1.0e-6;
/** Parameters and physical properties. */
Vec3d domain_lower_bound(-55.0 * length_scale, -75.0 * length_scale, -35.0 * length_scale);
Vec3d domain_upper_bound(35.0 * length_scale, 5.0 * length_scale, 35.0 * length_scale);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 45.0; /**< Initial particle spacing. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);

/** Material properties. */
Real rho0_s = 1.06e-3;
/** Active stress factor */
Real k_a = 100 * stress_scale;
Real a0[4] = {496.0 * stress_scale, 15196.0 * stress_scale, 3283.0 * stress_scale, 662.0 * stress_scale};
Real b0[4] = {7.209, 20.417, 11.176, 9.466};
/** reference stress to achieve weakly compressible condition */
Real poisson = 0.4995;
Real bulk_modulus = 2.0 * a0[0] * (1.0 + poisson) / (3.0 * (1.0 - 2.0 * poisson));
/** Electrophysiology parameters. */
Real diffusion_coff = 0.8;
Real bias_coff = 0.0;
/** Electrophysiology parameters. */
Real c_m = 1.0;
Real k = 8.0;
Real a = 0.01;
Real b = 0.15;
Real mu_1 = 0.2;
Real mu_2 = 0.3;
Real epsilon = 0.002;
/** Fibers and sheet. */
Vec3d fiber_direction(1.0, 0.0, 0.0);
Vec3d sheet_direction(0.0, 1.0, 0.0);

/** 
 * Define geometry and initial conditions of SPH bodies. 
 */
class HeartBody : public SolidBody
{
public:
	HeartBody(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		Vecd translation(-53.5 * length_scale, -70.0 * length_scale, -32.5 * length_scale);
		TriangleMeshShapeSTL triangle_mesh_heart_shape(full_path_to_stl_file, translation, length_scale);
		body_shape_.add<LevelSetShape>(this, triangle_mesh_heart_shape);
	}
};
/** Set diffusion relaxation method. */
class DiffusionRelaxation
	: public RelaxationOfAllDiffusionSpeciesRK2<
		  SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle,
		  RelaxationOfAllDiffussionSpeciesInner<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>,
		  BodyRelationInner>
{
public:
	explicit DiffusionRelaxation(BodyRelationInner &body_inner_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_inner_relation){};
	virtual ~DiffusionRelaxation(){};
};
/** Imposing diffusion boundary condition */
class DiffusionBCs
	: public ConstrainDiffusionBodyRegion<SolidBody, ElasticSolidParticles, BodySurface, LocallyOrthotropicMuscle>
{
protected:
	size_t phi_;
	virtual void Update(size_t index_i, Real dt = 0.0) override
	{
		Vecd dist_2_face = body_->body_shape_.findNormalDirection(pos_n_[index_i]);
		Vecd face_norm = dist_2_face / (dist_2_face.norm() + 1.0e-15);

		Vecd center_norm = pos_n_[index_i] / (pos_n_[index_i].norm() + 1.0e-15);

		Real angle = dot(face_norm, center_norm);
		if (angle >= 0.0)
		{
			species_n_[phi_][index_i] = 1.0;
		}
		else
		{
			if (pos_n_[index_i][1] < -body_->sph_adaptation_->ReferenceSpacing())
				species_n_[phi_][index_i] = 0.0;
		}
	};

public:
	DiffusionBCs(SolidBody &body, BodySurface &body_part)
		: ConstrainDiffusionBodyRegion<SolidBody, ElasticSolidParticles, BodySurface, LocallyOrthotropicMuscle>(body, body_part)
	{
		phi_ = material_->SpeciesIndexMap()["Phi"];
	};
	virtual ~DiffusionBCs(){};
};
/** Compute Fiber and Sheet direction after diffusion */
class ComputeFiberandSheetDirections
	: public DiffusionBasedMapping<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>
{
protected:
	size_t phi_;
	Real beta_epi_, beta_endo_;
	/** We define the centerline vector, which is parallel to the ventricular centerline and pointing  apex-to-base.*/
	Vecd center_line_;
	virtual void Update(size_t index_i, Real dt = 0.0) override
	{
		/**
		 * Ref: original doi.org/10.1016/j.euromechsol.2013.10.009
		 * 		Present  doi.org/10.1016/j.cma.2016.05.031
		 */
		/** Probe the face norm from Levelset field. */
		Vecd dist_2_face = body_->body_shape_.findNormalDirection(pos_n_[index_i]);
		Vecd face_norm = dist_2_face / (dist_2_face.norm() + 1.0e-15);
		Vecd center_norm = pos_n_[index_i] / (pos_n_[index_i].norm() + 1.0e-15);
		if (dot(face_norm, center_norm) <= 0.0)
		{
			face_norm = -face_norm;
		}
		/** Compute the centerline's projection on the plane orthogonal to face norm. */
		Vecd circumferential_direction = SimTK::cross(center_line_, face_norm);
		Vecd cd_norm = circumferential_direction / (circumferential_direction.norm() + 1.0e-15);
		/** The rotation angle is given by beta = (beta_epi - beta_endo) phi + beta_endo */
		Real beta = (beta_epi_ - beta_endo_) * species_n_[phi_][index_i] + beta_endo_;
		/** Compute the rotation matrix through Rodrigues rotation formulation. */
		Vecd f_0 = cos(beta) * cd_norm + sin(beta) * SimTK::cross(face_norm, cd_norm) +
				   dot(face_norm, cd_norm) * (1.0 - cos(beta)) * face_norm;

		if (pos_n_[index_i][1] < -body_->sph_adaptation_->ReferenceSpacing())
		{
			material_->local_f0_[index_i] = f_0 / (f_0.norm() + 1.0e-15);
			material_->local_s0_[index_i] = face_norm;
		}
		else
		{
			material_->local_f0_[index_i] = Vecd(0);
			material_->local_s0_[index_i] = Vecd(0);
		}
	};

public:
	explicit ComputeFiberandSheetDirections(SolidBody &body)
		: DiffusionBasedMapping<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>(body)
	{
		phi_ = material_->SpeciesIndexMap()["Phi"];
		center_line_ = Vecd(0.0, 1.0, 0.0);
		beta_epi_ = -(70.0 / 180.0) * M_PI;
		beta_endo_ = (80.0 / 180.0) * M_PI;
	};
	virtual ~ComputeFiberandSheetDirections(){};
};
//	define shape parameters which will be used for the constrained body part.
class MuscleBaseShapeParameters : public TriangleMeshShapeBrick::ShapeParameters
{
public:
	MuscleBaseShapeParameters() : TriangleMeshShapeBrick::ShapeParameters()
	{
		Real l = domain_upper_bound[0] - domain_lower_bound[0];
		Real w = domain_upper_bound[2] - domain_lower_bound[2];
		halfsize_ = Vec3d(0.5 * l, 1.0 * dp_0, 0.5 * w);
		resolution_ = 20;
		translation_ = Vec3d(-10.0 * length_scale, -1.0 * dp_0, 0.0);
	}
};
//	application dependent initial condition
class ApplyStimulusCurrentSI
	: public electro_physiology::ElectroPhysiologyInitialCondition
{
protected:
	size_t voltage_;

	void Update(size_t index_i, Real dt) override
	{
		if (-30.0 * length_scale <= pos_n_[index_i][0] && pos_n_[index_i][0] <= -15.0 * length_scale)
		{
			if (-2.0 * length_scale <= pos_n_[index_i][1] && pos_n_[index_i][1] <= 0.0)
			{
				if (-3.0 * length_scale <= pos_n_[index_i][2] && pos_n_[index_i][2] <= 3.0 * length_scale)
				{
					species_n_[voltage_][index_i] = 0.92;
				}
			}
		}
	};

public:
	explicit ApplyStimulusCurrentSI(SolidBody &muscle)
		: electro_physiology::ElectroPhysiologyInitialCondition(muscle)
	{
		voltage_ = material_->SpeciesIndexMap()["Voltage"];
	};
};
/**
 * application dependent initial condition 
 */
class ApplyStimulusCurrentSII
	: public electro_physiology::ElectroPhysiologyInitialCondition
{
protected:
	size_t voltage_;

	void Update(size_t index_i, Real dt) override
	{
		if (0.0 <= pos_n_[index_i][0] && pos_n_[index_i][0] <= 6.0 * length_scale)
		{
			if (-6.0 * length_scale <= pos_n_[index_i][1])
			{
				if (12.0 * length_scale <= pos_n_[index_i][2])
				{
					species_n_[voltage_][index_i] = 0.95;
				}
			}
		}
	};

public:
	explicit ApplyStimulusCurrentSII(SolidBody &muscle)
		: electro_physiology::ElectroPhysiologyInitialCondition(muscle)
	{
		voltage_ = material_->SpeciesIndexMap()["Voltage"];
	};
};
/**
 * define observer particle generator.
 */
class ObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	ObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		/** position and volume. */
		positions_volumes_.push_back(std::make_pair(Vecd(-45.0 * length_scale, -30.0 * length_scale, 0.0), 0.0));
		positions_volumes_.push_back(std::make_pair(Vecd(0.0, -30.0 * length_scale, 26.0 * length_scale), 0.0));
		positions_volumes_.push_back(std::make_pair(Vecd(-30.0 * length_scale, -50.0 * length_scale, 0.0), 0.0));
		positions_volumes_.push_back(std::make_pair(Vecd(0.0, -50.0 * length_scale, 20.0 * length_scale), 0.0));
		positions_volumes_.push_back(std::make_pair(Vecd(0.0, -70.0 * length_scale, 0.0), 0.0));
	}
};
/** 
 * The main program. 
 */
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	SPHSystem section
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, dp_0);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = false;
	/** Tag for reload initially relaxed particles. */
	system.reload_particles_ = true;
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;
//handle command line arguments
#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(ac, av);
#endif
	/** in- and output environment. */
	In_Output in_output(system);
	//----------------------------------------------------------------------
	//	SPHData section
	//----------------------------------------------------------------------
	/** create a SPH body, material and particles */
	HeartBody physiology_body(system, "ExcitationHeart");
	SharedPtr<ParticleGenerator> physiology_body_particle_generator = makeShared<ParticleGeneratorLattice>();
	if (!system.run_particle_relaxation_ && system.reload_particles_)
		physiology_body_particle_generator = makeShared<ParticleGeneratorReload>(in_output, physiology_body.getBodyName());
	AlievPanfilowModel muscle_reaction_model(k_a, c_m, k, a, b, mu_1, mu_2, epsilon);
	SharedPtr<LocalMonoFieldElectroPhysiology> myocardium_excitation =
		makeShared<LocalMonoFieldElectroPhysiology>(muscle_reaction_model, diffusion_coff, bias_coff, fiber_direction);
	ElectroPhysiologyParticles physiology_articles(physiology_body, myocardium_excitation, physiology_body_particle_generator);

	/** create a SPH body, material and particles */
	HeartBody mechanics_body(system, "ContractionHeart");
	SharedPtr<ParticleGenerator> mechanics_body_particle_generator = makeShared<ParticleGeneratorLattice>();
	if (!system.run_particle_relaxation_ && system.reload_particles_)
		mechanics_body_particle_generator = makeShared<ParticleGeneratorReload>(in_output, mechanics_body.getBodyName());
	SharedPtr<ActiveMuscle<LocallyOrthotropicMuscle>> myocardium_muscle =
		makeShared<ActiveMuscle<LocallyOrthotropicMuscle>>(rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0);
	ActiveMuscleParticles mechanics_particles(mechanics_body, myocardium_muscle, mechanics_body_particle_generator);

	/** check whether reload material properties. */
	if (!system.run_particle_relaxation_ && system.reload_particles_)
	{
		ReloadMaterialParameterIO read_muscle_fiber_and_sheet(in_output, myocardium_muscle);
		ReloadMaterialParameterIO read_myocardium_excitation_fiber(in_output, myocardium_excitation, myocardium_muscle->LocalParametersName());
		read_muscle_fiber_and_sheet.readFromFile();
		read_myocardium_excitation_fiber.readFromFile();
	}
	//----------------------------------------------------------------------
	//	SPH Particle relaxation section
	//----------------------------------------------------------------------
	/** check whether run particle relaxation for body fitted particle distribution. */
	if (system.run_particle_relaxation_)
	{
		HeartBody relax_body(system, "RelaxationHeart");
		StdVec<std::string> species_name_list{"Phi"};
		SharedPtr<DiffusionReaction<ElasticSolidParticles, LocallyOrthotropicMuscle>> relax_body_material =
			makeShared<DiffusionReaction<ElasticSolidParticles, LocallyOrthotropicMuscle>>(
				species_name_list, rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0);
		relax_body_material->initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", diffusion_coff);
		DiffusionReactionParticles<ElasticSolidParticles, LocallyOrthotropicMuscle> diffusion_particles(relax_body, relax_body_material);

		/** topology */
		BodyRelationInner relax_body_inner(relax_body);
		/** Random reset the relax solid particle position. */
		RandomizePartilePosition random_particles(relax_body);

		/**
		 * @brief 	Algorithms for particle relaxation.
		 */
		/** A  Physics relaxation step. */
		relax_dynamics::RelaxationStepInner relaxation_step_inner(relax_body_inner);
		/**
		 * Diffusion process.
		 */
		/** Time step for diffusion. */
		GetDiffusionTimeStepSize<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle> get_time_step_size(relax_body);
		/** Diffusion process for diffusion body. */
		DiffusionRelaxation diffusion_relaxation(relax_body_inner);
		/** Compute the fiber and sheet after diffusion. */
		ComputeFiberandSheetDirections compute_fiber_sheet(relax_body);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_relax_body_state_to_vtp(in_output, {relax_body});
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(in_output, {&relax_body, &relax_body}, {physiology_body.getBodyName(), mechanics_body.getBodyName()});
		/** Write material property to xml file. */
		ReloadMaterialParameterIO write_material_property(in_output, relax_body_material, myocardium_muscle->LocalParametersName());
		/**
		 * @brief 	Physics relaxation starts here.
		 */
		/** Relax the elastic structure. */
		random_particles.parallel_exec(0.25);
		relaxation_step_inner.surface_bounding_.parallel_exec();
		write_relax_body_state_to_vtp.writeToFile(0.0);
		/**
		 * From here the time stepping begins.
		 * Set the starting time.
		 */
		int ite = 0;
		int relax_step = 1000;
		int diffusion_step = 100;
		while (ite < relax_step)
		{
			relaxation_step_inner.parallel_exec();
			ite++;
			if (ite % 100 == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
				write_relax_body_state_to_vtp.writeToFile(Real(ite) * 1.0e-4);
			}
		}

		BodySurface surface_part(relax_body);
		/** constraint boundary condition for diffusion. */
		DiffusionBCs impose_diffusion_bc(relax_body, surface_part);
		impose_diffusion_bc.parallel_exec();

		write_relax_body_state_to_vtp.writeToFile(Real(ite) * 1.0e-4);

		Real dt = get_time_step_size.parallel_exec();
		while (ite <= diffusion_step + relax_step)
		{
			diffusion_relaxation.parallel_exec(dt);
			impose_diffusion_bc.parallel_exec();
			if (ite % 10 == 0)
			{
				std::cout << "Diffusion steps N=" << ite - relax_step << "	dt: " << dt << "\n";
				write_relax_body_state_to_vtp.writeToFile(Real(ite) * 1.0e-4);
			}
			ite++;
		}
		compute_fiber_sheet.exec();
		ite++;
		write_relax_body_state_to_vtp.writeToFile(Real(ite) * 1.0e-4);
		compute_fiber_sheet.parallel_exec();
		write_material_property.writeToFile(0);
		write_particle_reload_files.writeToFile(0);

		return 0;
	}
	//----------------------------------------------------------------------
	//	SPH Observation section
	//----------------------------------------------------------------------
	ObserverBody voltage_observer(system, "VoltageObserver");
	ObserverParticles observer_particles(voltage_observer, makeShared<ObserverParticleGenerator>());
	/** Define muscle Observer. */
	ObserverBody myocardium_observer(system, "MyocardiumObserver");
	ObserverParticles disp_observer_particles(myocardium_observer, makeShared<ObserverParticleGenerator>());
	//----------------------------------------------------------------------
	//	SPHBody relation (topology) section
	//----------------------------------------------------------------------
	BodyRelationInner physiology_body_inner(physiology_body);
	BodyRelationInner mechanics_body_inner(mechanics_body);
	BodyRelationContact physiology_body_contact(physiology_body, {&mechanics_body});
	BodyRelationContact mechanics_body_contact(mechanics_body, {&physiology_body});
	BodyRelationContact voltage_observer_contact(voltage_observer, {&physiology_body});
	BodyRelationContact myocardium_observer_contact(myocardium_observer, {&mechanics_body});
	//----------------------------------------------------------------------
	//	SPH Method section
	//----------------------------------------------------------------------
	// Corrected configuration.
	solid_dynamics::CorrectConfiguration correct_configuration_excitation(physiology_body_inner);
	// Time step size calculation.
	electro_physiology::GetElectroPhysiologyTimeStepSize get_physiology_time_step(physiology_body);
	// Diffusion process for diffusion body.
	electro_physiology::ElectroPhysiologyDiffusionRelaxationInner diffusion_relaxation(physiology_body_inner);
	// Solvers for ODE system.
	electro_physiology::ElectroPhysiologyReactionRelaxationForward reaction_relaxation_forward(physiology_body);
	electro_physiology::ElectroPhysiologyReactionRelaxationBackward reaction_relaxation_backward(physiology_body);
	//	Apply the Iron stimulus.
	ApplyStimulusCurrentSI apply_stimulus_s1(physiology_body);
	ApplyStimulusCurrentSII apply_stimulus_s2(physiology_body);
	// Active mechanics.
	solid_dynamics::CorrectConfiguration correct_configuration_contraction(mechanics_body_inner);
	observer_dynamics::CorrectInterpolationKernelWeights correct_kernel_weights_for_interpolation(mechanics_body_contact);
	/** Interpolate the active contract stress from electrophysiology body. */
	observer_dynamics::InterpolatingAQuantity<indexScalar, Real> active_stress_interpolation(mechanics_body_contact, "ActiveContractionStress", "ActiveContractionStress");
	/** Interpolate the particle position in physiology_body  from mechanics_body. */
	observer_dynamics::InterpolatingAQuantity<indexVector, Vecd> interpolation_particle_position(physiology_body_contact, "Position", "Position");
	/** Time step size calculation. */
	solid_dynamics::AcousticTimeStepSize get_mechanics_time_step(mechanics_body);
	/** active and passive stress relaxation. */
	solid_dynamics::StressRelaxationFirstHalf stress_relaxation_first_half(mechanics_body_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(mechanics_body_inner);
	/** Constrain region of the inserted body. */
	MuscleBaseShapeParameters muscle_base_parameters;
	TriangleMeshShapeBrick muscle_base_shape(muscle_base_parameters);
	BodyRegionByParticle muscle_base(mechanics_body, "Holder", muscle_base_shape);
	solid_dynamics::ConstrainSolidBodyRegion constrain_holder(mechanics_body, muscle_base);
	//----------------------------------------------------------------------
	//	SPH Output section
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_states(in_output, system.real_bodies_);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<indexScalar, Real>>
		write_voltage("Voltage", in_output, voltage_observer_contact);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<indexVector, Vecd>>
		write_displacement("Position", in_output, myocardium_observer_contact);
	//----------------------------------------------------------------------
	//	 Pre-simulation.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	correct_configuration_excitation.parallel_exec();
	correct_configuration_contraction.parallel_exec();
	correct_kernel_weights_for_interpolation.parallel_exec();
	/** Output initial states and observations */
	write_states.writeToFile(0);
	write_voltage.writeToFile(0);
	write_displacement.writeToFile(0);
	//----------------------------------------------------------------------
	//	 Physical parameters for main loop.
	//----------------------------------------------------------------------
	int screen_output_interval = 10;
	int ite = 0;
	int reaction_step = 2;
	Real End_Time = 100;
	Real Ouput_T = End_Time / 200.0;
	Real Observer_time = 0.01 * Ouput_T;
	Real dt = 0.0;	 /**< Default acoustic time step sizes for physiology. */
	Real dt_s = 0.0; /**< Default acoustic time step sizes for mechanics. */
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	std::cout << "Main Loop Starts Here : "
			  << "\n";
	/** Main loop starts here. */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		while (integration_time < Ouput_T)
		{
			Real relaxation_time = 0.0;
			while (relaxation_time < Observer_time)
			{
				if (ite % screen_output_interval == 0)
				{
					std::cout << std::fixed << std::setprecision(9) << "N=" << ite << "	Time = "
							  << GlobalStaticVariables::physical_time_
							  << "	dt = " << dt
							  << "	dt_s = " << dt_s << "\n";
				}
				/** Apply stimulus excitation. */
				if (0 <= GlobalStaticVariables::physical_time_ && GlobalStaticVariables::physical_time_ <= 0.5)
				{
					apply_stimulus_s1.parallel_exec(dt);
				}
				/** Single spiral wave. */
				// if( 60 <= GlobalStaticVariables::physical_time_
				// 	&&  GlobalStaticVariables::physical_time_ <= 65)
				// {
				// 	apply_stimulus_s2.parallel_exec(dt);
				// }
				/**Strong splitting method. */
				//forward reaction
				int ite_forward = 0;
				while (ite_forward < reaction_step)
				{
					reaction_relaxation_forward.parallel_exec(0.5 * dt / Real(reaction_step));
					ite_forward++;
				}
				/** 2nd Runge-Kutta scheme for diffusion. */
				diffusion_relaxation.parallel_exec(dt);

				//backward reaction
				int ite_backward = 0;
				while (ite_backward < reaction_step)
				{
					reaction_relaxation_backward.parallel_exec(0.5 * dt / Real(reaction_step));
					ite_backward++;
				}

				active_stress_interpolation.parallel_exec();

				Real dt_s_sum = 0.0;
				while (dt_s_sum < dt)
				{
					dt_s = get_mechanics_time_step.parallel_exec();
					if (dt - dt_s_sum < dt_s)
						dt_s = dt - dt_s_sum;
					stress_relaxation_first_half.parallel_exec(dt_s);
					constrain_holder.parallel_exec(dt_s);
					stress_relaxation_second_half.parallel_exec(dt_s);
					dt_s_sum += dt_s;
				}

				ite++;
				dt = get_physiology_time_step.parallel_exec();

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			write_voltage.writeToFile(ite);
			write_displacement.writeToFile(ite);
		}
		tick_count t2 = tick_count::now();
		interpolation_particle_position.parallel_exec();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
