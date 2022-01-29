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
* @file 	solid_dynamics.h
* @brief 	Here, we define the algorithm classes for solid dynamics. 
* @details 	We consider here a weakly compressible solids.   
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
*/

#ifndef SOLID_DYNAMICS_H
#define SOLID_DYNAMICS_H

#include "all_particle_dynamics.h"
#include "general_dynamics.h"
#include "elastic_solid.h"
#include "base_kernel.h"

namespace SPH
{
	template <int DataTypeIndex, typename VariableType>
	class BodySummation;
	template <int DataTypeIndex, typename VariableType>
	class BodyMoment;

	namespace solid_dynamics
	{
		//----------------------------------------------------------------------
		//		for general solid dynamics
		//----------------------------------------------------------------------
		typedef DataDelegateSimple<SolidBody, SolidParticles, Solid> SolidDataSimple;
		typedef DataDelegateInner<SolidBody, SolidParticles, Solid> SolidDataInner;
		typedef DataDelegateContact<SolidBody, SolidParticles, Solid, SolidBody, SolidParticles, Solid> ContactDynamicsData;

		/**
		 * @class SolidDynamicsInitialCondition
		 * @brief  set initial condition for solid fluid body
		 * This is a abstract class to be override for case specific initial conditions.
		 */
		class SolidDynamicsInitialCondition : public ParticleDynamicsSimple, public SolidDataSimple
		{
		public:
			explicit SolidDynamicsInitialCondition(SolidBody &solid_body)
				: ParticleDynamicsSimple(solid_body), SolidDataSimple(solid_body){};
			virtual ~SolidDynamicsInitialCondition(){};
		};

		/**
		* @class CorrectConfiguration
		* @brief obtain the corrected initial configuration in strong form
		*/
		class CorrectConfiguration : public InteractionDynamics, public SolidDataInner
		{
		public:
			explicit CorrectConfiguration(BaseBodyRelationInner &inner_relation);
			virtual ~CorrectConfiguration(){};

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Matd> &B_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ConstrainSolidBodyRegion
		 * @brief Constrain a solid body part with prescribed motion.
		 * Note the average values for FSI are prescirbed also.
		 */
		class ConstrainSolidBodyRegion : public PartSimpleDynamicsByParticle, public SolidDataSimple
		{
		public:
			//TODO: use only body part as argment
			ConstrainSolidBodyRegion(SPHBody &sph_body, BodyPartByParticle &body_part);
			virtual ~ConstrainSolidBodyRegion(){};

		protected:
			StdLargeVec<Vecd> &pos_n_, &pos_0_;
			StdLargeVec<Vecd> &n_, &n_0_;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_, &vel_ave_, &dvel_dt_ave_;
			virtual Vecd getDisplacement(Vecd &pos_0, Vecd &pos_n) { return pos_n; };
			virtual Vecd getVelocity(Vecd &pos_0, Vecd &pos_n, Vecd &vel_n) { return Vecd(0); };
			virtual Vecd getAcceleration(Vecd &pos_0, Vecd &pos_n, Vecd &dvel_dt) { return Vecd(0); };
			virtual SimTK::Rotation getBodyRotation(Vecd &pos_0, Vecd &pos_n, Vecd &dvel_dt) { return SimTK::Rotation(); }
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ConstrainSolidBodySurfaceRegion
		 * @brief Constrain the surface particles of a solid body part with prescribed motion.
		 */
		class ConstrainSolidBodySurfaceRegion :
			public PartSimpleDynamicsByParticle, public SolidDataSimple
		{
		public:
			ConstrainSolidBodySurfaceRegion(SPHBody& body, BodyPartByParticle& body_part);
			virtual ~ConstrainSolidBodySurfaceRegion() {};

			StdLargeVec<bool>& GetApplyConstrainToParticle(){ return apply_constrain_to_particle_; }
		protected:
			StdLargeVec<Vecd>& pos_n_, & pos_0_;
			StdLargeVec<Vecd>& vel_n_, & dvel_dt_;
			StdLargeVec<bool> apply_constrain_to_particle_;

			virtual Vecd getDisplacement(Vecd& pos_0, Vecd& pos_n) { return pos_n; };
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class PositionSolidBody
		 * @brief Moves the body into a defined position in a given time interval - position driven boundary condition
		 * Note the average values for FSI are prescirbed also.
		 */
		class PositionSolidBody : public PartSimpleDynamicsByParticle, public SolidDataSimple
		{
		public:
			PositionSolidBody(SPHBody &sph_body, BodyPartByParticle &body_part, Real start_time, Real end_time, Vecd pos_end_center);
			virtual ~PositionSolidBody(){};
			StdLargeVec<Vecd> &GetParticlePos0() { return pos_0_; };
			StdLargeVec<Vecd> &GetParticlePosN() { return pos_n_; };

		protected:
			StdLargeVec<Vecd> &pos_n_, &pos_0_;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_, &vel_ave_, &dvel_dt_ave_;
			Real start_time_, end_time_;
			Vecd pos_0_center_, pos_end_center_, translation_;
			Vecd getDisplacement(size_t index_i, Real dt);
			virtual Vecd getVelocity() { return Vecd(0); };
			virtual Vecd getAcceleration() { return Vecd(0); };
			virtual SimTK::Rotation getBodyRotation() { return SimTK::Rotation(); }
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class PositionScaleSolidBody
		 * @brief Scales the body in a given time interval - position driven boundary condition
		 * Note the average values for FSI are prescirbed also.
		 */
		class PositionScaleSolidBody : public PartSimpleDynamicsByParticle, public SolidDataSimple
		{
		public:
			PositionScaleSolidBody(SPHBody &sph_body, BodyPartByParticle &body_part, Real start_time, Real end_time, Real end_scale);
			virtual ~PositionScaleSolidBody(){};
			StdLargeVec<Vecd> &GetParticlePos0() { return pos_0_; };
			StdLargeVec<Vecd> &GetParticlePosN() { return pos_n_; };

		protected:
			StdLargeVec<Vecd> &pos_n_, &pos_0_;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_, &vel_ave_, &dvel_dt_ave_;
			Real start_time_, end_time_, end_scale_;
			Vecd pos_0_center_;
			Vecd getDisplacement(size_t index_i, Real dt);
			virtual Vecd getVelocity() { return Vecd(0); };
			virtual Vecd getAcceleration() { return Vecd(0); };
			virtual SimTK::Rotation getBodyRotation() { return SimTK::Rotation(); }
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class TranslateSolidBody
		 * @brief Translates the body in a given time interval -translation driven boundary condition; only moving the body; end position irrelevant;
		 * Note the average values for FSI are prescirbed also.
		 */
		class TranslateSolidBody : public PartSimpleDynamicsByParticle, public SolidDataSimple
		{
		public:
			TranslateSolidBody(SPHBody &sph_body, BodyPartByParticle &body_part, Real start_time, Real end_time, Vecd translation);
			virtual ~TranslateSolidBody(){};

		protected:
			StdLargeVec<Vecd> &pos_n_, &pos_0_;
			StdLargeVec<Vecd> pos_end_;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_;
			Real start_time_, end_time_;
			Vecd translation_;
			Vecd getDisplacement(size_t index_i, Real dt);
			virtual Vecd getVelocity() { return Vecd(0); };
			virtual Vecd getAcceleration() { return Vecd(0); };
			virtual SimTK::Rotation getBodyRotation() { return SimTK::Rotation(); }
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class TranslateSolidBodyPart
		 * @brief Translates the body in a given time interval -translation driven boundary condition; only moving the body; end position irrelevant;
		 * Only the particles in a given Bounding Box are translated. The Bounding Box is defined for the nondeformed shape.
		 * Note the average values for FSI are prescirbed also.
		 */
		class TranslateSolidBodyPart : public TranslateSolidBody
		{
		public:
			TranslateSolidBodyPart(SPHBody &sph_body, BodyPartByParticle &body_part, Real start_time, Real end_time, Vecd translation, BoundingBox bbox);
			virtual ~TranslateSolidBodyPart(){};

		protected:
			BoundingBox bbox_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ConstrainSolidBodyRegionVelocity
		 * @brief Constrain the velocity of a solid body part.
		 */
		class ConstrainSolidBodyRegionVelocity : public ConstrainSolidBodyRegion
		{
		public:
			ConstrainSolidBodyRegionVelocity(SPHBody &sph_body, BodyPartByParticle &body_part,
											 Vecd constrained_direction = Vecd(0))
				: solid_dynamics::ConstrainSolidBodyRegion(sph_body, body_part),
				  constrain_matrix_(Matd(1.0))
			{
				for (int k = 0; k != Dimensions; ++k)
					constrain_matrix_[k][k] = constrained_direction[k];
			};
			virtual ~ConstrainSolidBodyRegionVelocity(){};

		protected:
			Matd constrain_matrix_;
			virtual Vecd getVelocity(Vecd &pos_0, Vecd &pos_n, Vecd &vel_n)
			{
				return constrain_matrix_ * vel_n;
			};
		};

		/**
		 * @class SoftConstrainSolidBodyRegion
		 * @brief Soft the constrain of a solid body part
		 */
		class SoftConstrainSolidBodyRegion : public PartInteractionDynamicsByParticleWithUpdate,
											 public SolidDataInner
		{
		public:
			SoftConstrainSolidBodyRegion(BaseBodyRelationInner &inner_relation, BodyPartByParticle &body_part);
			virtual ~SoftConstrainSolidBodyRegion(){};

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_, &vel_ave_, &dvel_dt_ave_;
			StdLargeVec<Vecd> vel_temp_, dvel_dt_temp_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ClampConstrainSolidBodyRegion
		 * @brief Constrain a solid body part with prescribed motion and smoothing to mimic the clamping effect.
		 */
		class ClampConstrainSolidBodyRegion : public ParticleDynamics<void>
		{
		public:
			ConstrainSolidBodyRegion constrianing_;
			SoftConstrainSolidBodyRegion softing_;

			ClampConstrainSolidBodyRegion(BaseBodyRelationInner &inner_relation, BodyPartByParticle &body_part);
			virtual ~ClampConstrainSolidBodyRegion(){};

			virtual void exec(Real dt = 0.0) override;
			virtual void parallel_exec(Real dt = 0.0) override;
		};

		/**
		 * @class ConstrainSolidBodyMassCenter
		 * @brief Constrain the mass center of a solid body.
		 */
		class ConstrainSolidBodyMassCenter : public ParticleDynamicsSimple, public SolidDataSimple
		{
		public:
			explicit ConstrainSolidBodyMassCenter(SPHBody &sph_body, Vecd constrain_direction = Vecd(1.0));
			virtual ~ConstrainSolidBodyMassCenter(){};

		protected:
			virtual void setupDynamics(Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;

		private:
			Real total_mass_;
			Matd correction_matrix_;
			Vecd velocity_correction_;
			StdLargeVec<Vecd> &vel_n_;
			BodyMoment<indexVector, Vecd> compute_total_momentum_;
		};

		/**@class ImposeExternalForce
		 * @brief impose external force on a solid body part
		 * by add extra acceleration
		 */
		class ImposeExternalForce : public PartSimpleDynamicsByParticle, public SolidDataSimple
		{
		public:
			ImposeExternalForce(SolidBody &solid_body, SolidBodyPartForSimbody &body_part);
			virtual ~ImposeExternalForce(){};

		protected:
			StdLargeVec<Vecd> &pos_0_, &vel_n_, &vel_ave_;
			/**
			 * @brief acceleration will be specified by the application
			 */
			virtual Vecd getAcceleration(Vecd &pos) = 0;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};
		/**
		* @class SpringDamperConstraintParticleWise
		* @brief Exerts spring force and damping force in the form of acceleration to each particle.
		* The spring force is calculated based on the difference from the particle's initial position.
		* The damping force is calculated based on the particle's current velocity.
		* Only for 3D applications
		*/
		class SpringDamperConstraintParticleWise
			: public ParticleDynamicsSimple,
			  public SolidDataSimple
		{
		public:
			SpringDamperConstraintParticleWise(SolidBody &solid_body, Vecd stiffness, Real damping_ratio = 0.05);

		protected:
			StdLargeVec<Vecd> &pos_n_, &pos_0_, &vel_n_, &dvel_dt_prior_;
			Vecd stiffness_;
			Vecd damping_coeff_; // damping component parallel to the spring force component

			virtual Vecd getSpringForce(size_t index_i, Vecd &disp);
			virtual Vecd getDampingForce(size_t index_i);
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};
		/**
		* @class SpringNormalOnSurfaceParticles
		* @brief Exerts spring force force on the surface in normal direction in the form of acceleration to each particle.
		* The input stiffness should be defined in Pa/m. The stiffness is scaled by the surface area of the particle to get N/m
		* The force is applied to all the surface particles that can be seen (outer_surface = false)
		* or cannot be seen (outer_surface = true) from the source point.
		* Can be used for outer or inner surface of a shell structure ofr example.
		* The spring force is calculated based on the difference from the particle's initial position.
		* Only for 3D applications
		* Only for uniform surface particle size.
		*/
		class SpringNormalOnSurfaceParticles
			: public PartSimpleDynamicsByParticle,
			  public SolidDataSimple
		{
		public:
			SpringNormalOnSurfaceParticles(SolidBody &solid_body, BodyPartByParticle &body_part,
										   bool outer_surface, Vecd source_point, Real stiffness, Real damping_ratio = 0.05);

			StdLargeVec<bool> &GetApplySpringForceToParticle() { return apply_spring_force_to_particle_; }

		protected:
			StdLargeVec<Vecd> &pos_n_, &pos_0_, &n_, &n_0_, &vel_n_, &dvel_dt_prior_;
			StdLargeVec<Real> &mass_;
			Real stiffness_;
			Real damping_coeff_; // damping component parallel to the spring force component
			StdLargeVec<bool> apply_spring_force_to_particle_;

			virtual Vecd getSpringForce(size_t index_i, Vecd disp);
			virtual Vecd getDampingForce(size_t index_i);
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};
		/**
		* @class SpringOnSurfaceParticles
		* @brief Exerts spring force force on the surface in the form of acceleration to each particle.
		* The input stiffness should be defined in Pa/m. The stiffness is scaled by the surface area of the particle to get N/m
		* The force is applied to all the surface particles.
		* The spring force is calculated based on the difference from the particle's initial position.
		* Only for 3D applications
		* BodyPartByParticle define the ody part that the spring is applied to.
		* Only for uniform surface particle size.
		*/
		class SpringOnSurfaceParticles
			: public ParticleDynamicsSimple, public SolidDataSimple
		{
		public:
			SpringOnSurfaceParticles(SolidBody& body, Real stiffness, Real damping_ratio = 0.05);

			StdLargeVec<bool>& GetApplySpringForceToParticle(){ return apply_spring_force_to_particle_; }
		protected:
			StdLargeVec<Vecd>& pos_n_,& pos_0_,& vel_n_,& dvel_dt_prior_;
			StdLargeVec<Real>& mass_;
			Real stiffness_;
			Real damping_coeff_; // damping component parallel to the spring force component
			StdLargeVec<bool> apply_spring_force_to_particle_;

			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};
		/**
		* @class AccelerationForBodyPartInBoundingBox
		* @brief Adds acceleration to the part of the body that's inside a bounding box
		*/
		class AccelerationForBodyPartInBoundingBox
			: public ParticleDynamicsSimple,
			  public SolidDataSimple
		{
		public:
			AccelerationForBodyPartInBoundingBox(SolidBody &solid_body, BoundingBox &bounding_box, Vecd acceleration);
			virtual ~AccelerationForBodyPartInBoundingBox(){};

		protected:
			StdLargeVec<Vecd> &pos_n_, &dvel_dt_prior_;
			BoundingBox bounding_box_;
			Vecd acceleration_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ForceInBodyRegion
		* @brief ForceInBodyRegion, distributes the force vector as acceleration among the particles in a given body part
		*/
		class ForceInBodyRegion : public PartSimpleDynamicsByParticle, public SolidDataSimple
		{
		public:
			ForceInBodyRegion(SPHBody &sph_body, BodyPartByParticle &body_part, Vecd force, Real end_time);

		protected:
			StdLargeVec<Vecd> &pos_0_, &dvel_dt_prior_;
			Vecd acceleration_;
			Real end_time_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class SurfacePressureFromSource
		* @brief SurfacePressureFromSource, applies pressure on the surface particles coming from a source point
		*/
		class SurfacePressureFromSource : public PartSimpleDynamicsByParticle, public SolidDataSimple
		{
		public:
			SurfacePressureFromSource(SPHBody &sph_body, BodyPartByParticle &body_part,
									  Vecd source_point, StdVec<array<Real, 2>> pressure_over_time);

			StdLargeVec<bool> &GetApplyPressureToParticle() { return apply_pressure_to_particle_; }

		protected:
			StdLargeVec<Vecd> &pos_0_, &n_, &dvel_dt_prior_;
			StdLargeVec<Real> &mass_;
			StdVec<array<Real, 2>> pressure_over_time_;
			StdLargeVec<bool> apply_pressure_to_particle_;
			Real getPressure();
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		//----------------------------------------------------------------------
		//		for elastic solid dynamics
		//----------------------------------------------------------------------
		typedef DataDelegateSimple<SolidBody, ElasticSolidParticles, ElasticSolid> ElasticSolidDataSimple;
		typedef DataDelegateInner<SolidBody, ElasticSolidParticles, ElasticSolid> ElasticSolidDataInner;

		/**
		 * @class ElasticDynamicsInitialCondition
		 * @brief  set initial condition for a solid body with different material
		 * This is a abstract class to be override for case specific initial conditions.
		 */
		class ElasticDynamicsInitialCondition : public ParticleDynamicsSimple, public ElasticSolidDataSimple
		{
		public:
			explicit ElasticDynamicsInitialCondition(SolidBody &solid_body);
			virtual ~ElasticDynamicsInitialCondition(){};

		protected:
			StdLargeVec<Vecd> &pos_n_, &vel_n_;
		};

		/**
		* @class UpdateElasticNormalDirection
		* @brief update particle normal directions for elastic solid
		*/
		class UpdateElasticNormalDirection : public ParticleDynamicsSimple, public ElasticSolidDataSimple
		{
		public:
			explicit UpdateElasticNormalDirection(SolidBody &solid_body);
			virtual ~UpdateElasticNormalDirection(){};

		protected:
			StdLargeVec<Vecd> &n_, &n_0_;
			StdLargeVec<Matd> &F_;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class AcousticTimeStepSize
		* @brief Computing the acoustic time step size
		* computing time step size
		*/
		class AcousticTimeStepSize : public ParticleDynamicsReduce<Real, ReduceMin>,
									 public ElasticSolidDataSimple
		{
		public:
			explicit AcousticTimeStepSize(SolidBody &solid_body, Real CFL = 0.6);
			virtual ~AcousticTimeStepSize(){};

		protected:
			Real CFL_;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_;
			Real smoothing_length_;
			Real ReduceFunction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class DeformationGradientTensorBySummation
		* @brief computing deformation gradient tensor by summation
		*/
		class DeformationGradientTensorBySummation : public InteractionDynamics, public ElasticSolidDataInner
		{
		public:
			explicit DeformationGradientTensorBySummation(BaseBodyRelationInner &inner_relation);
			virtual ~DeformationGradientTensorBySummation(){};

		protected:
			StdLargeVec<Real> &Vol_;
			StdLargeVec<Vecd> &pos_n_;
			StdLargeVec<Matd> &B_, &F_;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class BaseElasticRelaxation
		* @brief base class for elastic relaxation
		*/
		class BaseElasticRelaxation
			: public ParticleDynamics1Level,
			  public ElasticSolidDataInner
		{
		public:
			explicit BaseElasticRelaxation(BaseBodyRelationInner &inner_relation);
			virtual ~BaseElasticRelaxation(){};

		protected:
			StdLargeVec<Real> &Vol_, &rho_n_, &mass_;
			StdLargeVec<Vecd> &pos_n_, &vel_n_, &dvel_dt_;
			StdLargeVec<Matd> &B_, &F_, &dF_dt_;
		};

		/**
		* @class StressRelaxationFirstHalf
		* @brief computing stress relaxation process by verlet time stepping
		* This is the first step
		*/
		class StressRelaxationFirstHalf : public BaseElasticRelaxation
		{
		public:
			explicit StressRelaxationFirstHalf(BaseBodyRelationInner &inner_relation);
			virtual ~StressRelaxationFirstHalf(){};

		protected:
			Real rho0_, inv_rho0_;
			StdLargeVec<Vecd> &dvel_dt_prior_, &force_from_fluid_;
			StdLargeVec<Matd> &stress_PK1_;
			Real numerical_dissipation_factor_;
			Real smoothing_length_;
			Real inv_W0_ = 1.0 / body_->sph_adaptation_->getKernel()->W0(Vecd(0));

			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class KirchhoffParticleStressRelaxationFirstHalf
		*/
		class KirchhoffParticleStressRelaxationFirstHalf : public StressRelaxationFirstHalf
		{
		public:
			explicit KirchhoffParticleStressRelaxationFirstHalf(BaseBodyRelationInner &inner_relation);
			virtual ~KirchhoffParticleStressRelaxationFirstHalf(){};

		protected:
			const Real one_over_dimensions_ = 1.0 / (Real)Dimensions;

			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class KirchhoffStressRelaxationFirstHalf
		* @brief Decompose the stress into particle stress includes isotropic stress 
		* and the stress due to non-homogeneous material properties.
		* The preliminary shear stress is introduced by particle pair to avoid 
		* spurious stress and deformation. 
		* Note that, for the shear stress term, 
		* due to the mismatch of the divergence contribution between 
		* the pair-wise second-order derivative Laplacian formulation
		* and particle-wise first-order gradient formulation, 
		* a correction factor slight large than one is introduced.       
		*/
		class KirchhoffStressRelaxationFirstHalf : public StressRelaxationFirstHalf
		{
		public:
			explicit KirchhoffStressRelaxationFirstHalf(BaseBodyRelationInner &inner_relation);
			virtual ~KirchhoffStressRelaxationFirstHalf(){};

		protected:
			StdLargeVec<Real> J_to_minus_2_over_dimension_;
			StdLargeVec<Matd> stress_on_particle_, inverse_F_T_;
			const Real one_over_dimensions_ = 1.0 / (Real)Dimensions;
			const Real correction_factor_ = 1.05;

			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class StressRelaxationSecondHalf
		* @brief computing stress relaxation process by verlet time stepping
		* This is the second step
		*/
		class StressRelaxationSecondHalf : public BaseElasticRelaxation
		{
		public:
			explicit StressRelaxationSecondHalf(BaseBodyRelationInner &inner_relation)
				: BaseElasticRelaxation(inner_relation){};
			virtual ~StressRelaxationSecondHalf(){};

		protected:
			virtual void Initialization(size_t index_i, Real dt = 0.0) override;
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
			virtual void Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class ConstrainSolidBodyPartBySimBody
		 * @brief Constrain a solid body part from the motion
		 * computed from Simbody.
		 */
		class ConstrainSolidBodyPartBySimBody : public ConstrainSolidBodyRegion
		{
		public:
			ConstrainSolidBodyPartBySimBody(SolidBody &solid_body,
											SolidBodyPartForSimbody &body_part,
											SimTK::MultibodySystem &MBsystem,
											SimTK::MobilizedBody &mobod,
											SimTK::Force::DiscreteForces &force_on_bodies,
											SimTK::RungeKuttaMersonIntegrator &integ);
			virtual ~ConstrainSolidBodyPartBySimBody(){};

		protected:
			SimTK::MultibodySystem &MBsystem_;
			SimTK::MobilizedBody &mobod_;
			SimTK::Force::DiscreteForces &force_on_bodies_;
			SimTK::RungeKuttaMersonIntegrator &integ_;
			const SimTK::State *simbody_state_;
			Vec3d initial_mobod_origin_location_;

			virtual void setupDynamics(Real dt = 0.0) override;
			void virtual Update(size_t index_i, Real dt = 0.0) override;
		};

		/**
		 * @class TotalForceOnSolidBodyPartForSimBody
		 * @brief Compute the force acting on the solid body part
		 * for applying to simbody forces latter
		 */
		class TotalForceOnSolidBodyPartForSimBody
			: public PartDynamicsByParticleReduce<SimTK::SpatialVec, ReduceSum<SimTK::SpatialVec>>,
			  public SolidDataSimple
		{
		public:
			TotalForceOnSolidBodyPartForSimBody(SolidBody &solid_body,
												SolidBodyPartForSimbody &body_part,
												SimTK::MultibodySystem &MBsystem,
												SimTK::MobilizedBody &mobod,
												SimTK::Force::DiscreteForces &force_on_bodies,
												SimTK::RungeKuttaMersonIntegrator &integ);
			virtual ~TotalForceOnSolidBodyPartForSimBody(){};

		protected:
			StdLargeVec<Vecd> &force_from_fluid_, &contact_force_, &pos_n_;
			SimTK::MultibodySystem &MBsystem_;
			SimTK::MobilizedBody &mobod_;
			SimTK::Force::DiscreteForces &force_on_bodies_;
			SimTK::RungeKuttaMersonIntegrator &integ_;
			const SimTK::State *simbody_state_;
			Vec3d current_mobod_origin_location_;

			virtual void SetupReduce() override;
			virtual SimTK::SpatialVec ReduceFunction(size_t index_i, Real dt = 0.0) override;
		};
	}
}
#endif //SOLID_DYNAMICS_H