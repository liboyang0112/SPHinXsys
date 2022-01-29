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
* @file 	contact_dynamics.h
* @brief 	Here, we define the algorithm classes for solid contact dynamics. 
* @details 	We consider here a weakly compressible solids.   
* @author	Chi Zhang and Xiangyu Hu
*/

#ifndef CONTACT_DYNAMICS_H
#define CONTACT_DYNAMICS_H

#include "solid_dynamics.h"

namespace SPH
{
	namespace solid_dynamics
	{
		/**
		* @class SelfContactDensitySummation
		* @brief Computing the summation density due to solid self-contact model.
		*/
		class SelfContactDensitySummation : public PartInteractionDynamicsByParticle, public SolidDataInner
		{
		public:
			explicit SelfContactDensitySummation(SolidBodyRelationSelfContact &self_contact_relation);
			virtual ~SelfContactDensitySummation(){};

		protected:
			StdLargeVec<Real> &mass_, &contact_density_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ContactDensitySummation
		* @brief Computing the summation density due to solid-solid contact model.
		*/
		class ContactDensitySummation : public PartInteractionDynamicsByParticle, public ContactDynamicsData
		{
		public:
			explicit ContactDensitySummation(SolidBodyRelationContact &solid_body_contact_relation);
			virtual ~ContactDensitySummation(){};

		protected:
			StdLargeVec<Real> &mass_, &contact_density_;
			StdVec<StdLargeVec<Real> *> contact_mass_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class SelfContactForce
		* @brief Computing the self-contact force.
		*/
		class SelfContactForce : public PartInteractionDynamicsByParticle, public SolidDataInner
		{
		public:
			explicit SelfContactForce(SolidBodyRelationSelfContact &self_contact_relation);
			virtual ~SelfContactForce(){};

		protected:
			StdLargeVec<Real> &mass_, &contact_density_, &Vol_;
			StdLargeVec<Vecd> &dvel_dt_prior_, &contact_force_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ContactForce
		* @brief Computing the contact force.
		*/
		class ContactForce : public PartInteractionDynamicsByParticle, public ContactDynamicsData
		{
		public:
			explicit ContactForce(SolidBodyRelationContact &solid_body_contact_relation);
			virtual ~ContactForce(){};

		protected:
			StdLargeVec<Real> &contact_density_, &Vol_, &mass_;
			StdLargeVec<Vecd> &dvel_dt_prior_, &contact_force_;
			StdVec<StdLargeVec<Real> *> contact_contact_density_, contact_Vol_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class DynamicContactForce
		* @brief Computing the contact force for problems dominated by the contact dynamic process itself.
		* For example, the high speed impact problems in which the detailed contact behavior is crucial for 
		* physical sound solutions. Therefore, for simple low speed problem in which contact force is 
		* used merely prevent penetration. We can still use the simple formulation in the class ContactForce. 
		* The idea is to introduce conact force based on Riemann problem like formulation, 
		* in which the artificial dissipation is the main interaction force to prevent
		* penetration. Furthermore, a penalty type force is used as supplementary to prevent penetration 
		* when the contact velocity is small.
		*/
		class DynamicContactForce : public PartInteractionDynamicsByParticle, public ContactDynamicsData
		{
		public:
			explicit DynamicContactForce(SolidBodyRelationContact &solid_body_contact_relation, Real penalty_strength = 1.0);
			virtual ~DynamicContactForce(){};

		protected:
			StdLargeVec<Real> &Vol_, &mass_;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_prior_, &contact_force_;
			StdVec<StdLargeVec<Real> *> contact_Vol_;
			StdVec<StdLargeVec<Vecd> *> contact_vel_n_;
			Real penalty_strength_;
			StdVec<Real> contact_impedance_, contact_reference_pressure_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		/**
		* @class ContactForceWithWall
		* @brief Computing the contact force with a rigid wall.
		*  Note that the body surface of the wall should be
		*  updated before computing the contact force.
		*/
		class ContactForceWithWall : public PartInteractionDynamicsByParticle, public ContactDynamicsData
		{
		public:
			explicit ContactForceWithWall(SolidBodyRelationContact &solid_body_contact_relation, Real penalty_strength = 1.0);
			virtual ~ContactForceWithWall(){};

		protected:
			StdLargeVec<Real> &Vol_, &mass_;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_prior_, &contact_force_;
			StdVec<StdLargeVec<Real> *> contact_Vol_;
			StdVec<StdLargeVec<Vecd> *> contact_vel_n_, contact_n_;
			Real penalty_strength_;
			Real impedance_, reference_pressure_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};

		class DynamicSelfContactForce : public PartInteractionDynamicsByParticle, public SolidDataInner
		{
		public:
			explicit DynamicSelfContactForce(SolidBodyRelationSelfContact &self_contact_relation, Real penalty_strength = 1.0);
			virtual ~DynamicSelfContactForce(){};

		protected:
			StdLargeVec<Real> &Vol_, &mass_;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_prior_, &contact_force_;
			Real penalty_strength_;
			Real particle_spacing_j1_, particle_spacing_ratio2_;
			Real contact_impedance_, contact_reference_pressure_;

			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
	}
}
#endif //CONTACT_DYNAMICS_H