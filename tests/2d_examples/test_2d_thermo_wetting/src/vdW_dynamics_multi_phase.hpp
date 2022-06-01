/**
 * @file 	vdW_dynamics_multi_phase.hpp
 * @author	Boyang Li
 */

#pragma once

#include "vdW_dynamics_multi_phase.h"
//=================================================================================================//

namespace SPH
{
//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		template<class PressureRelaxationInnerType>
        vdWPressureRelaxationMultiPhase<PressureRelaxationInnerType>::
            vdWPressureRelaxationMultiPhase(BaseBodyRelationInner &inner_relation,
				BaseBodyRelationContact &contact_relation) :
				BasePressureRelaxationMultiPhase<PressureRelaxationInnerType>(inner_relation,
					contact_relation){}
		//=================================================================================================//
		template<class PressureRelaxationInnerType>
        vdWPressureRelaxationMultiPhase<PressureRelaxationInnerType>::
            vdWPressureRelaxationMultiPhase(ComplexBodyRelation &complex_relation) :
				vdWPressureRelaxationMultiPhase(complex_relation.inner_relation_, complex_relation.contact_relation_) {}
		//=================================================================================================//
		template<class PressureRelaxationInnerType>
        void vdWPressureRelaxationMultiPhase<PressureRelaxationInnerType>::Interaction(size_t index_i, Real dt)
		{
			PressureRelaxationInnerType::Interaction(index_i, dt);
			Real pres = this->rho_n_[index_i] * this->temperature_[index_i] * k_B;
			FluidState state_i(this->rho_n_[index_i], this->vel_n_[index_i], pres);
			Vecd acceleration(0.0);
			Vecd vel_derivative(0);
			Vecd acceleration_i(0);
			Real temp_change_rate = 0;
			for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
			{
				StdLargeVec<Real>& Vol_k = *(this->contact_Vol_[k]);
				StdLargeVec<Real>& rho_k = *(this->contact_rho_n_[k]);
				StdLargeVec<Real>& p_k = *(this->contact_p_[k]);
				StdLargeVec<Vecd>& vel_k = *(this->contact_vel_n_[k]);
				StdLargeVec<Real>& temp_k = *getParticleTemperature(this->contact_particles_[k]);
				decltype(PressureRelaxationInnerType::riemann_solver_)& riemann_solver_k = this->riemann_solvers_[k];
				Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
				Real mu_j = this->contact_material_[k]->ReferenceViscosity();
				for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
				{
					size_t index_j = contact_neighborhood.j_[n];
					Vecd& e_ij = contact_neighborhood.e_ij_[n];
					Real dW_ij = contact_neighborhood.dW_ij_[n];
					Real r_ij = contact_neighborhood.r_ij_[n];
					Vecd vij = this->vel_n_[index_i] - vel_k[index_j];
					Real p_j = rho_k[index_j] * k_B * temp_k[index_j];
					FluidState state_j(rho_k[index_j], vel_k[index_j], p_j);
					Real p_star = riemann_solver_k.getPStar(state_i, state_j, e_ij);
					acceleration_i = -2.0 * p_star * e_ij * Vol_k[index_j] * dW_ij / state_i.rho_;
					vel_derivative = 2.0 * (state_i.vel_ - vel_k[index_j]) /
									 (r_ij + 0.01 * this->smoothing_length_);
					Real mu_ij = 2.0 * this->material_->ReferenceViscosity() * mu_j / (this->material_->ReferenceViscosity() + mu_j);
					acceleration_i += 2.0 * mu_ij * vel_derivative *
									dW_ij * Vol_k[index_j] / state_i.rho_;
					acceleration += acceleration_i;	
					temp_change_rate -= 0.5*dot(acceleration_i, vij);
				}
			}
			this->dvel_dt_[index_i] += acceleration;
			this->temperature_[index_i] += temp_change_rate * 2 / dof / k_B * dt;
		}
       //=================================================================================================//
	}
//=================================================================================================//
}
//=================================================================================================//