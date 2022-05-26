#ifndef __PHOTON_DYNAMICS
#define __PHOTON_DYNAMICS
#include <cassert>
#include "photon_particles.h"
#include "photon.h"
/**
 * @file 	photon_dynamics.h
 * @brief 	laser heating
 * @author 	Boyang Li
 */
/**
 * @brief 	SPHinXsys Library.
 */
/**
 * @brief Namespace cite here.
 */
using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */

class LaserField{
public:
	LaserField(){};
	virtual Vecd getDirection(Vecd pos){};
	virtual Real getIntensity(Vecd pos){};
};

class PlainWave : public LaserField{
public:
	std::string name = "PlainWave";
	Vecd direction_;
	Real intensity_;
	PlainWave(Vecd direction, Real intensity):LaserField(){
		intensity_ = intensity;
		direction_ = direction;
	}
	Vecd getDirection(Vecd pos){
		return direction_;
	}
	Real getIntensity(Vecd pos){
		return intensity_;
	}
};
class PhotonInitialization // Similar to TimeStepInitialization
	: public ParticleDynamicsSimple, public DataDelegateSimple<SPHBody,PhotonParticles,BaseMaterial>
{

public:
	PhotonInitialization(SPHBody &sph_body, LaserField &laser_field, SPHBody *contact_body = 0);
	virtual ~PhotonInitialization(){};
	virtual void exec(Real dt = 0.0) override{
		particles_->total_real_particles_ = total_photons;
		ParticleDynamicsSimple::exec(dt);
	};

protected:
	StdLargeVec<Vecd> &pos_n_;
	StdLargeVec<Vecd> pos_backup_;
	StdLargeVec<Vecd> &vel_n_;
	LaserField *laser_field_;
	SPHBody *contact_body_;
	size_t total_photons;
	virtual void setupDynamics(Real dt = 0.0) override {particles_->total_ghost_particles_ = 0;};
	virtual void Update(size_t index_i, Real dt = 0.0) override;
};

PhotonInitialization::PhotonInitialization(SPHBody &sph_body, LaserField &laser_field, SPHBody *contact_body)
	: ParticleDynamicsSimple(sph_body), DataDelegateSimple<SPHBody,PhotonParticles,BaseMaterial>(sph_body),
	  pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_),
	  laser_field_(&laser_field),
	  contact_body_(contact_body),
	  total_photons(particles_->total_real_particles_),
	  pos_backup_(particles_->pos_n_) {

	  }

void PhotonInitialization::Update(size_t index_i, Real dt)
{
	pos_n_[index_i] = pos_backup_[index_i];
	vel_n_[index_i] = laser_field_->getDirection(pos_n_[index_i]);
	particles_->intensity_n_[index_i] = laser_field_->getIntensity(pos_n_[index_i]);
	particles_->inside_body_n_[index_i] = contact_body_;
	particles_->rho_gradiant_prev_n_[index_i] = 0;
}

class PhotonPropagation // Similar to TimeStepInitialization
	: public InteractionDynamicsWithUpdate, 
	public DataDelegateContact<SPHBody,PhotonParticles,Photon,SPHBody,BaseParticles,BaseMaterial>
{

protected:
	StdLargeVec<Vecd> &pos_n_;
	StdLargeVec<Vecd> &vel_n_;
	StdLargeVec<Real> &intensity_n_;
	StdLargeVec<void*> &inside_body_n_;
	StdLargeVec<Real> &rho_gradiant_prev_n_;
	StdLargeVec<Real> &rho_gradiant_prev_max_n_;
	StdVec<StdLargeVec<Real>*> temperatures_;
	StdLargeVec<size_t> particle_kill_list_;
	StdLargeVec<size_t> particle_copy_list_;
	std::atomic_size_t n_particles_to_kill_;
	std::atomic_size_t n_particles_;
//	StdVec<Real> &refractive_index;
	SPHBody &sph_body_;
public:
	PhotonPropagation(BaseBodyRelationContact &body_relation_contact_):
	InteractionDynamicsWithUpdate(*body_relation_contact_.sph_body_),
	DataDelegateContact<SPHBody,PhotonParticles,Photon,SPHBody,BaseParticles,BaseMaterial>(body_relation_contact_),
	sph_body_(*body_relation_contact_.sph_body_), pos_n_(particles_->pos_n_), vel_n_(particles_->vel_n_),
	intensity_n_(particles_->intensity_n_),
	inside_body_n_(particles_->inside_body_n_),
	rho_gradiant_prev_n_(particles_->rho_gradiant_prev_n_),
	rho_gradiant_prev_max_n_(particles_->rho_gradiant_prev_max_n_),
	n_particles_(particles_->total_real_particles_){
		for(auto &body : contact_bodies_){
			temperatures_.push_back(getParticleTemperature(body->base_particles_));
		}
		particle_kill_list_.resize(particles_->total_real_particles_*2);
		particle_copy_list_.resize(particles_->total_real_particles_*2);
	};
	virtual ~PhotonPropagation(){};
	void updateList(){
		particles_->total_real_particles_ = n_particles_.load();
		if(n_particles_to_kill_.load()) killParticles(*particles_, n_particles_to_kill_.load(), particle_kill_list_, particle_copy_list_);
		particles_->total_real_particles_-=n_particles_to_kill_;
		n_particles_ = particles_->total_real_particles_;
		n_particles_to_kill_ = 0;
	}
	void exec(Real dt = 0.0){
		n_particles_ = particles_->total_real_particles_;
		InteractionDynamicsWithUpdate::exec(dt);
		updateList();
	}
	void parallel_exec(Real dt = 0.0){
		n_particles_ = particles_->total_real_particles_;
		InteractionDynamicsWithUpdate::parallel_exec(dt);
		updateList();
	}
protected:
	virtual void Interaction(size_t index_i, Real dt = 0.0) override;
	virtual void Update(size_t index_i, Real dt = 0.0) override;
};

void reflectParticle(size_t reflected_index_, size_t index_i, BaseParticles &particles_, Vecd surfacenorm){
	particles_.copyFromAnotherParticle(reflected_index_,index_i);
	particles_.vel_n_[reflected_index_] = particles_.vel_n_[index_i]-surfacenorm*2.*dot(surfacenorm,particles_.vel_n_[index_i])/pow(surfacenorm.norm(),2);
}

void PhotonPropagation::Interaction(size_t index_i, Real dt){
	for (size_t k = 0; k < contact_configuration_.size(); ++k)
	{
		Neighborhood &inside_neighbors = (*contact_configuration_[k])[index_i];
		Real densum = 0;
		Vecd rhograd(0);
		for (size_t n = 0; n != inside_neighbors.current_size_; ++n){
			size_t index_j = inside_neighbors.j_[n];
			densum+= inside_neighbors.W_ij_[n];
			rhograd += inside_neighbors.dW_ij_[n]*inside_neighbors.e_ij_[n];
		}
		Real gradnorm = dot(rhograd,vel_n_[index_i])/contact_bodies_[k]->base_particles_->base_material_->ReferenceDensity();
		if(contact_bodies_[k] == inside_body_n_[index_i]){
			Real heat_absorbed = intensity_n_[index_i]*material_->getLightSpeed()*dt/material_->getAttenuation(densum);
			for (size_t n = 0; n != inside_neighbors.current_size_; ++n){
				size_t index_j = inside_neighbors.j_[n];
				(*temperatures_[k])[index_j]+= inside_neighbors.W_ij_[n]*heat_absorbed/densum;
			}
			if((intensity_n_[index_i]-=heat_absorbed)<=material_->getKillThreshold()){
				particle_kill_list_[n_particles_to_kill_++] = index_i;
				return; //the particle is killed already, nothing to do next
			}
			//check if the photon is leaving the body
			if(-gradnorm > material_->getGradThresholdReflect()){
				if(gradnorm>rho_gradiant_prev_n_[index_i]){
					if(material_->getReflectivity() * intensity_n_[index_i] > material_->getKillThreshold()){
						Real reflected_index_ = n_particles_++;
						reflectParticle(reflected_index_, index_i, *particles_, rhograd);
						intensity_n_[reflected_index_] = material_->getReflectivity()*intensity_n_[index_i];
					}
					intensity_n_[index_i]*=1-material_->getReflectivity();
					inside_body_n_[index_i] = 0;
					gradnorm = 0;
					Interaction(index_i, dt); //check if entering another body
				}else{
					rho_gradiant_prev_n_[index_i] = gradnorm;
				}
			}
		}else {
			//check if the photon is entering the body
			if(gradnorm!=0){
				int halt = 1;
			}
			if(gradnorm > material_->getGradThresholdReflect()){
				if(gradnorm<rho_gradiant_prev_max_n_[index_i]-1e-5){
					if(inside_body_n_[index_i] ==0){
						Real reflectivity = particles_->vel_n_[index_i][1]<0.0?material_->getReflectivity():0;
						if(reflectivity * intensity_n_[index_i] > material_->getKillThreshold()){
							Real reflected_index_ = n_particles_++;
							reflectParticle(reflected_index_, index_i, *particles_, rhograd);
							rho_gradiant_prev_max_n_[reflected_index_] = 0;
							intensity_n_[reflected_index_] = reflectivity*intensity_n_[index_i];
						}
						intensity_n_[index_i]*=1-reflectivity;
						inside_body_n_[index_i] = contact_bodies_[k];
					}
				}else{
					rho_gradiant_prev_max_n_[index_i] = gradnorm;
				}
			}
		}
	}
}
void PhotonPropagation::Update(size_t index_i, Real dt){
	pos_n_[index_i]+=vel_n_[index_i]*dt*material_->getLightSpeed();
}
//=================================================================================================//

/*
	template <class BasePressureRelaxationType>
	class LaserHeating : public RelaxationWithWall<BasePressureRelaxationType>
	{
	public:
		// template for different combination of constructing body relations
		template <class BaseBodyRelationType>
		LaserHeating(BaseBodyRelationType &base_body_relation,
						   BaseBodyRelationContact &wall_contact_relation);
		virtual ~LaserHeating(){};
	protected:
		virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		virtual Vecd computeNonConservativeAcceleration(size_t index_i) override;
	};
	template <class BasePressureRelaxationType>
	template <class BaseBodyRelationType>
	LaserHeating<BasePressureRelaxationType>::
		LaserHeating(BaseBodyRelationType &base_body_relation,
						   BaseBodyRelationContact &wall_contact_relation)
		: RelaxationWithWall<BasePressureRelaxationType>(base_body_relation, wall_contact_relation) {}
	//=================================================================================================//
	template <class BasePressureRelaxationType>
	void LaserHeating<BasePressureRelaxationType>::Interaction(size_t index_i, Real dt)
	{
		BasePressureRelaxationType::Interaction(index_i, dt);
		Real pi_att = -((vdWFluid*)this->material_)->getAlpha()*this->rho_n_[index_i]*this->rho_n_[index_i];
		Real pi = this->p_[index_i]-pi_att;
		FluidState state_i(this->rho_n_[index_i], this->vel_n_[index_i], pi);
		FluidState state_i_att(this->rho_n_[index_i], this->vel_n_[index_i], pi_att);
		Vecd dvel_dt_prior_i = computeNonConservativeAcceleration(index_i);
		Vecd acceleration(0.0);
		Vecd acceleration_i;
		Real internalEnergyIncrease = 0;
		Vecd vel_derivative(0);
		for (size_t k = 0; k < FluidWallData::contact_configuration_.size(); ++k)
		{
			StdLargeVec<Real> &Vol_k = *(this->wall_Vol_[k]);
			StdLargeVec<Vecd> &vel_ave_k = *(this->wall_vel_ave_[k]);
			StdLargeVec<Vecd> &dvel_dt_ave_k = *(this->wall_dvel_dt_ave_[k]);
			StdLargeVec<Vecd> &n_k = *(this->wall_n_[k]);
			Neighborhood &wall_neighborhood = (*FluidWallData::contact_configuration_[k])[index_i];
			//Real rho_in_wall = 1.6;// + dp / this->material_->getSoundSpeed(state_i.p_ + 0.5 * dp, state_i.rho_); // sound speed = dp/drho
			//Real pj_att = ((vdWFluid*)this->material_)->getAlpha()*rho_in_wall*rho_in_wall;
			Real p_in_wall = 1.6;
			for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
			{
				size_t index_j = wall_neighborhood.j_[n];
				Vecd &e_ij = wall_neighborhood.e_ij_[n];
				Real dW_ij = wall_neighborhood.dW_ij_[n];
				//Real dW_ij_att = 0;
				Real dW_ij_att = wall_neighborhood.dW_ij_n_[0][n];
				Real r_ij = wall_neighborhood.r_ij_[n];
				Real face_wall_external_acceleration = dot((dvel_dt_prior_i - dvel_dt_ave_k[index_j]), -e_ij);
				Vecd vel_in_wall = 2.0 * vel_ave_k[index_j] - state_i.vel_;
				Vecd vij = state_i.vel_ - vel_in_wall;
				Real dp = state_i.rho_ * r_ij * SMAX(0.0, face_wall_external_acceleration);
				Real rho_in_wall = state_i.rho_;// + dp / this->material_->getSoundSpeed(state_i.p_ + 0.5 * dp, state_i.rho_); // sound speed = dp/drho
				//Real rho_in_wall = 1.6;// + dp / this->material_->getSoundSpeed(state_i.p_ + 0.5 * dp, state_i.rho_); // sound speed = dp/drho
				Real pj_att = -((vdWFluid*)this->material_)->getAlpha()*rho_in_wall*rho_in_wall;
				//Real p_in_wall = state_i.p_-pj_att;
				FluidState state_j(rho_in_wall, vel_in_wall, p_in_wall);
				//FluidState state_j_att(this->rho_n_[index_j], this->vel_n_[index_j], pj_att);
				//Real p_star = this->riemann_solver_.getPStar(state_i, state_j, n_k[index_j]);
				//Real p_star_att = this->riemann_solver_.getPStar(state_i_att, state_j_att, e_ij);
				vel_derivative = vij / (r_ij + 0.01 * this->smoothing_length_);
				//acceleration_i = -2.0 * (p_star * dW_ij) * e_ij * Vol_k[index_j] / state_i.rho_
				acceleration_i = -5.0 * dW_ij * e_ij
				+ 2.0 * this->material_->getViscosity(this->rho_n_[index_i], this->temperature_[index_i]) * vel_derivative * Vol_k[index_j] * dW_ij / state_i.rho_;
				internalEnergyIncrease -= 0.5*dot(acceleration_i, vij);
				acceleration += acceleration_i + 4 * e_ij * dW_ij_att;
			}
		}
		this->temperature_[index_i] += internalEnergyIncrease * 2 / dof / k_B * dt;
		this->dvel_dt_[index_i] += acceleration;
	}
	//=================================================================================================//
	template <class BasePressureRelaxationType>
	Vecd LaserHeating<BasePressureRelaxationType>::computeNonConservativeAcceleration(size_t index_i)
	{
		return this->dvel_dt_prior_[index_i];
	}
*/
#endif