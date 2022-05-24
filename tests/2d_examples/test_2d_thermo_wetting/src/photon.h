#ifndef PHOTON
#define PHOTON
#include "base_material.h"
#include "photon_particles.h"
#include "SPHconfig.h"

namespace SPH
{
	class Photon : public BaseMaterial
	{
	protected:
		Real c0, lambda_; // reference sound speed, wavelength. 
		Real threshold_kill_; // When the intensity is lower than threshold, the particle will be killed.
		Real grad_threshold_reflect_; // Reflect only when the density gradiant is higher than threshold.
		Real rho_threshold_reflect_; // Reflect only when the density is higher than threshold.
		PhotonParticles *photon_particles_;

	public:
		explicit Photon(libconfig::Setting* photonConfig)
			: BaseMaterial(0), photon_particles_(nullptr)
		{
			material_type_ = "Photon";
			photonConfig->lookupValue("c0",c0);
			photonConfig->lookupValue("threshold_kill_",threshold_kill_);
			photonConfig->lookupValue("grad_threshold_reflect_",grad_threshold_reflect_);
			photonConfig->lookupValue("rho_threshold_reflect_",rho_threshold_reflect_);
		};
		virtual ~Photon(){};

		void assignPhotonParticles(PhotonParticles *photon_particles)
		{
			photon_particles_ = photon_particles;
		};

		virtual Real getLightSpeed(Real n=1) {return n*c0; };
		virtual Real getWaveLength(Real n=1) {return lambda_/n;};
		virtual Real getReflectivity(Real par1=0, Real par2=0) { return 0.9; };
		virtual Real getAttenuation(Real par1=0, Real par2=0) { return 0.3; };
		virtual Real getKillThreshold(Real par1=0, Real par2=0) { return 1e-5; };
		virtual Real getGradThresholdReflect(Real par1=0, Real par2=0) { return 0.01; };
		virtual Photon *ThisObjectPtr() override { return this; };
	};
}
#endif