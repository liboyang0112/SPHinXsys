#include "SPHconfig.h"
#include "sphinxsys.h"

template<class lowTType, class highTType>
class phaseTransition : public lowTType, public highTType
{
public:
	Real transitionT_ = 0.5;
	Real epsilon_ = 1e-3;
	Real latent_heat_ = 3;
	phaseTransitionMaterial(libconfig::Setting* phaseTransitionConfig) : 
	DiffusionReaction<lowTParticlesType, lowTType>((*phaseTransitionConfig)["lowTphase"]),
	DiffusionReaction<highTParticlesType, highTType>((*phaseTransitionConfig)["highTphase"])
	{
		(*phaseTransitionConfig)["transitionConfig"].lookupValue("transitionT",transitionT_);
		(*phaseTransitionConfig)["transitionConfig"].lookupValue("epsilon_",epsilon_);
		(*phaseTransitionConfig)["transitionConfig"].lookupValue("latent_heat_",latent_heat_);
	};
}


template<class lowTParticlesType, class highTParticlesType, class phase1BodyType, class phase2BodyType>
class phaseTransitionParticles
{
public:
	unique_ptr<lowTParticlesType> lowTParticles_;
	unique_ptr<highTParticlesType> highTParticles_;
	shared_ptr<phase1BodyType> sph_body_phase1_;
	shared_ptr<phase2BodyType> sph_body_phase2_;
	Real T0;
	shared_ptr<phaseTransitionConfig>
	phaseTransitionParticles(SPHSystem &sph_system, SPHBody &sph_body_phase1, libconfig::Setting* phaseTransitionConfig) : 
	{
		sph_body_phase1_ = &sph_body_phase1;
		sph_body_phase2_ = make_shared<phase2BodyType>(sph_system,"TransitionPhaseBody", sph_body_phase1.sph_adaptation_);
		if(T>transitionT_) {
			lowTParticles_ = make_unique<lowTParticlesType>(*sph_body_phase2_,(*phaseTransitionConfig)["lowTphase"])
			highTParticles_ = make_unique<lowTParticlesType>(*sph_body_phase1_,(*phaseTransitionConfig)["highTphase"])
		}else{
			lowTParticles_ = make_unique<lowTParticlesType>(*sph_body_phase1_,(*phaseTransitionConfig)["lowTphase"])
			highTParticles_ = make_unique<lowTParticlesType>(*sph_body_phase2_,(*phaseTransitionConfig)["highTphase"])
		}
		(*phaseTransitionConfig)["transitionConfig"].lookupValue("T0",transitionT_);

	};
	getFluidBody
}

phaseTransitionParticles

using phaseTransitionMaterial<
							DiffusionReaction<ElasticSolidParticles, ElasticSolid>,
							DiffusionReaction<CompressibleFluidParticles, vdWFluid>
						> = vdWphaseTransitionMaterial

class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &sph_system, const string &body_name, shared_ptr<SPHAdaptation> adp)
		: FluidBody(sph_system, body_name, adp)
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};


phaseTransitionDynamics<	
							DiffusionReactionParticles<SolidParticles, Solid>,
							DiffusionReactionParticles<CompressibleFluidParticles, vdWFluid>
						> (p1, p2, phaseTransitionConfig){
							phaseTransitionDynamicsLowT d1;
							phaseTransitionDynamicsHighT d1;
						}

phaseTransitionDynamics::exec(){
	d1.exec();
	d2.exec();
}

static enum PhaseIndex{lowT=-1,T0=0,highT=1};
class BasePhaseTransitionDynamics{
protected:
	processLatentHeat(int index_i, Real dT, PhaseIndex pi);
}
class phaseTransitionDynamicsLowT: public BasePhaseTransitionDynamics, public InteractionDynamics<void>{
	interaction(int index_i);
}
class phaseTransitionDynamicsHighT: public BasePhaseTransitionDynamics, public InteractionDynamics<void>{
	interaction(int index_i);
}

BasePhaseTransitionDynamics::processLatentHeat(int index_i, PhaseIndex pi){
	Real dT = (transfer_T_-temperature_[index_i])*pi;
	if(dT<0){
		if(latent_heat_storage_[index_i]==0){
			return;
		}else if(latent_heat_storage_[index_i]>-dT){
			latent_heat_storage_[index_i]+=dT;
			temperature_[index_i] = transfer_T_;
		}else {
			temperature_[index_i]-=pi*latent_heat_storage_[index_i];
			latent_heat_storage_[index_i] = 0;
		}
	}else{
		if((latent_heat_storage_[index_i]+=dT) >= latent_heat_){
			temperature_p2_[particleTransfer(index_i)] += pi*(latent_heat_storage_[index_i]-latent_heat_);
		}
	}
}
size_t BasePhaseTransitionDynamics::particleTransfer(int index_i){
	size_t new_index =　copyParticleToP2();
	destroyParticle();
	return new_index;
}
phaseTransitionDynamicsLowT::interaction(int index_i){
		processLatentHeat(index_i,lowT);
}
phaseTransitionDynamicshighT::interaction(int index_i){
		processLatentHeat(index_i,highT);
}