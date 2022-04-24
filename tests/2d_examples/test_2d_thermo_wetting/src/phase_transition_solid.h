#include "SPHconfig.h"
#include "sphinxsys.h"

template<class lowTType, class highTType>
class phaseTransition : public lowTType, public highTType
{
public:
	Real transitionT_ = 0.5;
	Real epsilon_ = 1e-3;
	Real latentHeat = 3;
	phaseTransitionMaterial(libconfig::Setting* phaseTransitionConfig) : 
	DiffusionReaction<lowTParticlesType, lowTType>((*phaseTransitionConfig)["lowTphase"]),
	DiffusionReaction<highTParticlesType, highTType>((*phaseTransitionConfig)["highTphase"])
	{
		(*phaseTransitionConfig)["transitionConfig"].lookupValue("transitionT",transitionT_);
		(*phaseTransitionConfig)["transitionConfig"].lookupValue("epsilon_",epsilon_);
		(*phaseTransitionConfig)["transitionConfig"].lookupValue("latentHeat",latentHeat);
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
							DiffusionReaction<SolidParticles, Solid>,
							DiffusionReaction<CompressibleFluidParticles, vdWFluid>
						> = vdWphaseTransitionMaterial

phaseTransitionParticles<	
							DiffusionReactionParticles<SolidParticles, Solid>,
							DiffusionReactionParticles<CompressibleFluidParticles, vdWFluid>
						> (water_block, phaseTransitionConfig)

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