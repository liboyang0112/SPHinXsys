#include "SPHconfig.h"
#include "sphinxsys.h"
#include <atomic>
using namespace SPH;
template <template <int DataTypeIndex, typename VariableType> typename OperationType>
class DualParticleDataOperation
{
public:
	OperationType<indexScalar, Real> scalar_operation;
	OperationType<indexVector, Vecd> vector_operation;
	OperationType<indexMatrix, Matd> matrix_operation;
	OperationType<indexInteger, int> integer_operation;
	OperationType<indexPointer, void*> pointer_operation;
	template<class ParticleType1, class ParticleType2 >
	DualParticleDataOperation(ParticleType1 &p1, ParticleType2 &p2) : 
	scalar_operation(p1,p2),vector_operation(p1,p2),
	matrix_operation(p1,p2),integer_operation(p1,p2),
	pointer_operation(p1,p2){}
	template <typename... ParticleArgs>
	void operator()(ParticleArgs... particle_args)
	{
		scalar_operation(particle_args...);
		vector_operation(particle_args...);
		matrix_operation(particle_args...);
		integer_operation(particle_args...);
		pointer_operation(particle_args...);
	}
};
template <int DataTypeIndex, typename VariableType>
class copyAParticleDataValueFromOtherParticles
{
public:
	ParticleData &particle_data_p1_;
	ParticleData &particle_data_p2_;
	StdVec<std::pair<size_t,size_t>> variableDataMapping;
	template<class ParticleType1, class ParticleType2 >
	copyAParticleDataValueFromOtherParticles(ParticleType1 &p1, ParticleType2 &p2):
	particle_data_p1_(p1.all_particle_data_),
	particle_data_p2_(p2.all_particle_data_){
		for(auto itr : p1.all_variable_maps_[DataTypeIndex]){
			auto ele = p2.all_variable_maps_[DataTypeIndex].find(itr.first);
			if(ele!=p2.all_variable_maps_[DataTypeIndex].end()){
				variableDataMapping.emplace_back(itr.second,ele->second);
			}
		}
	}
	void operator()(size_t this_index, size_t another_index) const;
};
template <int DataTypeIndex, typename VariableType>
void copyAParticleDataValueFromOtherParticles<DataTypeIndex, VariableType>::
operator()(size_t this_index, size_t another_index) const
{
    for (auto itr : variableDataMapping)
        (*std::get<DataTypeIndex>(particle_data_p1_)[itr.first])[this_index] =
        (*std::get<DataTypeIndex>(particle_data_p2_)[itr.second])[another_index];
}

enum PhaseIndex{lowT=-1,T0=0,highT=1};

void killParticles(BaseParticles &p, size_t n_particles_to_kill_, StdLargeVec<size_t> &particle_kill_list_, StdLargeVec<size_t> &particle_copy_list_, size_t offset = 0){
	if(n_particles_to_kill_<=0) return;
	size_t ntot = p.total_real_particles_;
	size_t remaining = ntot-n_particles_to_kill_;
	size_t killing = 0;
	size_t copying = ntot-1;
	size_t tailkilled = n_particles_to_kill_-1;
	std::sort(particle_kill_list_.begin()+offset,particle_kill_list_.begin()+n_particles_to_kill_+offset);
	for (size_t copied=0; particle_kill_list_[copied]<remaining && copied<n_particles_to_kill_;)
	{
		if(particle_kill_list_[tailkilled]==copying){
			tailkilled--;
			copying--;
			continue;
		}
		particle_copy_list_[copied++] = copying--;
	}
	
	parallel_for(
		blocked_range<size_t>(0, n_particles_to_kill_),
		[&](const blocked_range<size_t> &r)
		{
			for (size_t i = r.begin(); i != r.end(); ++i)
			{
				p.copyFromAnotherParticle(particle_kill_list_[i],particle_copy_list_[i]);
			}
		},
		ap
	);
}

template<class Phase1Particles, class Phase2Particles>
class BasePhaseTransitionDynamics : public InteractionDynamics{
	Phase1Particles &p1_;
	StdLargeVec<Real> &temperature_;
	Phase2Particles &p2_;
	StdLargeVec<Real> &temperature_p2_;
	Real transfer_T_;
	Real latent_heat_;
	StdLargeVec<Real> latent_heat_storage_;
	std::atomic_size_t n_particles_to_kill_;
	std::atomic_size_t n_particles_p2_;
	StdLargeVec<size_t> &particle_kill_list_;
	StdLargeVec<size_t> &particle_copy_list_;
	StdVec<Vecd> displacements;
	char* transitionParName;
	std::unique_ptr<DualParticleDataOperation<copyAParticleDataValueFromOtherParticles>> copy_a_particle_value_;

public:
	BasePhaseTransitionDynamics(Phase1Particles &p1, Phase2Particles &p2, StdLargeVec<size_t> &particle_kill_list, StdLargeVec<size_t> &particle_copy_list, libconfig::Setting* phaseTransitionConfig, char* transitionParName = "Temperature") : 
	InteractionDynamics(*p1.getSPHBody()), p1_(p1), p2_(p2),
	transitionParName(transitionParName),
	temperature_((*(std::get<indexScalar>(p1.all_particle_data_)[p1.all_variable_maps_[indexScalar][transitionParName]]))),
	temperature_p2_((*(std::get<indexScalar>(p2.all_particle_data_)[p2.all_variable_maps_[indexScalar][transitionParName]]))),
	n_particles_p2_(p2.total_real_particles_),
	particle_kill_list_(particle_kill_list),
	particle_copy_list_(particle_copy_list),
	splitting(1)
	{
		phaseTransitionConfig->lookupValue("PhaseChangeTemperature",transfer_T_);
		phaseTransitionConfig->lookupValue("LatentHeat",latent_heat_);
		latent_heat_storage_.resize(p1.real_particles_bound_, Real(0));
		p1.template registerAVariable<indexScalar, Real>(latent_heat_storage_, "LatentHeatStorage");
		copy_a_particle_value_ = make_unique<DualParticleDataOperation<copyAParticleDataValueFromOtherParticles>>(p2,p1);
	}
	void initializeDisplacements(){
		Real d = this->sph_adaptation_->ReferenceSpacing();
		displacements.resize(pow(splitting,Dimensions),Vecd(0.5*d/splitting-0.5*d));
		for(int i = 0; i<splitting ; i++){
			displacements[i][0] += i*d/splitting;
		}
		int length = splitting;
		for(int N = 1; N < Dimensions; N++){
			for(int i = 1; i<splitting ; i++){
				for(int j = 0; j<length ; j++){
					displacements[length*i+j] = displacements[j];
					displacements[length*i+j][N] += i*d/splitting;
				}
			}
			length *= splitting;
		}
	}
	size_t particleTransfer(size_t index_i, Real newT);
	void updateList(){
		if(n_particles_to_kill_.load()) killParticles(p1_, n_particles_to_kill_.load(), particle_kill_list_, particle_copy_list_);
		p1_.total_real_particles_ -= n_particles_to_kill_;
		p2_.total_real_particles_ = n_particles_p2_;
		n_particles_to_kill_ = 0;
	}
	void exec(Real dt = 0){
		n_particles_p2_ = p2_.total_real_particles_;
		InteractionDynamics::exec(dt);
		updateList();
	}
	void parallel_exec(Real dt = 0){
		n_particles_p2_ = p2_.total_real_particles_;
		InteractionDynamics::parallel_exec(dt);
		updateList();
	}
protected:
	int splitting;
	void processLatentHeat(size_t index_i, PhaseIndex pi);
};

template<class Phase1Particles, class Phase2Particles>
void BasePhaseTransitionDynamics<Phase1Particles,Phase2Particles>::processLatentHeat(size_t index_i, PhaseIndex pi){
	Real dT = (transfer_T_-temperature_[index_i])*pi;
	if(dT<0){
		if(latent_heat_storage_[index_i]==0){
			return;
		}else if(latent_heat_storage_[index_i]+dT>0){
			latent_heat_storage_[index_i]+=dT;
			temperature_[index_i] = transfer_T_;
		}else {
			temperature_[index_i]-=pi*latent_heat_storage_[index_i];
			latent_heat_storage_[index_i] = 0;
		}
	}else{
		if((latent_heat_storage_[index_i]+=dT) >= latent_heat_){
			Real extra = latent_heat_storage_[index_i]-latent_heat_;
			latent_heat_storage_[index_i] = 0;
			particleTransfer(index_i, transfer_T_-pi*extra);
		}else{
			temperature_[index_i] = transfer_T_;
		}
	}
}

template<class Phase1Particles, class Phase2Particles>
class phaseTransitionDynamicsLowT: public BasePhaseTransitionDynamics<Phase1Particles,Phase2Particles>{
public:
	phaseTransitionDynamicsLowT(Phase1Particles &p1, Phase2Particles &p2, StdLargeVec<size_t> &particle_kill_list_, StdLargeVec<size_t> &particle_copy_list_, libconfig::Setting* phaseTransitionConfig):
	BasePhaseTransitionDynamics<Phase1Particles,Phase2Particles>(p1, p2, particle_kill_list_, particle_copy_list_, phaseTransitionConfig){
		phaseTransitionConfig->lookupValue("splittingLowT",this->splitting);
		if(this->splitting > 1) this->initializeDisplacements();
	}
protected:
	virtual void Interaction(size_t index_i, Real dt = 0.0) override {
		this->processLatentHeat(index_i,lowT);
	}
};

template<class Phase1Particles, class Phase2Particles>
class phaseTransitionDynamicsHighT: public BasePhaseTransitionDynamics<Phase2Particles,Phase1Particles>{
public:
	phaseTransitionDynamicsHighT(Phase1Particles &p1, Phase2Particles &p2, StdLargeVec<size_t> &particle_kill_list, StdLargeVec<size_t> &particle_copy_list, libconfig::Setting* phaseTransitionConfig):
	BasePhaseTransitionDynamics<Phase2Particles,Phase1Particles>(p2, p1, particle_kill_list, particle_copy_list, phaseTransitionConfig){
		phaseTransitionConfig->lookupValue("splittingHighT",this->splitting);
		if(this->splitting > 1) this->initializeDisplacements();
	}
protected:
	virtual void Interaction(size_t index_i, Real dt = 0.0) override {
		this->processLatentHeat(index_i,highT);
	};
};

template<class Phase1Particles, class Phase2Particles>
size_t BasePhaseTransitionDynamics<Phase1Particles,Phase2Particles>::particleTransfer(size_t index_i, Real newT){
	int splitN = pow(splitting,Dimensions);
	size_t ret = n_particles_p2_.fetch_add(splitN);
	for(int i = 0; i< splitN; i++){
		size_t target = ret+i;
		(*copy_a_particle_value_)(target, index_i);
		p2_.mass_[target] /= splitN;
		temperature_p2_[target] = newT;
		if(splitN>1) p2_.pos_n_[target] += displacements[i];
	}
	particle_kill_list_[n_particles_to_kill_++] = index_i;
	return ret;
}

template<class Phase1Particles, class Phase2Particles>
class phaseTransitionDynamics{
public:
	phaseTransitionDynamicsLowT<Phase1Particles,Phase2Particles> d1;
	phaseTransitionDynamicsHighT<Phase1Particles,Phase2Particles> d2;
	StdLargeVec<size_t> particle_kill_list_;
	StdLargeVec<size_t> particle_copy_list_;
	phaseTransitionDynamics(Phase1Particles &p1, Phase2Particles &p2, libconfig::Setting* phaseTransitionConfig):
	d1(p1, p2, particle_kill_list_, particle_copy_list_, phaseTransitionConfig), d2(p1, p2, particle_kill_list_, particle_copy_list_, phaseTransitionConfig){
		size_t total_particles_p1 = p1.total_real_particles_;
		p1.addBufferParticles(p2.total_real_particles_);
		p2.addBufferParticles(total_particles_p1);
		particle_kill_list_.resize(p1.real_particles_bound_);
		particle_copy_list_.resize(p1.real_particles_bound_/2);
	}//Create cache owned only by phase transion dynamics, not recommanded.
	phaseTransitionDynamics(Phase1Particles &p1, Phase2Particles &p2, libconfig::Setting* phaseTransitionConfig, StdLargeVec<size_t> &cache1, StdLargeVec<size_t> &cache2):
	d1(p1, p2, cache1, cache2, phaseTransitionConfig), d2(p1, p2, cache1, cache2, phaseTransitionConfig)
	{
		size_t total_particles_p1 = p1.total_real_particles_;
		p1.addBufferParticles(p2.total_real_particles_);
		p2.addBufferParticles(total_particles_p1);
		if(cache1.size()<p1.real_particles_bound_){
			printf("Warning: cache1 size %d less than the sum of particle bound %d, this may cause crash when all particles do phase transition at the same time\n", cache1.size(), p1.real_particles_bound_);
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
		}
		if(cache2.size()<p1.real_particles_bound_){
			printf("Warning: cache2 size %d less than the sum of particle bound %d, this may cause crash when all particles do phase transition at the same time\n", cache2.size(), p1.real_particles_bound_);
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
		}
	}//The caches can be shared by other dynamics, recommanded.
	void exec(Real dt = 0){
		d1.exec();
		d2.exec();
	}
	void parallel_exec(Real dt = 0){
		d1.parallel_exec();
		d2.parallel_exec();
	}
};
