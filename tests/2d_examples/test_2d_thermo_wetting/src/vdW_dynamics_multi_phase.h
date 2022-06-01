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
* @file fluid_dynamics_multi_phase.h
* @brief Here, we define the algorithm classes for the dynamics involving multiple fluids.   
* @author	Boyang Li
*/


#ifndef VDW_DYNAMICS_MULTI_PHASE_H
#define VDW_DYNAMICS_MULTI_PHASE_H

#include "/home/boyang/softwares/SPHinXsys/SPHinXsys/SPHINXsys/src/shared/particle_dynamics/fluid_dynamics/fluid_dynamics_multi_phase.hpp"
#include "sphinxsys.h"
namespace SPH
{
	namespace fluid_dynamics
	{
		template<class PressureRelaxationInnerType>
		class vdWPressureRelaxationMultiPhase : public BasePressureRelaxationMultiPhase<PressureRelaxationInnerType>
		{
		public:
			vdWPressureRelaxationMultiPhase(BaseBodyRelationInner &inner_relation,
				BaseBodyRelationContact &contact_relation);
			explicit vdWPressureRelaxationMultiPhase(ComplexBodyRelation &complex_relation);
			virtual ~vdWPressureRelaxationMultiPhase() {};
		protected:
			virtual void Interaction(size_t index_i, Real dt = 0.0) override;
		};
	}
}
#endif //FLUID_DYNAMICS_MULTI_PHASE_H