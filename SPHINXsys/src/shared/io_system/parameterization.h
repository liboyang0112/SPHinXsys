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
 * @file 	parameterization.h
 * @brief 	This is the base classes for introducing the parameterization 
 *			of a class or method.
 * @author	Xiangyu Hu
 */

#ifndef SPHINXSYS_PARAMETERIZATION_H
#define SPHINXSYS_PARAMETERIZATION_H

#include "base_data_package.h"
#include "in_output.h"
#include "xml_engine.h"

#include <string>

namespace SPH 
{
	class ParameterizationIO
	{
	public:
		XmlEngine xml_paremeters_;
		std::string filefullpath_;

		ParameterizationIO(In_Output& in_output);
		~ParameterizationIO() {};

		void writeProjectParameters();
	};

	template<class BaseClassType>
	class BaseParameterization : public BaseClassType
	{
	public:
		template<typename... ConstructorArgs>
		explicit BaseParameterization(ParameterizationIO& parameterization_io, ConstructorArgs... constructor_args) :
			BaseClassType(constructor_args...), xml_paremeters_(parameterization_io.xml_paremeters_),
			filefullpath_(parameterization_io.filefullpath_) {};
		~BaseParameterization() {};
	protected:
		XmlEngine& xml_paremeters_;
		std::string filefullpath_;

		template<typename VariableType>
		void getAParameter(const std::string& element_name, const std::string& variable_name, VariableType& variable_addrs)
		{
			SimTK::Xml::element_iterator ele_ite = 
				xml_paremeters_.root_element_.element_begin(element_name);
			if(ele_ite != xml_paremeters_.root_element_.element_end())
			{
				xml_paremeters_.getRequiredAttributeValue(ele_ite, variable_name, variable_addrs);
			}
			else {
				std::cout << "\n Error: the variable '" << variable_name << "' is given not in project_parameters.dat !" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}
		};
	
		template<typename VariableType>
		void setAParameter(const std::string& element_name, const std::string& variable_name, VariableType& variable_addrs)
		{
			SimTK::Xml::element_iterator ele_ite =
				xml_paremeters_.root_element_.element_begin(element_name);
			if (ele_ite == xml_paremeters_.root_element_.element_end())
			{
				xml_paremeters_.addElementToXmlDoc(element_name);
				ele_ite = xml_paremeters_.root_element_.element_begin(element_name);
			}
			xml_paremeters_.setAttributeToElement(ele_ite, variable_name, variable_addrs);
		};
	};
}
#endif //SPHINXSYS_PARAMETERIZATION_H
