#ifndef INITIALIZEPARAMETERS_HPP
#define INITIALIZEPARAMETERS_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include "Kernel.hpp"
#include <tinyxml2.h>

void ParseXMLFile(const std::string &filename, Parameters &argParameters, AuxiliarParameters &argFilenames);

void ComputeKernelVolume(Parameters &argParameters);

void InitializeConstants(Parameters &argParameters, AuxiliarParameters &argFilenames);

#endif // INITIALIZEPARAMETERS_HPP
