#ifndef INITIALIZEPARAMETERS_HPP
#define INITIALIZEPARAMETERS_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include <tinyxml2.h>

void ParseXMLFile(const std::string &filename, Parameters &argParameters, AuxiliarParameters &argFilenames);

void InitializeConstants( Parameters &argParameters, AuxiliarParameters &argFilenames);

#endif // INITIALIZEPARAMETERS_HPP
