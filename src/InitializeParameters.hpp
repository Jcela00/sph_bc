#ifndef INITIALIZEPARAMETERS_HPP
#define INITIALIZEPARAMETERS_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"
#include <tinyxml2.h>

void ParseXMLFile(const std::string &filename, Parameters &argParameters);

void InitializeConstants(Vcluster<> &v_cl, Parameters &argParameters);

void WriteParameters(const Parameters &argParameters);

#endif // INITIALIZEPARAMETERS_HPP
