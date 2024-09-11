#include "InitializeParameters.hpp"

void ParseXMLFile(const std::string &filename, Parameters &argParameters)
{
    tinyxml2::XMLDocument doc;
    tinyxml2::XMLError eResult = doc.LoadFile(filename.c_str());

    if (eResult != tinyxml2::XML_SUCCESS)
    {
        std::cerr << "Error: Unable to load XML file: " << filename << std::endl;
        return;
    }

    // Parse <simulation> element
    tinyxml2::XMLElement *simulationElement = doc.FirstChildElement("configuration")->FirstChildElement("simulation");
    if (simulationElement)
    {
        const char *scenario_str = simulationElement->Attribute("scenario");
        const char *bc_type_str = simulationElement->Attribute("bcType");
        const char *density_type_Str = simulationElement->Attribute("densityType");
        const char *writer_str = simulationElement->Attribute("writerType");

        if (strcmp(scenario_str, "Poiseuille") == 0)
            argParameters.SCENARIO = POISEUILLE;
        else if (strcmp(scenario_str, "Couette") == 0)
            argParameters.SCENARIO = COUETTE;
        else if (strcmp(scenario_str, "Hydrostatic") == 0)
            argParameters.SCENARIO = HYDROSTATIC;
        else if (strcmp(scenario_str, "CylinderArray") == 0)
            argParameters.SCENARIO = CYLINDER_ARRAY;
        else if (strcmp(scenario_str, "CylinderLattice") == 0)
            argParameters.SCENARIO = CYLINDER_LATTICE;
        else if (strcmp(scenario_str, "Square") == 0)
            argParameters.SCENARIO = SQUARE;
        else if (strcmp(scenario_str, "Triangle") == 0)
            argParameters.SCENARIO = TRIANGLE;
        else if (strcmp(scenario_str, "TriangleEquilateral") == 0)
            argParameters.SCENARIO = TRIANGLE_EQUILATERAL;
        else if (strcmp(scenario_str, "Cavity") == 0)
            argParameters.SCENARIO = CAVITY;
        else if (strcmp(scenario_str, "Step") == 0)
            argParameters.SCENARIO = STEP;
        else if (strcmp(scenario_str, "TaylorCouette") == 0)
            argParameters.SCENARIO = TAYLOR_COUETTE;
        else if (strcmp(scenario_str, "MovingObstacle") == 0)
            argParameters.SCENARIO = MOVING_OBSTACLE;
        else if (strcmp(scenario_str, "Ellipse") == 0)
            argParameters.SCENARIO = ELLIPSE;
        else
        {
            throw std::runtime_error("Unknown scenario");
        }

        argParameters.filename = scenario_str;

        if (strcmp(bc_type_str, "old") == 0)
            argParameters.BC_TYPE = NO_SLIP;
        else if (strcmp(bc_type_str, "new") == 0)
            argParameters.BC_TYPE = NEW_NO_SLIP;

        argParameters.filename += "_";
        argParameters.filename += bc_type_str;

        if (strcmp(density_type_Str, "Summation") == 0)
            argParameters.DENSITY_TYPE = DENSITY_SUMMATION;
        else if (strcmp(density_type_Str, "Differential") == 0)
            argParameters.DENSITY_TYPE = DENSITY_DIFFERENTIAL;

        argParameters.filename += "_";
        argParameters.filename += density_type_Str;

        if (strcmp(writer_str, "VTK") == 0)
            argParameters.WRITER = VTK_WRITER;
        else if (strcmp(writer_str, "CSV") == 0)
            argParameters.WRITER = CSV_WRITER;

        simulationElement->QueryDoubleAttribute("time", &(argParameters.t_end));
        simulationElement->QueryDoubleAttribute("write_const", &(argParameters.write_const));
        const char *custom_str = simulationElement->Attribute("custom_string");

        if (strcmp(custom_str, "") != 0)
        {
            argParameters.filename += "_";
            argParameters.filename += custom_str;
        }
    }

    // Parse <geometry> element
    tinyxml2::XMLElement *geometryElement = doc.FirstChildElement("configuration")->FirstChildElement("geometry");
    if (geometryElement)
    {
        geometryElement->QueryDoubleAttribute("lengthScale", &argParameters.LengthScale);
        geometryElement->QueryDoubleAttribute("rf", &argParameters.rf);
        geometryElement->QueryAttribute("sizeX", &(argParameters.Nfluid[0]));
        geometryElement->QueryAttribute("sizeY", &(argParameters.Nfluid[1]));

        if constexpr (DIM == 3)
            geometryElement->QueryAttribute("sizeZ", &argParameters.Nfluid[2]);
    }

    argParameters.filename += "_" + std::to_string(argParameters.Nfluid[0]) + "_" + std::to_string(argParameters.Nfluid[1]) + "_" + std::to_string((int)argParameters.rf) + "rf";

    // Parse <physics> element
    tinyxml2::XMLElement *physicsElement = doc.FirstChildElement("configuration")->FirstChildElement("physics");
    if (physicsElement)
    {

        physicsElement->QueryDoubleAttribute("rho0", &argParameters.rho_zero);
        physicsElement->QueryDoubleAttribute("nu", &argParameters.nu);
        physicsElement->QueryDoubleAttribute("Bfactor", &argParameters.Bfactor);
        physicsElement->QueryDoubleAttribute("gamma", &argParameters.gamma_);
        physicsElement->QueryDoubleAttribute("umax", &argParameters.umax);
        physicsElement->QueryDoubleAttribute("gravityX", &(argParameters.gravity_vector.get(0)));
        physicsElement->QueryDoubleAttribute("gravityY", &(argParameters.gravity_vector.get(1)));
        if constexpr (DIM == 3)
            physicsElement->QueryDoubleAttribute("gravityZ", &(argParameters.gravity_vector.get(2)));
        physicsElement->QueryDoubleAttribute("vTopX", &(argParameters.vw_top.get(0)));
        physicsElement->QueryDoubleAttribute("vTopY", &(argParameters.vw_top.get(1)));
        if constexpr (DIM == 3)
            physicsElement->QueryDoubleAttribute("vTopZ", &(argParameters.vw_top.get(2)));

        physicsElement->QueryDoubleAttribute("vBotX", &(argParameters.vw_bottom.get(0)));
        physicsElement->QueryDoubleAttribute("vBotY", &(argParameters.vw_bottom.get(1)));
        if constexpr (DIM == 3)
            physicsElement->QueryDoubleAttribute("vBotZ", &(argParameters.vw_bottom.get(2)));
    }

    // Parse <obstacle> element
    tinyxml2::XMLElement *obstacleElement = doc.FirstChildElement("configuration")->FirstChildElement("obstacle");
    if (obstacleElement)
    {

        // optional parameters ObstacleBase, ObstacleHeight, tilt
        argParameters.ObstacleBase = 0.0;
        argParameters.ObstacleHeight = 0.0;
        const char *obstacleBaseStr = obstacleElement->Attribute("ObstacleBase");
        if (obstacleBaseStr)
        {
            argParameters.ObstacleBase = std::stod(obstacleBaseStr); // Convert to float
        }
        const char *obstacleHeightStr = obstacleElement->Attribute("ObstacleHeight");
        if (obstacleHeightStr)
        {
            argParameters.ObstacleHeight = std::stod(obstacleHeightStr); // Convert to float
        }
        const char *tiltStr = obstacleElement->Attribute("tilt");
        if (tiltStr)
        {
            argParameters.ObstacleTilt = std::stod(tiltStr); // Convert to float
        }

        obstacleElement->QueryDoubleAttribute("velX", &(argParameters.ObstacleVelocity.get(0)));
        obstacleElement->QueryDoubleAttribute("velY", &(argParameters.ObstacleVelocity.get(1)));
        if constexpr (DIM == 3)
            obstacleElement->QueryDoubleAttribute("velZ", &(argParameters.ObstacleVelocity.get(2)));

        obstacleElement->QueryDoubleAttribute("omega", &argParameters.ObstacleOmega);
    }

    // Parse <TaylorCouette> element
    tinyxml2::XMLElement *taylorCouetteElement = doc.FirstChildElement("configuration")->FirstChildElement("TaylorCouette");
    if (taylorCouetteElement)
    {
        taylorCouetteElement->QueryDoubleAttribute("Rin", &argParameters.Rin);
        taylorCouetteElement->QueryDoubleAttribute("Rout", &argParameters.Rout);
        taylorCouetteElement->QueryDoubleAttribute("Win", &argParameters.Win);
        taylorCouetteElement->QueryDoubleAttribute("Wout", &argParameters.Wout);
    }
}
void InitializeConstants(Vcluster<> &v_cl, Parameters &argParameters)
{
    // Given the scenario and xml parameters, we set the constants of argParameters

    // Set boundary conditions periodic or non periodic
    size_t Nbound = (argParameters.BC_TYPE == NEW_NO_SLIP) ? 1 : 3;

    if (argParameters.SCENARIO == POISEUILLE ||
        argParameters.SCENARIO == COUETTE ||
        argParameters.SCENARIO == CYLINDER_ARRAY ||
        argParameters.SCENARIO == SQUARE ||
        argParameters.SCENARIO == TRIANGLE ||
        argParameters.SCENARIO == TRIANGLE_EQUILATERAL)
    {
        argParameters.Nboundary[0] = 0;
        argParameters.Nboundary[1] = Nbound;
        argParameters.bc[0] = PERIODIC;
        argParameters.bc[1] = NON_PERIODIC;
    }
    else if (argParameters.SCENARIO == HYDROSTATIC ||
             argParameters.SCENARIO == CAVITY ||
             argParameters.SCENARIO == MOVING_OBSTACLE ||
             argParameters.SCENARIO == TAYLOR_COUETTE)
    {
        argParameters.Nboundary[0] = Nbound;
        argParameters.Nboundary[1] = Nbound;
        argParameters.bc[0] = NON_PERIODIC;
        argParameters.bc[1] = NON_PERIODIC;
    }
    else if (argParameters.SCENARIO == CYLINDER_LATTICE || argParameters.SCENARIO == ELLIPSE)
    {
        argParameters.Nboundary[0] = 0;
        argParameters.Nboundary[1] = 0;
        argParameters.bc[0] = PERIODIC;
        argParameters.bc[1] = PERIODIC;
    }

    // Set particle spacing, definition depends on scenario
    if (argParameters.SCENARIO != CYLINDER_ARRAY && argParameters.SCENARIO != CYLINDER_LATTICE)
    {
        argParameters.dp = argParameters.LengthScale / argParameters.Nfluid[1];
    }
    else
    {
        if (argParameters.SCENARIO == CYLINDER_ARRAY) // chanel height is 4 times the cylinder radius
            argParameters.dp = 4.0 * argParameters.LengthScale / (double)argParameters.Nfluid[1];
        else if (argParameters.SCENARIO == CYLINDER_LATTICE) // chanel height is 5 times the cylinder radius
            argParameters.dp = 5.0 * argParameters.LengthScale / (double)argParameters.Nfluid[1];
    }

    argParameters.length[0] = argParameters.dp * (argParameters.Nfluid[0] - argParameters.bc[0]);
    argParameters.length[1] = argParameters.dp * (argParameters.Nfluid[1] - argParameters.bc[1]);

    // WE NO LONGER SET THE MAXIMUM VELOCITY; WE READ IT FROM THE XML FILE, BUT THIS VALUES ARE USEFUL
    // // Set maximum velocity
    // if (argParameters.SCENARIO == POISEUILLE)
    //     argParameters.umax = argParameters.gravity_vector.get(0) * argParameters.LengthScale * argParameters.LengthScale / (8.0 * argParameters.nu);
    // else if (argParameters.SCENARIO == COUETTE)
    //     argParameters.umax = argParameters.vw_top.get(0);
    // else if (argParameters.SCENARIO == HYDROSTATIC)
    //     argParameters.umax = 0.1;
    // else if (argParameters.SCENARIO == CYLINDER_ARRAY)
    //     argParameters.umax = 2.5e-3;
    // else if (argParameters.SCENARIO == CYLINDER_LATTICE)
    //     argParameters.umax = 1.2e-4; // umax = 5.77 * 1e-5; // (morris, to get c=5.77*1e-4)
    // else if (argParameters.SCENARIO == SQUARE)
    //     argParameters.umax = 4.1 * 1e-1;
    // else if (argParameters.SCENARIO == TRIANGLE)
    //     argParameters.umax = 3.5 * 1e-1; // from previous simulations for nu = 0.01
    // else if (argParameters.SCENARIO == TRIANGLE_EQUILATERAL)
    //     argParameters.umax = 4e-1; // from previous simulations for nu = 0.01
    // else if (argParameters.SCENARIO == CAVITY)
    //     argParameters.umax = argParameters.vw_top.get(0);
    // else if (argParameters.SCENARIO == MOVING_OBSTACLE)
    //     argParameters.umax = 1.5; // from previous simulations for nu = 0.01
    // else if (argParameters.SCENARIO == TAYLOR_COUETTE)
    //     argParameters.umax = 4.0 * abs(argParameters.Wout * argParameters.Rout);

    // Enable probes in some scenarios
    if (argParameters.SCENARIO == CAVITY || argParameters.SCENARIO == CYLINDER_LATTICE)
    {
        argParameters.PROBES_ENABLED = 1;
    }
    // Set general parameters
    argParameters.H = argParameters.Hconst * argParameters.dp;
    argParameters.r_threshold = 3.0 * argParameters.H;
    argParameters.Kquintic = (DIM == 3) ? 1.0 / 120.0 / M_PI / argParameters.H / argParameters.H / argParameters.H : 7.0 / 478.0 / M_PI / argParameters.H / argParameters.H;
    argParameters.MassFluid = argParameters.rho_zero * (DIM == 3 ? argParameters.dp * argParameters.dp * argParameters.dp : argParameters.dp * argParameters.dp);
    argParameters.MassBound = argParameters.rho_zero * (DIM == 3 ? argParameters.dp * argParameters.dp * argParameters.dp : argParameters.dp * argParameters.dp);
    argParameters.cbar = argParameters.coeff_sound * argParameters.umax;
    argParameters.B = argParameters.rho_zero * argParameters.cbar * argParameters.cbar / argParameters.gamma_;
    argParameters.Pbackground = argParameters.Bfactor * argParameters.B;
    argParameters.eta = argParameters.nu * argParameters.rho_zero;
    argParameters.Re = argParameters.umax * argParameters.LengthScale / argParameters.nu;
    argParameters.gravity = getVectorNorm(argParameters.gravity_vector);
    argParameters.ObstacleCenter = {(argParameters.length[0] + argParameters.bc[0] * argParameters.dp) / 2.0, (argParameters.length[1] + argParameters.bc[1] * argParameters.dp) / 2.0};

    // add the number of processors to the filename
    const long int Nprc = v_cl.getProcessingUnits();
    argParameters.filename += "_" + std::to_string(Nprc) + "prc";
}

// void SetFilename(std::string &filename, const long int Nprc, const Parameters &params)
// {
//     filename = params.filename;
//     // Scenario name
//     if (params.SCENARIO == POISEUILLE)
//         filename = "Poiseuille";
//     else if (params.SCENARIO == COUETTE)
//         filename = "Couette";
//     else if (params.SCENARIO == HYDROSTATIC)
//         filename = "Hydrostatic";
//     else if (params.SCENARIO == CYLINDER_ARRAY)
//         filename = "CylinderArray";
//     else if (params.SCENARIO == CYLINDER_LATTICE)
//         filename = "CylinderLattice";
//     else if (params.SCENARIO == SQUARE)
//         filename = "Square";
//     else if (params.SCENARIO == TRIANGLE)
//         filename = "Triangle";
//     else if (params.SCENARIO == CAVITY)
//         filename = "Cavity";
//     else if (params.SCENARIO == STEP)
//         filename = "Step";
//     else if (params.SCENARIO == TAYLOR_COUETTE)
//         filename = "TaylorCouette";
//     else if (params.SCENARIO == TRIANGLE_EQUILATERAL)
//         filename = "TriangleEquilateral";
//     else if (params.SCENARIO == MOVING_OBSTACLE)
//         filename = "MovingObstacle";

//     // BC name
//     if (params.BC_TYPE == NO_SLIP)
//         filename += "_OLD_BC";
//     else if (params.BC_TYPE == NEW_NO_SLIP)
//         filename += "_NEW_BC";

//     // Density name
//     if (params.DENSITY_TYPE == DENSITY_SUMMATION)
//         filename += "_Summation";
//     else if (params.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
//         filename += "_Differential";

//     // Add the size of the simulation, the number of processors and the refinement factor
//     std::string size_proc_name = std::to_string(params.Nfluid[0]) + "_" + std::to_string(params.Nfluid[1]) + "_" + std::to_string(Nprc) + "prc";
//     if (params.BC_TYPE == NEW_NO_SLIP)
//         size_proc_name += "_" + std::to_string((int)params.rf) + "rf";
//     filename += "_" + size_proc_name;
//     filename += params.custom_string;
// }
