#include "InitializeParameters.hpp"

void ParseXMLFile(const std::string &filename, Parameters &argParameters, AuxiliarParameters &argAuxParameters)
{
    tinyxml2::XMLDocument doc;
    tinyxml2::XMLError eResult = doc.LoadFile(filename.c_str());

    if (eResult != tinyxml2::XML_SUCCESS)
    {
        std::cerr << "Error: Unable to load XML file: " << filename << std::endl;
        return;
    }
    // to read double values and later cast to real_number
    double tmpdouble;

    // Parse <simulation> element
    tinyxml2::XMLElement *simulationElement = doc.FirstChildElement("configuration")->FirstChildElement("simulation");
    if (simulationElement)
    {

        // read scenario
        const char *scenario_str = simulationElement->Attribute("scenario");
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
        else if (strcmp(scenario_str, "DamBreak") == 0)
            argParameters.SCENARIO = DAM_BREAK;
        else if (strcmp(scenario_str, "TriangleTest") == 0)
            argParameters.SCENARIO = TRIANGLE_TEST;
        else
        {
            std::cout << "Unknown Scenario, defaulting to Poiseuille" << std::endl;
            argParameters.SCENARIO = POISEUILLE;
        }
        argAuxParameters.filename = scenario_str;

        // read BC type (old or new)
        const char *bc_type_str = simulationElement->Attribute("bcType");
        if (strcmp(bc_type_str, "old") == 0)
            argParameters.BC_TYPE = NO_SLIP;
        else if (strcmp(bc_type_str, "new") == 0)
            argParameters.BC_TYPE = NEW_NO_SLIP;

        argAuxParameters.filename += "_";
        argAuxParameters.filename += bc_type_str;

        // read BC periodic or non periodic
        const char *bcX = simulationElement->Attribute("bcX");
        if (strcmp(bcX, "periodic") == 0)
            argParameters.bc[0] = PERIODIC;
        else if (strcmp(bcX, "non_periodic") == 0)
            argParameters.bc[0] = NON_PERIODIC;

        const char *bcY = simulationElement->Attribute("bcY");
        if (strcmp(bcY, "periodic") == 0)
            argParameters.bc[1] = PERIODIC;
        else if (strcmp(bcY, "non_periodic") == 0)
            argParameters.bc[1] = NON_PERIODIC;

        if constexpr (DIM == 3)
        {
            const char *bcZ = simulationElement->Attribute("bcZ");
            if (strcmp(bcZ, "periodic") == 0)
                argParameters.bc[2] = PERIODIC;
            else if (strcmp(bcZ, "non_periodic") == 0)
                argParameters.bc[2] = NON_PERIODIC;
        }

        // read density type
        const char *density_type_Str = simulationElement->Attribute("densityType");
        if (strcmp(density_type_Str, "summation") == 0)
            argParameters.DENSITY_TYPE = DENSITY_SUMMATION;
        else if (strcmp(density_type_Str, "differential") == 0)
            argParameters.DENSITY_TYPE = DENSITY_DIFFERENTIAL;
        else
        {
            std::cout << "Unknown density type, using default: summation" << std::endl;
            argParameters.DENSITY_TYPE = DENSITY_SUMMATION;
        }

        argAuxParameters.filename += "_";
        argAuxParameters.filename += density_type_Str;

        // read writer type
        const char *writer_str = simulationElement->Attribute("writerType");
        if (strcmp(writer_str, "VTK") == 0)
            argParameters.WRITER = VTK_WRITER;
        else if (strcmp(writer_str, "CSV") == 0)
            argParameters.WRITER = CSV_WRITER;

        // read time parameters
        simulationElement->QueryDoubleAttribute("time", &tmpdouble);
        argParameters.t_end = static_cast<real_number>(tmpdouble);
        simulationElement->QueryDoubleAttribute("write_const", &tmpdouble);
        argParameters.write_const = static_cast<real_number>(tmpdouble);
        simulationElement->QueryDoubleAttribute("CFL", &tmpdouble);
        argParameters.CFLnumber = static_cast<real_number>(tmpdouble);

        // read custom string (optional) for custom output names
        const char *custom_str = simulationElement->Attribute("custom_string");
        if (strcmp(custom_str, "") != 0)
        {
            argAuxParameters.filename += "_";
            argAuxParameters.filename += custom_str;
        }
    }

    // Parse <geometry> element
    tinyxml2::XMLElement *geometryElement = doc.FirstChildElement("configuration")->FirstChildElement("geometry");
    if (geometryElement)
    {

        geometryElement->QueryDoubleAttribute("lengthScale", &tmpdouble);
        argParameters.LengthScale = static_cast<real_number>(tmpdouble);
        geometryElement->QueryDoubleAttribute("rf", &tmpdouble);
        argParameters.rf = static_cast<real_number>(tmpdouble);
        geometryElement->QueryAttribute("sizeX", &(argParameters.Nfluid[0]));
        geometryElement->QueryAttribute("sizeY", &(argParameters.Nfluid[1]));

        if constexpr (DIM == 3)
            geometryElement->QueryAttribute("sizeZ", &argParameters.Nfluid[2]);
    }

    argAuxParameters.filename += "_" + std::to_string(argParameters.Nfluid[0]) + "_" + std::to_string(argParameters.Nfluid[1]) + "_" + std::to_string((int)argParameters.rf) + "rf";

    // Parse <physics> element
    tinyxml2::XMLElement *physicsElement = doc.FirstChildElement("configuration")->FirstChildElement("physics");
    if (physicsElement)
    {

        physicsElement->QueryDoubleAttribute("rho0", &tmpdouble);
        argParameters.rho0 = static_cast<real_number>(tmpdouble);
        physicsElement->QueryDoubleAttribute("nu", &tmpdouble);
        argParameters.nu = static_cast<real_number>(tmpdouble);
        physicsElement->QueryDoubleAttribute("Bfactor", &tmpdouble);
        argParameters.Bfactor = static_cast<real_number>(tmpdouble);
        physicsElement->QueryDoubleAttribute("gamma", &tmpdouble);
        argParameters.gamma = static_cast<real_number>(tmpdouble);
        physicsElement->QueryDoubleAttribute("umax", &tmpdouble);
        argParameters.umax = static_cast<real_number>(tmpdouble);
        physicsElement->QueryDoubleAttribute("gravityX", &tmpdouble);
        argParameters.gravity_vector[0] = static_cast<real_number>(tmpdouble);
        physicsElement->QueryDoubleAttribute("gravityY", &tmpdouble);
        argParameters.gravity_vector[1] = static_cast<real_number>(tmpdouble);
        if constexpr (DIM == 3)
        {
            physicsElement->QueryDoubleAttribute("gravityZ", &tmpdouble);
            argParameters.gravity_vector[2] = static_cast<real_number>(tmpdouble);
        }

        physicsElement->QueryDoubleAttribute("vTopX", &tmpdouble);
        argParameters.vw_top[0] = static_cast<real_number>(tmpdouble);
        physicsElement->QueryDoubleAttribute("vTopY", &tmpdouble);
        argParameters.vw_top[1] = static_cast<real_number>(tmpdouble);
        if constexpr (DIM == 3)
        {
            physicsElement->QueryDoubleAttribute("vTopZ", &tmpdouble);
            argParameters.vw_top[2] = static_cast<real_number>(tmpdouble);
        }

        physicsElement->QueryDoubleAttribute("vBotX", &tmpdouble);
        argParameters.vw_bottom[0] = static_cast<real_number>(tmpdouble);
        physicsElement->QueryDoubleAttribute("vBotY", &tmpdouble);
        argParameters.vw_bottom[1] = static_cast<real_number>(tmpdouble);
        if constexpr (DIM == 3)
        {
            physicsElement->QueryDoubleAttribute("vBotZ", &tmpdouble);
            argParameters.vw_bottom[2] = static_cast<real_number>(tmpdouble);
        }
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

        obstacleElement->QueryDoubleAttribute("velX", &tmpdouble);
        argParameters.ObstacleVelocity[0] = static_cast<real_number>(tmpdouble);
        obstacleElement->QueryDoubleAttribute("velY", &tmpdouble);
        argParameters.ObstacleVelocity[1] = static_cast<real_number>(tmpdouble);
        if constexpr (DIM == 3)
        {
            obstacleElement->QueryDoubleAttribute("velZ", &tmpdouble);
            argParameters.ObstacleVelocity[2] = static_cast<real_number>(tmpdouble);
        }

        obstacleElement->QueryDoubleAttribute("omega", &tmpdouble);
        argParameters.ObstacleOmega = static_cast<real_number>(tmpdouble);
    }

    // Parse <TaylorCouette> element
    tinyxml2::XMLElement *taylorCouetteElement = doc.FirstChildElement("configuration")->FirstChildElement("TaylorCouette");
    if (taylorCouetteElement)
    {
        taylorCouetteElement->QueryDoubleAttribute("Rin", &tmpdouble);
        argParameters.Rin = static_cast<real_number>(tmpdouble);
        taylorCouetteElement->QueryDoubleAttribute("Rout", &tmpdouble);
        argParameters.Rout = static_cast<real_number>(tmpdouble);
        taylorCouetteElement->QueryDoubleAttribute("Win", &tmpdouble);
        argParameters.Win = static_cast<real_number>(tmpdouble);
        taylorCouetteElement->QueryDoubleAttribute("Wout", &tmpdouble);
        argParameters.Wout = static_cast<real_number>(tmpdouble);
    }

    // Parse <DamBreak> element
    tinyxml2::XMLElement *DamBreakElement = doc.FirstChildElement("configuration")->FirstChildElement("dambreak");
    if (DamBreakElement)
    {
        DamBreakElement->QueryDoubleAttribute("WaterHeigth", &tmpdouble);
        argParameters.waterH = static_cast<real_number>(tmpdouble);
        DamBreakElement->QueryDoubleAttribute("WaterBase", &tmpdouble);
        argParameters.waterB = static_cast<real_number>(tmpdouble);
    }
}
void InitializeConstants(Parameters &argParameters, AuxiliarParameters &argAuxParameters)
{
    // Given the scenario and xml parameters, we set the constants of argParameters

    Vcluster<> &v_cl = create_vcluster();

    if (argParameters.SCENARIO == TAYLOR_COUETTE)
    {
        // the input value is the number of particles in the radial direction
        // by multiplying times 8 we get number of particles in the square grid containing the cylinders
        argParameters.Nfluid[0] = argParameters.Nfluid[0] * 8;
        argParameters.Nfluid[1] = argParameters.Nfluid[1] * 8;
    }

    // First set fixed values
    argParameters.Hconst = 1.0;
    argParameters.coeff_sound = 10.0;
    argParameters.xi = 0.0;
    argParameters.PROBES_ENABLED = 0;

    // Set boundary conditions periodic or non periodic
    size_t Nbound = (argParameters.BC_TYPE == NEW_NO_SLIP) ? 1 : 3;

    if (argParameters.bc[0] == PERIODIC)
    {
        if (argParameters.SCENARIO == CAVITY)
        {
            argParameters.Nboundary[0] = Nbound;
        }
        else
        {
            argParameters.Nboundary[0] = 0;
        }
    }
    else
        argParameters.Nboundary[0] = Nbound;

    if (argParameters.bc[1] == PERIODIC)
        argParameters.Nboundary[1] = 0;
    else
        argParameters.Nboundary[1] = Nbound;

    // Set particle spacing, definition depends on scenario, for most scenarios length scale gives height of the box, so dp is height divided by Nfluid[1]
    if (argParameters.SCENARIO != CYLINDER_ARRAY && argParameters.SCENARIO != CYLINDER_LATTICE && argParameters.SCENARIO != TAYLOR_COUETTE)
    {
        argParameters.dp = argParameters.LengthScale / argParameters.Nfluid[1];
    }
    else
    {
        if (argParameters.SCENARIO == CYLINDER_ARRAY) // length scale is cylinder radius, chanel height is 4 times the cylinder radius
            argParameters.dp = 4.0 * argParameters.LengthScale / (real_number)argParameters.Nfluid[1];
        else if (argParameters.SCENARIO == CYLINDER_LATTICE) // length scale is cylinder radius, chanel height is 5 times the cylinder radius
            argParameters.dp = 5.0 * argParameters.LengthScale / (real_number)argParameters.Nfluid[1];
        else if (argParameters.SCENARIO == TAYLOR_COUETTE) // Length scale is inter-cylinder distance
            argParameters.dp = 8.0 * argParameters.LengthScale / argParameters.Nfluid[1];
    }

    if (argParameters.SCENARIO == CAVITY)
    {
        argParameters.length[0] = argParameters.dp * (argParameters.Nfluid[0]);
        argParameters.length[1] = argParameters.dp * (argParameters.Nfluid[1]);
    }
    else
    {
        argParameters.length[0] = argParameters.dp * (argParameters.Nfluid[0] - argParameters.bc[0]);
        argParameters.length[1] = argParameters.dp * (argParameters.Nfluid[1] - argParameters.bc[1]);
    }

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
    argParameters.r_cut = 3.0 * argParameters.H;
    argParameters.r_cut2 = argParameters.r_cut * argParameters.r_cut;
    argParameters.Kquintic = (DIM == 3) ? 1.0 / 120.0 / M_PI / argParameters.H / argParameters.H / argParameters.H : 7.0 / 478.0 / M_PI / argParameters.H / argParameters.H;
    argParameters.MassFluid = argParameters.rho0 * (DIM == 3 ? argParameters.dp * argParameters.dp * argParameters.dp : argParameters.dp * argParameters.dp);
    argParameters.MassBound = argParameters.rho0 * (DIM == 3 ? argParameters.dp * argParameters.dp * argParameters.dp : argParameters.dp * argParameters.dp);
    argParameters.gravity = getVectorNorm(argParameters.gravity_vector);
    real_number c, cu, cnu, cg;
    cu = argParameters.coeff_sound * argParameters.umax;
    cnu = argParameters.coeff_sound * sqrt(argParameters.nu * argParameters.umax / argParameters.LengthScale);
    cg = argParameters.coeff_sound * sqrt(argParameters.gravity * argParameters.LengthScale);
    c = max(max(cu, cnu), cg);
    argParameters.cbar = c;
    if (c == cu)
    {
        std::cout << "Sound speed is set by umax c = " << c << std::endl;
    }
    else if (c == cnu)
    {
        std::cout << "Sound speed is set by nu c = " << c << std::endl;
    }
    else if (c == cg)
    {
        std::cout << "Sound speed is set by gravity c = " << c << std::endl;
    }

    // argParameters.cbar = argParameters.coeff_sound * argParameters.umax;
    argParameters.B = argParameters.rho0 * argParameters.cbar * argParameters.cbar / argParameters.gamma;
    argParameters.Pbackground = argParameters.Bfactor * argParameters.B;
    argParameters.eta = argParameters.nu * argParameters.rho0;
    argParameters.Re = argParameters.umax * argParameters.LengthScale / argParameters.nu;
    argParameters.ObstacleCenter[0] = (argParameters.length[0] + (real_number)argParameters.bc[0] * argParameters.dp) / 2.0f;
    argParameters.ObstacleCenter[1] = (argParameters.length[1] + (real_number)argParameters.bc[1] * argParameters.dp) / 2.0f;

    // add the number of processors to the filename
    const long int Nprc = v_cl.getProcessingUnits();
    argAuxParameters.filename += "_" + std::to_string(Nprc) + "prc";
}
