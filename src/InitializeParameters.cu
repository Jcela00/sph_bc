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
        else if (strcmp(scenario_str, "Ellipse") == 0)
            argParameters.SCENARIO = ELLIPSE;
        else if (strcmp(scenario_str, "DamBreak") == 0)
            argParameters.SCENARIO = DAM_BREAK;
        else if (strcmp(scenario_str, "DamBreakAdj") == 0)
            argParameters.SCENARIO = DAM_BREAK_ADJ;
        else if (strcmp(scenario_str, "TriangleTest") == 0)
            argParameters.SCENARIO = TRIANGLE_TEST;
        else if (strcmp(scenario_str, "Sphere") == 0)
            argParameters.SCENARIO = SPHERE;
        else if (strcmp(scenario_str, "Custom") == 0)
            argParameters.SCENARIO = CUSTOM;
        else if (strcmp(scenario_str, "Flower") == 0)
            argParameters.SCENARIO = FLOWER;
        else if (strcmp(scenario_str, "PoiseuilleTank") == 0)
            argParameters.SCENARIO = POISEUILLE_TANK;
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
        else if (strcmp(writer_str, "BINARY") == 0)
            argParameters.WRITER = FORMAT_BINARY;

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
    if constexpr (DIM == 2)
        argAuxParameters.filename += "_" + std::to_string(argParameters.Nfluid[0]) + "_" + std::to_string(argParameters.Nfluid[1]) + "_" + std::to_string((int)argParameters.rf) + "rf";
    else if constexpr (DIM == 3)
        argAuxParameters.filename += "_" + std::to_string(argParameters.Nfluid[0]) + "_" + std::to_string(argParameters.Nfluid[1]) + "_" + std::to_string(argParameters.Nfluid[2]) + "_" + std::to_string((int)argParameters.rf) + "rf";

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

        // optional parameters ObstacleBase, ObstacleHeight, tilt, obstacle centre
        argParameters.ObstacleBase = 0.0;
        argParameters.ObstacleHeight = 0.0;
        argParameters.ObstacleTilt = 0.0;
        argParameters.CustomObstacle = 0;

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
        const char *centerX = obstacleElement->Attribute("centerX");
        const char *centerY = obstacleElement->Attribute("centerY");
        const char *Vinflowx = obstacleElement->Attribute("VinflowX");
        const char *Vinflowy = obstacleElement->Attribute("VinflowY");
        if (centerX && centerY)
        {
            argParameters.ObstacleCenter[0] = std::stod(centerX); // Convert to float
            argParameters.ObstacleCenter[1] = std::stod(centerY); // Convert to float
        }
        if (Vinflowx && Vinflowy)
        {
            argParameters.Vinflow[0] = std::stod(Vinflowx); // Convert to float
            argParameters.Vinflow[1] = std::stod(Vinflowy); // Convert to float
        }
        if constexpr (DIM == 3)
        {
            const char *centerZ = obstacleElement->Attribute("centerZ");
            if (centerZ)
            {
                argParameters.ObstacleCenter[2] = std::stod(centerZ); // Convert to float
            }
            const char *Vinflowz = obstacleElement->Attribute("VinflowZ");
            if (Vinflowz)
            {
                argParameters.Vinflow[2] = std::stod(Vinflowz); // Convert to float
            }
        }

        const char *customObstacleStr = obstacleElement->Attribute("custom");
        if (customObstacleStr)
        {
            if (strcmp(customObstacleStr, "Cylinder") == 0)
            {
                argParameters.CustomObstacle = CYLINDER_LATTICE;
                argAuxParameters.filename += "_Cylinder";
            }
            else if (strcmp(customObstacleStr, "Square") == 0)
            {
                argParameters.CustomObstacle = SQUARE;
                argAuxParameters.filename += "_Square";
            }
            else if (strcmp(customObstacleStr, "Triangle") == 0)
            {
                argParameters.CustomObstacle = TRIANGLE;
                argAuxParameters.filename += "_Triangle";
            }
            else if (strcmp(customObstacleStr, "TriangleEquilateral") == 0)
            {
                argParameters.CustomObstacle = TRIANGLE_EQUILATERAL;
                argAuxParameters.filename += "_TriangleEquilateral";
            }
            else if (strcmp(customObstacleStr, "Ellipse") == 0)
            {
                argParameters.CustomObstacle = ELLIPSE;
                argAuxParameters.filename += "_Ellipse";
            }
            else
            {
                std::cout << "Unknown custom obstacle, defaulting to Ellipse " << std::endl;
                argParameters.CustomObstacle = ELLIPSE;
                argAuxParameters.filename += "_Ellipse";
            }
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
    // Parse <DamBreak> element
    tinyxml2::XMLElement *FlowerElement = doc.FirstChildElement("configuration")->FirstChildElement("flower");
    if (FlowerElement)
    {
        FlowerElement->QueryDoubleAttribute("a", &tmpdouble);
        argParameters.flowerA = static_cast<real_number>(tmpdouble);
        FlowerElement->QueryDoubleAttribute("b", &tmpdouble);
        argParameters.flowerB = static_cast<real_number>(tmpdouble);
        FlowerElement->QueryDoubleAttribute("k", &tmpdouble);
        argParameters.flowerK = static_cast<real_number>(tmpdouble);
        FlowerElement->QueryDoubleAttribute("m", &tmpdouble);
        argParameters.flowerM = static_cast<real_number>(tmpdouble);
        FlowerElement->QueryDoubleAttribute("nlobes", &tmpdouble);
        argParameters.flowerNlobes = static_cast<real_number>(tmpdouble);
        FlowerElement->QueryDoubleAttribute("alternate", &tmpdouble);
        argParameters.flowerAlternate = static_cast<int>(tmpdouble);
        FlowerElement->QueryDoubleAttribute("propRadius", &tmpdouble);
        argParameters.flowerPropRadius = static_cast<real_number>(tmpdouble);
        FlowerElement->QueryDoubleAttribute("propRadius2", &tmpdouble);
        argParameters.flowerPropRadius2 = static_cast<real_number>(tmpdouble);
    }
}

void ComputeKernelVolume(Parameters &argParameters)
{
    // compute the kernel volume as 1/sum(Wij) for a cartesian lattice
    if constexpr (DIM == 2)
    {
        std::vector<Point<2, real_number>> points;
        // place particles in a square from (-4dp,-4dp) to (4dp,4dp)
        for (int i = -4; i < 5; i++)
        {
            for (int j = -4; j < 5; j++)
            {
                Point<2, real_number> current_point = {i * argParameters.dp, j * argParameters.dp};
                points.push_back(current_point);
            }
        }

        // iterate vector and compute sum(Wij) for particle at (0,0)
        real_number sumWij = 0.0;
        Point<2, real_number> CentralParticle = {0.0, 0.0};
        for (auto point : points)
        {
            real_number r = getVectorNorm(point - CentralParticle);
            real_number Wij = Wab(r, argParameters.H, argParameters.Kquintic);
            sumWij += Wij;
        }

        argParameters.Vp = 1.0 / sumWij;
        // printf("Kernel volume: %.17g\n", argParameters.Vp);
        // printf("Dx*Dx: %.17g\n", argParameters.dp * argParameters.dp);
    }
    else if constexpr (DIM == 3)
    {
        std::vector<Point<3, real_number>> points;
        // place particles in a square from (-4dp,-4dp) to (4dp,4dp)
        for (int i = -4; i < 5; i++)
        {
            for (int j = -4; j < 5; j++)
            {
                for (int k = -4; k < 5; k++)
                {
                    Point<3, real_number> current_point = {i * argParameters.dp, j * argParameters.dp, k * argParameters.dp};
                    points.push_back(Point<3, real_number>(current_point));
                }
            }
        }

        // iterate vector and compute sum(Wij) for particle at (0,0)
        real_number sumWij = 0.0;
        Point<3, real_number> CentralParticle = {0.0, 0.0, 0.0};
        for (auto point : points)
        {
            real_number r = getVectorNorm(point - CentralParticle);
            real_number Wij = Wab(r, argParameters.H, argParameters.Kquintic);
            sumWij += Wij;
        }

        argParameters.Vp = 1.0 / sumWij;
        printf("Kernel volume: %.17g\n", argParameters.Vp);
        printf("Dx*Dx*Dx: %.17g\n", argParameters.dp * argParameters.dp * argParameters.dp);
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

    if constexpr (DIM == 3)
    {
        if (argParameters.bc[2] == PERIODIC)
            argParameters.Nboundary[2] = 0;
        else
            argParameters.Nboundary[2] = Nbound;
    }

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
        if constexpr (DIM == 3)
        {
            argParameters.length[2] = argParameters.dp * (argParameters.Nfluid[2] - argParameters.bc[2]);
        }
    }

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

    ComputeKernelVolume(argParameters);
    // argParameters.MassFluid = argParameters.rho0 * (DIM == 3 ? argParameters.dp * argParameters.dp * argParameters.dp : argParameters.dp * argParameters.dp);
    // argParameters.MassBound = argParameters.rho0 * (DIM == 3 ? argParameters.dp * argParameters.dp * argParameters.dp : argParameters.dp * argParameters.dp);

    argParameters.MassFluid = argParameters.rho0 * argParameters.Vp;
    argParameters.MassBound = argParameters.rho0 * argParameters.Vp;

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

    // cosmax
    argParameters.thetamax = 89.0 * M_PI / 180.0;
    argParameters.cosmax = cos(argParameters.thetamax);

    if (argParameters.SCENARIO != CUSTOM && argParameters.SCENARIO)
    {
        argParameters.ObstacleCenter[0] = (argParameters.length[0] + static_cast<real_number>(argParameters.bc[0]) * argParameters.dp) / 2.0;
        argParameters.ObstacleCenter[1] = (argParameters.length[1] + static_cast<real_number>(argParameters.bc[1]) * argParameters.dp) / 2.0;
        if constexpr (DIM == 3)
        {
            argParameters.ObstacleCenter[2] = (argParameters.length[2] + static_cast<real_number>(argParameters.bc[2]) * argParameters.dp) / 2.0;
        }
        // std::cout << "Obstacle center: " << argParameters.ObstacleCenter[0] << " " << argParameters.ObstacleCenter[1] << " " << argParameters.ObstacleCenter[2] << std::endl;
    }

    // add the number of processors to the filename
    const long int Nprc = v_cl.getProcessingUnits();
    argAuxParameters.filename += "_" + std::to_string(Nprc) + "prc";
}
