# SPH New Boundary Condition Simulation Code
![](Header.png)

This project is a simulation code for Smoothed Particle Hydrodynamics (SPH) used to study the new boundary condition developed in the master thesis work "A single layer ghost particle boundary method for Smoothed Particle Hydrodynamics". The code is written in C++ and uses the OpenFMP library.

## Installation
First follow the instructions from the [OpenFPM library](https://openfpm.mpi-cbg.de/building/) to install it. If you wish to enable the CUDA support, you will need to set GPU_CUDA_SUPPORT=1 during the OpenFPM installation. 

Run the openFPM tests to make sure the library is correctly installed. Make sure you are able to run the example inside example/Vector/7_SPH_dlb and  example/Vector/7_SPH_dlb_gpu

Clone this repository inside the example/Vector/ directory of the OpenFPM library. 

The GPU version has been tested for CUDA toolkit release 12.0 in Ubuntu 22.04.5.

## Usage
Compile the code using make.

Run the code using following command:

```bash
./sph input_file.xml
```

In XML/example/ you can find example input files that can be used to run the code.

You can modify the input parameters in the XML file to change the simulation parameters.

To create whole new simulation scenarios, modify the code inside /src/CreateParticleGeometry.hpp and /src/CreateParticleGeometry.cu. If you want to read additional parameters, modify /src/InitializeParameters.hpp and /src/InitializeParameters.cu.

