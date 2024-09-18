# SHARK
<u>**S**tructural **H**ealth **A**ssessment and **R**eal-time **K**inematics</u> observer (SHARK) is a real-time state estimation platform for fixed-bottom and floating wind turbines. Using a [semi-linear dynamics model](https://www.overleaf.com/read/vwhkfsmnwypk#2b01ba), SHARK computes a real-time state estimate using either an Extended Kalman filter or Square-Root Unscented Kalman filter (CITE).

# Preparation Steps

To ensure SHARK can work correctly, do the following:

1. Ensure that the entire SHARK directory is added to the Matlab path.

# Running TurbSim

To run Turbsim and generate a turbulent wind file, first prepare the corresponding folder under OpenFAST/Turbsim. The folder name must match the input file name. In this folder, prepare the Turbsim input file for Turbsim v2.0. An example is provided for the Gulf of Maine Monhegan test site [CITE]. 

To run Turbsim and generate a turbulent wind file, do the following:

1. Create a folder under OpenFAST/Turbsim
2. Create a Turbsim input file in this folder with the same name.
3. Navigate to the home SHARK directory and open the "TurbSimDriver.m" script.
4. Adjust the variable "home_dir" as needed. This should be the current SHARK directory on your computer.
5. Replace the "input_case" value with the name of the folder/file you just created and run.

# Running OpenFAST

To generate OpenFAST results for comparison, do the following:

1. Open the "OpenFAST_Driver.m" script and adjust the "home_dir" variable to the current SHARK directory location.
2. Change "sim_folder" to the desired simulation output name. You do not need to create this folder manually.
3. Change "model" to the name of the OpenFAST model folder.
4. Adjust the OpenFAST files in the corresponding model folder to the desired simulation. Minimum outputs required for any simulation are TTDspFA, RootMyb1-3, BldPitch1-3, Wind1VelX, RotSpeed, Azimuth.
5. Models build in other versions of OpenFAST are supported. Simply adjust the "version" and "bin_name" variables and place the compiled OpenFAST executables in the corresponding locations.
6. Run the OpenFAST driver. Select "yeah, sure" to create the new directory, or "OMG, NO" if there is an error in the output folder name.

# Fixed (Land-Based) Turbine Simulation

To run the fixed wind turbine case, simply run the "FixedTurbine_Driver.m" file. Within this driver file you can adjust the simulation time by adjusting the time vector max value. Note that dt **must** match the OpenFAST output dt.

## Adjusting the Observer



