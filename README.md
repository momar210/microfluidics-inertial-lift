# microfluidics-inertial-lift
 Tools for simulating particle migration, focusing, and separation in inertial microfluidic devices with ANSYS Fluent

### General Strategy for Inertial Microfluidic Simulation
Inertial microfluidic devices often have several repeated segments to modify fluid flow and particle positions within a channel. Meshing and simulating the entire device would be computationally expensive and applying periodic boundary conditions assumes the device is infinitely long.

To simulate particle flow in an inertial microfluidic device with a finite number of repeated segments, we first release particles at the inlet from an initial DPM (Discrete Phase Model) injection definition. Particles are then sampled at the outlet plane and repeatededly reinjected in subsequent segments.

### Calculating Inertial Lift Forces
Calculation of the inertial lift forces is performed using a DEFINE_DPM_BODY_FORCE UDF (User Defined Function). The lift forces are calculated based on the position and diameter of the particle and the cross sectional area and Reynold's number of the channel. The function that relates these parameters to the lift force was calculated using the data and MATLAB neural network trainer from Su et al. (2021) [1]. The neural network was trained in MATLAB as described in the paper and the resulting shallow neural network was deployed with `genFunction` and ported to C with `codegen` for use in the UDF.

### Currently running with:
- Windows 10 Enterprise (64 bit, OS Build 19041.1415)
- ANSYS Academic Student 2021 R2
- Microsoft Visual Studio 2017

### Current workflow:
1. Generate geometry of one repeated segment of a periodically repeating inertial microfluidic device
2. Generate mesh
3. Export mesh as Fluent mesh into `".\microfluidic-simulation_files\user_files\{filename}.msh"`
4. Open Fluent with case `".\microfluidic-simulation_files\user_files\microfluidic-simulation.cas.h5"`
5. Import mesh
6. Set flow rates and channel parameters
7. Initialize and Calculate
8. Simulate particle release and flow through repeated segments by Reading Scheme (File > Read > Scheme), `".\microfluidic-simulation_files\user_files\particle-tracking-repeated-steps-with-bootstrap.scm"`
9. Display results

### Future updates:
- Include DPM Body Force UDF for calculation of inertial lift
- Include UDF for precomputing inertial lift at each cell centroid for inspection of inertial force field
- System for generative design of channel cross sections
- System for description of channel cross sections with functions
- System for optimization of channel geometry given a desired particle distribution


## References
> 1. Su, J., Chen, X., Zhu, Y. & Hu, G. Machine learning assisted fast prediction of inertial lift in microchannels. Lab Chip 21, 2544–2556 (2021). [doi:10.1039/D1LC00225B](https://doi.org/10.1039/D1LC00225B)