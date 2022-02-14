# microfluidics-inertial-lift
 Tools for simulating particle migration, focusing, and separation in inertial microfluidic devices with ANSYS Fluent

### Currently running on:
Windows 10 Enterprise (OS Build 19041.1415)
ANSYS Academic Student 2021 R2
Visual Studio 2017

### Current workflow:
1. Generate geometry of one repeated segment of a periodically repeating inertial microfluidic device
2. Generate mesh
3. Export mesh as Fluent mesh into `".\microfluidic-simulation_files\user_files\{filename}.msh"`
4. Open Fluent with case `".\microfluidic-simulation_files\user_files\microfluidic-simulation.cas.h5"`
5. Import mesh
6. Set flow rates and channel parameters
7. Initialize and Calculate
8. Simulate particle release and flow through repeated segments by Reading Scheme, `".\microfluidic-simulation_files\user_files\particle-tracking-repeated-steps-with-bootstrap.scm"`
9. Display results