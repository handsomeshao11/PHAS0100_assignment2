PHAS0100Assignment2
------------------

[![Build Status](https://travis-ci.com/[USERNAME]/PHAS0100Assignment2.svg?branch=master)](https://travis-ci.com/[USERNAME]/PHAS0100Assignment2)
[![Build Status](https://ci.appveyor.com/api/projects/status/[APPVEYOR_ID]/branch/master)](https://ci.appveyor.com/project/[USERNAME]/PHAS0100Assignment2)


Purpose
-------

This project serves as a starting point for the PHAS0100 Assignment 2 Gravitational N-body Simulation coursework. It has a reasonable folder structure for [CMake](https://cmake.org/) based projects,
that use [CTest](https://cmake.org/) to run unit tests via [Catch](https://github.com/catchorg/Catch2). 

Further information on the specific project is left as an exercise for the student.


Credits
-------

This project is maintained by [Dr. Jim Dobson](https://www.ucl.ac.uk/physics-astronomy/people/dr-jim-dobson). It is based on [CMakeCatch2](https://github.com/UCL/CMakeCatch2.git) that was originally developed as a teaching aid for UCL's ["Research Computing with C++"](http://rits.github-pages.ucl.ac.uk/research-computing-with-cpp/)
course developed by [Dr. James Hetherington](http://www.ucl.ac.uk/research-it-services/people/james)
and [Dr. Matt Clarkson](https://iris.ucl.ac.uk/iris/browse/profile?upi=MJCLA42).


Build Instructions
------------------
The function of files:
Code/lib/nbsimBasicTypes.cpp and Code/lib/nbsimBasicTypes.h  -- The class and some functions of initialize a body.
Code/CommandLineApps/solarSystemSimulator.cpp -- Run the solar system body with a timestep and a number of year.
Code/CommandLineApps/nbsimMyFirstApp.cpp -- Run 2000 bodies by initialized randomly with a timestep, a number of year and a number of machine processors.
TestingnbsimBasicTest.cpp -- Check the nbsimBasicTypes.cpp and nbsimBasicTypes.h can be used correctly.

Open the PHAS0100Assignment2 folder.

```Bash
mkdir build
cd build
cmake ..
make
```

Run Instructions
----------------
To run the solarSystemSimulator.cpp:

Check you are in the build folder now.

```Bash
./bin/solarSystemSimulator -help
```

Then watch the print content in the terminal, and input the command as help message.
Example:
```Bash
./bin/solarSystemSimulator 0.0000274 1.0
```

To run the nbsimMyFirstApp.cpp:

Check you are in the build folder now.

```Bash
./bin/nbsimMyFirstApp -help
```

Then watch the print content in the terminal, and input the command as help.
Example:
```Bash
./bin/nbsimMyFirstApp 0.000274 1 8
```

To run the nbsimBasicTest.cpp:

Check you are in the build folder now.

```Bash
./bin/nbsimBasicTest ctest -V
```

Then watch the print content in the terminal, and you can see how many test has passed.


Energy and Benchmark result:
---------------------------
System initial energy are : kinetic energy is 0.187358, potential energy is -0.355517, total_energy is -0.168158

There are cases with different time steps and 100 years:

case 0: timestep 0.0000274  num_time 100

System final energy are : kinetic energy is 0.135274, potential energy is -0.303163, total_energy is -0.16789
CPU time used: 335077.56 ms
Wall clock time has passed: 335077.56 ms

case 1: timestep 0.0001 num_time 100

System final energy are : kinetic energy is 0.149917, potential energy is -0.316279, total_energy is -0.166362
CPU time used: 478510.46 ms
Wall clock time has passed: 478510.46 ms

case 2: timestep 0.001 num_time 100

System final energy are : kinetic energy is 0.182478, potential energy is -0.324498, total_energy is -0.14202
CPU time used: 507169.78 ms
Wall clock time has passed: 507169.78 ms

case 3: timestep 0.00274 num_time 100

System final energy are : kinetic energy is 0.0876923, potential energy is -0.182417, total_energy is -0.0947247
CPU time used: 398754.43 ms
Wall clock time has passed: 398754.43 ms

case 4: timestep 0.0274 num_time 100

System final energy are : kinetic energy is 0.0439279, potential energy is -0.0566199, total_energy is -0.012692
CPU time used: 354186.91 ms
Wall clock time has passed: 354186.91 ms

From the result, the case 0 with timestep 0.0000274, number of year 100 has the best accuracy and a great time. 
More result you can use the solarSystemSimulator.cpp project.


Energy and Benchmark result with the parallel computing:
-------------------------------------------------------
Non-parallel: timestep 0.000274  num_time 1 number_proccesors 8

System initial energy are : kinetic energy is 81951.3, potential energy is -1.04652e+09, total_energy is -1.04643e+09
System final energy are : kinetic energy is 1.18048e+20, potential energy is -1.54182e+08, total_energy is 1.18048e+20

CPU time used: 105826.13 ms
Wall clock time has passed: 105828.79 ms
Processor number: 8

parallel: timestep 0.000274  num_time 1 number_proccesors 8

System initial energy are : kinetic energy is 81951.3, potential energy is -1.04652e+09, total_energy is -1.04643e+09
System final energy are : kinetic energy is 1.18048e+20, potential energy is -1.54182e+08, total_energy is 1.18048e+20

CPU time used: 186243.53 ms
Wall clock time has passed: 25038.43 ms
Processor number: 8

From the result, we can see there a significant improvement in program computing.
