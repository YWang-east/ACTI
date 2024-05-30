<!-- omit in toc -->
# Adaptive Conservative Time Integration (ACTI)

This code use ACTI to solve the 2D compressible Navier-Stokes equation:

$$
    \frac{\partial \boldsymbol{U}}{\partial t} +
    \nabla \cdot (\boldsymbol{F}^c - \boldsymbol{F}^v)
    =\mathbf{0},
$$

where

$$
    \boldsymbol{U} = \left(\begin{array}{c}
    \rho \\
    \rho u \\
    \rho v \\
    \rho E
    \end{array}\right),
    \quad
    \boldsymbol{F}^c = \left(\begin{array}{cc}
    \rho u & \rho v\\
    \rho u u+p & \rho v u\\
    \rho u v & \rho v v+p \\
    u(\rho E+p) & v(\rho E+p)
    \end{array}\right),
$$

$$
    \boldsymbol{F}^v = \mu\left(\begin{array}{cc}
    0 & 0\\
    4/3 u_x - 2/3v_y & v_x + u_y\\
    u_y + v_x & 4/3 v_y - 2/3u_x \\
    u(4/3 u_x - 2/3v_y) + v(u_y + v_x) + \frac{\lambda}{\mu} T_x & 
    v(4/3 v_y - 2/3u_x) + u(v_x + u_y) + \frac{\lambda}{\mu} T_y
    \end{array}\right)
$$

<!-- omit in toc -->
## Table of Contents
- [Getting Started](#getting-started)
- [Architecture](#architecture)
- [Usage](#usage)
  - [Create your own case](#create-your-own-case)
  - [Add new input parameters](#add-new-input-parameters)

## Getting Started

After cloning the repository, make sure you have a compiler (e.g. `gcc`) and `CMake` installed. To compile code

```bash
cd $PATH_TO_PROJECT
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

Or alternatively, if you wish to debug the code

```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
```

In order to run, the code needs an input file named `control.txt` and an argument which specifies the path to the input file. You can copy one of the template input files located under `examples/`. We suggest doing the following

```bash
cd $PATH_TO_PROJECT
mkdir run
cp examples/control_mixing_layer.txt run/control.txt
./build/ACTI run
```

The result files will also be saved under `run`. To visualize the results, you can use the provided python script `plot.py`.

## Architecture
The codes contain three main parts
- `control`: handles reading input files
- `cell`: defines classes used in the simulation
- `driver`: contains all computing functions

The `Cell` class defined in [`cell.h`](src/cell.h) stores all solutions. Here are the key variables

`double _u[7]`: solution array
- `_u[0]`: density $\rho$
- `_u[1]`: momentum $\rho u$
- `_u[2]`: momentum $\rho v$
- `_u[3]`: total energy $\rho E$
- `_u[4]`: pressure $p$
- `_u[5]`: temperature $T$
- `_u[6]`: speed of sound $c$

`double _fl[4], _fr[4]`: fluxes on the left and right of the cell

`double _gl[4], _gr[4]`: fluxes on the top and bottom of the cell

## Usage

### Create your own case
We provide three test cases located under [`examples`](examples/). To create your own case, please follow the steps:

1. Create a control file (you can use one of the examples as a template) and assign a new case number, i.e. modify the `testCase` in the control file.
2. Find the function `void Driver::setupCase` in [`driver.cpp`](src/driver.cpp), define your initial conditions under the section
    ```cpp
    switch (testCase) {
        case $YOUR_CASE_NUMBER: {
            // define your initial condition
        }
    }
    ```
3. Modify the control inputs in `control.txt` as you like.

### Add new input parameters
You can add new input parameters on top of the existing ones in the examples. Just add a new line in `control.txt` with
```
YOUR_NEW_PARAMETER      VALUE   // some commments
```
In the code, it can be accessed using the `Control` class
```cpp
// instantiate a control object with name and path
Control ctr("control.txt", dir);  
// read from the control file  
ctr.read();
// declare your variable which takes value in the control file
double VAR_NAME = ctr.getDouble("YOUR_NEW_PARAMETER");
```
Similarly, you can also declare `int` and `bool` variables
```cpp
int  VAR_INT  = ctr.getInt("YOUR_INT_PARAMETER");
bool VAR_BOOL = ctr.getBool("YOUR_BOOL_PARAMETER");
```

