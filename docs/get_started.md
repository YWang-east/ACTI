# Getting Started

After cloning the repository, make sure you have a compiler (e.g. `gcc`) and `CMake` installed. To compile code

```bash
$cd $PATH_TO_PROJECT
$mkdir build
$cd build
$cmake -DCMAKE_BUILD_TYPE=Release ..
$make
```

Or alternatively, if you wish to debug the code

```bash
$cmake -DCMAKE_BUILD_TYPE=Debug ..
$make
```

In order to run, the code needs an input file named `control.txt` and an argument which specifies the path to the input file. You can copy one of the template input files located under `examples/`. We suggest doing the following

```bash
$cd $PATH_TO_PROJECT
$mkdir run
$cp examples/control_mixing_layer.txt run/control.txt
$./build/ACTI run
```

The result files will also be saved under `run`. To visualize the results, you can use the provided python script `plot.py`.

   