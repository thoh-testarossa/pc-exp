# pc-exp

This is the repo that stores the coding works of Parallel Computing Course at 2017-2018 Spring Semester

## System Requirement

You should install CMake with the version at least 3.9

You should install a reasonable C++11-supported and openmp-supported compiler (I use gcc-7)

## Config

Maybe you should check the compiler in the CMakeLists.txt and the script placed in src/ which I will mention behind to make sure these config files can works well.

## Installation (For Unix-Like OS)

### For OpenMP programs

1. Clone this project to your local machine and move to the root directory of this project

```
cd [path]/[to]/pc-exp
```

2. Create a build directory

```
mkdir build

cd build
```

3. Execute cmake

```
cmake ..
```

4. Compile this project

```
make
```

After that, some executables will appear in the directory "build".

### For MPI programs

Since the fxxking cmake can only use one compiler at one time and MPI codes can only be compiled by `mpicc` or `mpicxx`, so these codes should be compiled manually. I place a script for that in the src/ for easier compilation

## Usage

### For OpenMP programs

It's okay to just execute them without doing nothing. By default the OpenMP programs will use only 1 thread as the config to run. To change the concurrent thread number, add a number behind the program name as the parameter which determines thread number.

### For MPI programs

To run MPI programs, use instructions in this form:

```
mpirun -np [num-of-processors] [name-of-program]
```

which -np options determines the processors you want to use in an MPI model.

## Parts finished

## Todo List

No idea
