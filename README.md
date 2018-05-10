# pc-exp

This is the repo that stores the coding works of Parallel Computing Course at 2017-2018 Spring Semester

## System Requirement

You should install CMake with the version at least 3.9

You should install a reasonable C++11-supported and openmp-supported compiler (I use gcc-7)

## Config

## Installation (For Unix-Like OS)

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

However, since the fxxking cmake can only use one compiler at one time, so some code should be compiled manually (mainly for mpi code). I place a script for that in the src/

## Usage

Just execute them.

## Parts finished

## Todo List

No idea
