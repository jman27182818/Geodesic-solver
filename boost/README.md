

# Compilation

Make -jx

where x is the number of processors you want to use to compile.


# Running:

```
g++ -I /usr/local/boost/ geodesicSolver.cpp
./geodesicSolver -metric <path to metric file> -
```

# Dependencies

On ubuntu run:

```
sudo apt-get install libboost-all-dev
```