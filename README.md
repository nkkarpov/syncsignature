The code based 
To build the project, execute the following from the root directory (this will build all targets).
```bash
mkdir build
cd build
cmake ..
make
```
In the ``build`` directory you find the binary ``ted`` that executes the algorithms from command line.

The ted program has the following interface

```./ted <FILE> <ID> <K> <C> <R> <SIM> <SEED> [<NUMBER_OF_THREADS]```

<FILE> -- path to file
<ID> -- id of algorithm (0, 1, 2)
<K> -- distance threshold
<C> -- parameter for control size of parititions
<R> -- number of repetitions
<SIM> -- similarity threshold
<SEED> -- seed of randomness 
<NUMBER_OF_THREADS> -- number of threads for multi-thread version

