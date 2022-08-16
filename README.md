The code of SyncSignature algorithm. 

To build the project, execute the following from the root directory.
```bash
mkdir build
cd build
cmake ..
make ted
```
In the ``build`` directory you find the binary ``ted`` that executes the algorithms from command line.

The ted program has the following interface

```./ted <FILE> <ID> <K> <C> <R> <SIM> <SEED> <NUMBER_OF_THREADS <CUT>```

``FILE`` - path to input file. File should contain trees in format ``{a{b}{c}}`` (see datasets for specific examples).

``ID`` - id of algorithm, it can be set to 0 (``EJoin``), 1 (``BJoin``), or 2 (``TJoin``).

``K`` - distance threshold.

``C`` - neighbourhood resolution (equal to ``0.3`` for most of experiments).

``R`` - number of repetitions.

``SIM`` - similarity threshold (equal to ``K/5``)

``SEED`` - seed of randomness 

``NUMBER_OF_THREADS`` - number of threads for multi-thread version

``CUT`` - lower bound on the size of tree (1000 for most of experiments)

Example of usage (after unziping datasets):
```./ted ../datasets/python_sorted.bracket 0 10 0.3 11 2 239 1 1000```
