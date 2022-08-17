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

``FILE`` - path to input file. File should contain n lines, each line represents a tree in format ``{rootLabel{child1subtree}{child2subree}...{childksubtree}}``. For example, the string ``{a{b{c}{d}{e}}{f{g}}}`` corresponds to tree 

``
    a

   / \

  b   f

 /|\   \

c d e   g
``


``ID`` - id of algorithm, it can be set to 0 (``EJoin``), 1 (``BJoin``), or 2 (``TJoin``).

``K`` - distance threshold.

``C`` - neighbourhood resolution (equal to ``0.3`` for most experiments).

``R`` - number of repetitions.

``SIM`` - similarity threshold (equal to ``K/5``)

``SEED`` - seed of randomness 

``NUMBER_OF_THREADS`` - number of threads for multi-thread version

``CUT`` - lower bound on the size of tree (1000 for most experiments)

Example of usage (after unziping datasets):
```./ted ../datasets/python_sorted.bracket 0 10 0.3 11 2 239 1 1000```

Datasets : https://indiana-my.sharepoint.com/:f:/g/personal/nkarpov_iu_edu/Ek7a2RkL6CNBnPvVn16SsYYB12IHnvjMzPlQU4FJlYZuhQ?e=t3YH7a
