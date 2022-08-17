The code for the SyncSignature algorithm. 

To build the project, execute the following from the root directory.
```bash
mkdir build
cd build
cmake ..
make ted
```
In the ``build`` directory you can find the binary file ``ted`` that executes the algorithms from command line.

The ted program has the following interface

```./ted <FILE> <ID> <K> <C> <R> <SIM> <SEED> <NUMBER_OF_THREADS> <CUT>```

``FILE`` - path to input file. File should contain n lines. Each line represents a tree in the format of ``{rootLabel{child1subtree}{child2subtree}...{childksubtree}}``. For example, the string ``{a{b{c}{d}{e}}{f{g}}}`` corresponds to tree 
~~~
                   a 
                  / \ 
                 b   f   
                /|\   \    
               c d e   g 
~~~


``ID`` - id of the algorithm. It can be set to 0 (``EJoin``), 1 (``BJoin``), or 2 (``TJoin``).

``K`` - distance threshold.

``C`` - neighbourhood resolution (we set it to be ``0.3`` by default).

``R`` - number of repetitions.

``SIM`` - similarity threshold (we set it to be ``K/5`` by default).

``SEED`` - seed of randomness.

``NUMBER_OF_THREADS`` - number of threads for the multi-thread version.

``CUT`` - minimum size of a tree (we set it to be ``1000`` by default).

Example of usage (after unziping datasets):
```./ted ../datasets/python_sorted.bracket 0 10 0.3 11 2 239 1 1000```

Datasets : https://indiana-my.sharepoint.com/:f:/g/personal/nkarpov_iu_edu/Ek7a2RkL6CNBnPvVn16SsYYB12IHnvjMzPlQU4FJlYZuhQ?e=t3YH7a
