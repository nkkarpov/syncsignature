# SyncSignature

The code for the SyncSignature algorithm.

### Build Instructions

To build the project, execute the following from the root directory.

```bash
mkdir build
cd build
cmake ..
make ted
```

In the ``build`` directory you can find the binary file ``ted`` that executes the algorithms from command line.

### Run Instructions

The ted program has the following interface

```./ted <FILE> <ID> <K> <C> <R> <SIM> <SEED> <NUMBER_OF_THREADS> <CUT>```

``FILE`` - path to input file. File should contain n lines. Each line represents a tree in the format
of ``{rootLabel{child1subtree}{child2subtree}...{childksubtree}}``. For example, the string ``{a{b{c}{d}{e}}{f{g}}}``
corresponds to tree

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

### License

All code in this repository is currently under the [MIT licence](https://opensource.org/licenses/MIT).

### Datasets

We use 6 datasets for experiments in the main paper.

- **Swiss** : protein sequence data ([raw data](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz))
- **Python** : abstract syntax trees of Python files ([raw data](http://files.srl.inf.ethz.ch/data/py150.tar.gz))
- **JScript** : abstract syntax trees of JavaScript files ([raw data](https://files.sri.inf.ethz.ch/data/js_dataset.tar.gz))

The following three datasets obtained from the previous after filtering out trees of size less than 1 ,000.

- **Swiss1K**
- **Python1K**
- **JScript1K**

#### Statistics on Datasets

| name      | # trees | min. size | max. size | avg. size | 
|-----------|---------|-----------|-----------|-----------|
| Swiss     | 565,254 | 105       | 48,286    | 917       |
| Python    | 148,270 | 1         | 46,481    | 948       |
| JScript   | 142,373 | 4         | 1,716,813 | 2,619     |
| Swiss1K   | 122,772 | 1,000     | 48,286    | 1,902     |
| Python1K  | 35,754  | 1,000     | 46,481    | 3,016     |
| JScript1K | 39,110  | 1,000     | 1,716,813 | 9,006     |

#### Link to Datasets

https://indiana-my.sharepoint.com/:f:/g/personal/nkarpov_iu_edu/Ek7a2RkL6CNBnPvVn16SsYYB12IHnvjMzPlQU4FJlYZuhQ?e=t3YH7a
