# Angular Multi-Index Hashing (AMIH)
This is an implementation of the Angular Mutli-Index Hashing which performs exact angular *K* nearest neighbors in the binary space. It is built on top of multi-index hashign proposed in [1] [(mih)](https://github.com/norouzi/mih)
## Build Instructions
The requirments of building the project are, *cmake, make, hdf5 lib* and *hdf5-dev*. To build it, run:

```
cmake CMakeLists.txt
make
```
if successful, two files are created: *amih* which executes amih technique and *linscan* for executing linear scan.

## Dataset

A sample data of 1M sift binary codes is given in the *data* folder. It contains two matrices, *B* which contains the binary codes of the base dataset and *Q* which contains the query set. In each matrix, each 8 bits are stored an unsigned integer.

## Running AMIH
The sysntax of running AMIH is as follows:
```
./amih <input file> <output file> -N #points -B #bits -Q #queries  -K #NN

```
where *NN* is the number of nearest neighbors to retrive.

Example:
```
./amih lsh_64_1M.mat out.h5 -N 1000000 -B 64 -Q 100 -K 1

```

The syntax of running linear scan is as follows:

```
./linscan <input file> <output file> -N #points -B #bits -Q #queries -K #NN

```
Example:

```
./linscan <input file> <output file> -N 1000000 -B 64 -Q 100 -K 1
```

[1] Fast Exact Search in Hamming Space with Multi-Index Hashing, M. Norouzi, A. Punjani, D. J. Fleet, IEEE TPAMI 2014
