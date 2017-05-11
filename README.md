# Angular Multi-Index Hashing (AMIH)
This is an implementation of the Angular Mutli-Index Hashing which performs exact angular *K* nearest neighbors in the binary space. It is built on top of multi-index hashign proposed in [1] [mih](https://www.google.com)
## Compile
The requirments of building the project are, *cmake, make, hdf5 lib* and *hdf5-dev*. To build it, run:

```
cmake CMakeLists.txt
make
```
if successful, two files are created: *amih* which executes amih technique and *linscan* for executing linear scan 
