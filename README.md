# SPIM Single Cell Analysis

### Single cell analysis of multicolor SPIM Data. The data is assumed to be in HDF5 form and has gone through [initial aligning and processing](https://github.com/alonyan/bigstitchparallel)

This is code implemented in MATLAB to translate a multicolor HDF5 SPIM dataset into a table of cell positions and attributes. The cells are than tracked using the algorithm described in [Jaqaman et al. 2008](https://www.nature.com/articles/nmeth.1237)
