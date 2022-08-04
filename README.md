# SFC
Space Filling Curve sparse matrix reordering implementations
Stephan Frickenhaus(1,*), Annika Fuchs, Natalja Rakowsky(1) (2022)

(1) Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research, Bremerhaven, Germany

(*) Contact: stephan.frickenhaus@awi.de

This is a collection of codes for Space Filling Curve renumbering of sparse matrices.
Main code of SFC-sonumbering is based on a work by Rakowsky and Fuchs (2011).

There is a C++-interface to R, and a demontrator of improved cache-efficiency for Matrix Vector product (SparseM, CSR-fornmat) with a pseudo FE-matrix created on the fly from a FE-like 2D-mesh. 
I also include an R-script applying SFC on the matrix structure, depicted in a 2D-graph representation, rather than on a FE-like mesh. This is achieved by a very fast graph-layout in 2D by pivotMDS.

To demonstrate cache effects in a multi-core approach, e.g., where multiple threads share a common cache, an R script is included for parallel matrix-vector-products on smaller pseudo-FE-matrices. 


Findings:
For large N (1.5e6) a clearly improved cacche-efficiency on sparse matrixvector product is demonstrated for both, the mesh-reordered as well as the matrix-layout reordered matrix.
We propose to compute statistics of standard deviations of column indices across all rows to show the effect of cache-reuse in accessing nearby memory positions in indirect addressing, which is common to sparse matrix-vector products.


Outlook:
Other sparse matrix formats, like ELLPACK-R, SELL-P, SELL-C-sigma should be implemented and benchmarked as well, e.g. to show cache-effects of SFC-sorting  on GPUs under parallel processing (openACC, native CUDA).

Codes:

resortgrid_SFC.c : original code from Rakowsky & Fuchs (2011);
resortgrid_SFC_Rcpp.cpp : R-interface to the SFC function from Rakowsky & Fuchs (2011);
benchmarks.R : R-script to benchmark sparse mat-vec-product in CSR-matrix-format, NO column-indices sorted;
benchmarks_column_sort.R :  R-script to benchmark sparse mat-vec-product in CSR-matrix-format, WITH column-indices sorted;
benchmarks_column_sort_parallel.R :  R-script to benchmark sparse mat-vec-product in CSR-matrix-format on muultiple threads, WITH column-indices sorted

Standard config of the pseudo FE matrix:
N=1.5e6 nodes (dimension of vector)
n=15 non-zeros per row

Order of execution:
1 bencharks.R  (also saves the original matrix in rdat-file)
2 bencharks_colmn_sort.R (also saves the original matrix in rdat-file)

3 resort_matrix_graph.R (uses saved original matrix from benchmarks.R)
4 resort_matrix_graph_column_sort.R (uses saved original matrix from benchmarks_colmn_sort.R)


References

Rakowsky N. & Fuchs A. Efficient local resorting techniques with space filling curves applied
to the tsunami simulation model TsunAWI. In IMUM 2011 - The 10th International Workshop
on Multiscale (Un-)structured Mesh Numerical Modelling for coastal, shelf and global ocean
dynamics, August 2011. http://hdl.handle.net/10013/epic.39576.d001

used R Packages:

RANN             for fast nearest neighbor search on a random 2D mesh;
Rcpp             for the R-interface to sfc-code;
SparseM          for CSR matrix-vector product etc.; 
Matrix           for interconversion to adjacency graph format;
igraph           for adjacency graph;
graphlayouts     for pivotMDS fast graph layout giving 2D coordinates for SFC without Mesh;
microbenchmark   for bencharking runtimes;
parallel         for parallel tests on multi-core / CPUs with shared cache;

Tested with R.version :
platform       x86_64-w64-mingw32          
arch           x86_64                      
os             mingw32                     
system         x86_64, mingw32             
status                                     
major          4                           
minor          0.5                         
year           2021                        
month          03                          
day            31                          
svn rev        80133                       
language       R                           
version.string R version 4.0.5 (2021-03-31)
nickname       Shake and Throw   
