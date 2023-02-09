Fast and accurate protein substructure searching with simulated annealing and GPUs
# Fast and accurate protein substructure searching with simulated annealing and GPUs

## Software

Imported from https://sites.google.com/site/alexdstivala/home/satabsearch

This software can be used freely for any purpose, modified, redistributed, etc.
with no restrictions, except where noted that NVIDIA SDK sample code
is used by "This software contains source code provided by NVIDIA Corporation.",
namely, for the Mersenne Twister pseudorandom number generator. In this
case the only restriction is that the preceding NVIDIA statement is retained
(as it is here and in comments in the relevant code fragments), as specified
in the NVIDIA CUDA SDK End User License Agreement.

However we would appreciate it if you acknowledge
your use of it, and in particular if you would cite our paper
in any publication that makes use of it.

#### Current implementation

(New, as of August 2021).

The code in nvcc_src_current/ has been updated to work with more recent
versions of CUDA and newer NVIDIA GPUs. Specifically, it has been
built and tested with the NVIDIA nvcc compiler and CUDA SDK 
release 11.1, V11.1.74 on a 64-bit Linux system and a
NVIDIA GeForce GTX 1080 (Founders Edition 8GB GDDR5X 2560 CUDA cores) GPU.


#### CUDA 5.0 implementation

(New, as of November 2013).

SA Tableau Search was originally developed on CUDA 2.3, which did not
even have a device pseudorandom number generator library, so the Mersenne
Twister example from the SDK was adapted to provide a RNG on the device.
Later versions provide the [CURAND library](http://docs.nvidia.com/cuda/curand/), and CUDA 5.0 changed some APIs,
so the code has been updated to run on more modern GPUs with CUDA 5.0.
This source code (and compiled executable) was built and tested
on a 64-bit Linux system with CUDA 5.0.35 and a Tesla M2070 GPU.

### Web server

SA Tableau Search is also available for online use via the [Pro-origami web server](http://munk.cis.unimelb.edu.au/pro-origami).
However, the server is does not have a GPU, so
SA Tableau Search runs on the server CPU and is therefore orders of magnitude slower
than running on a modern GPU.

## Tableaux database

The tableaux database was built with 
`buildtableauxdb.py -t dssp -3 -5 -p none` 
The distance matrices were built with 
`buildtableauxdb.py -t dssp -3 -5 -p none -d` and the combined
tableaux and distance matrix database in ASCII format created with
the convdb2.py script, and sorted by size with the 
sort_tableux_db_pbs_script.sh script.
These scripts are included in the scripts directory.

## Usage

Usage: `cudaSaTabsearch [-c] [-q dbfile] [-r restarts] < inputfile`

`-c `: run on host CPU not GPU card

`-q` : query list mode: instead of reading query data on stdin
     just as in the original Fortran version tlocsd (https://github.com/stivalaa/qptabsearch), a list
     of query sids to be read from the database is read on stdin (one per
     line),
     and db filenaame is specified on command
     line. In this mode options are assumed as LORDER=T, LTYPE=T,
     LSOLN=N. The output is still to stdout, but each query following
     immediately from the previous (can parse using the  header comment
     information lines as separators).

`-r restarts`: number of restarts (iterations of cooling schedule).
             Should be a multiple of blocksize. Defaults to 128.

The 'database' to search is an ASCII file of tableaux
(Omega matrices) in format described in rdtabd.f (https://github.com/stivalaa/qptabsearch).

The results are printed to stdout as 

`name rawscore norm2score z-score p-value`

Both the name of the database file to read, and the actual
query tableau are read from stdin. 
The first line is the name
of the database file.
The second line is for options. There are currently 3 logical
options, for SSE type constraint (only allow SSEs of same type to
match) and ordering constraint (disallow out of sequence order 
matches). The third is to output not just the scores but also solution
vector values.
They are single character logical values (T or F).
First is type, second is order, third is solution output,
separated by one space.

The subsequent lines are multiple query structures (tableau and
distance matrix for each), separated by a blank line. I.e. the
same format as the database.

The tableau is in the same format as
each tableau entry in the database i.e.:

The first line of an entry is the identifier and
The first line of an entry is the identifier and
order of tableau (i.e. dimension of square array), then
each subsequent row is a row of the tableau, lower triangle
only (since it is symmetric).
The diagonal entries are meaningless (self-angle) in tableaux,
and are included instead to specify the SSE type, with
the following codes:


```
e     beta strand
xa    alpha helix
xi    pi helix
xg    3_10 helix
```

Width of identifier is 8 chars, blank padded on right,
width of order is 4 digits, blank padded on left.
There is a single space between identifier and order.
Each entry in tableau is two characters, with a space betwen
each on a line, and one line
per row of matrix.

Following the tableau is the distance matrix.
Each row is a row of the distance matrix, lower triangle
only (since it is symmetric).
The diagonal entries are meaningless (self-distance)
and are included instead to specify the SSE type, with
the following codes:

```
0.000 beta strand
1.000 alpha helix
2.000 pi helix
3.000 3_10 helix
```

Each entry in matrix is in Angstroms format
F6.3 with a space between each on a line, and one line
per row of matrix.

E.g.:

```
/local/charikar/astivala/tableauxdb/astral/tableauxdistmatrixdb.ascii
T T F
D1UBIA_    8
e  
OT e  
LE RT xa 
PD OS RD xg 
RT LE RT LS e  
LE RD LE LS OT e  
RT LS LS RD PE OS xg 
PE RT LE RD OT PE RT e  
 0.000 
 4.501  0.000 
 1.662 10.386  1.000 
16.932 17.644  9.779  3.000 
10.588 13.738 11.815 10.527  0.000 
15.025 18.692 17.143 15.341  6.466  0.000 
15.298 17.276 16.276 20.075 13.264 11.610  3.000 
 7.549 11.072 12.248 12.446  4.583  9.903 15.689  0.000

D1AE6H1   13
e
PD e
OT OS e
LS LS RD xg
LE LE RT RT e
RT RT LE LE OT e
RT RT LE LE OT PE e
LE LE RT RT PE OT PE e
RT OT PE LS OT LE PE OT e
PE PE OT LS PE RT OT PE OT e
OT OT PE RD OS LE PE OT PE OT xg
OT RT PE RD OT PE PE OT PE OT PD e
PD PE OT LS LE RT RT LE OT PE OT OT e
 0.000
19.130  0.000
 8.850 13.371  0.000
14.608 29.221 15.945  3.000
12.469 19.135 11.231 16.008  0.000
18.479 21.128 16.959 21.982  6.730  0.000
16.153 22.704 13.140 13.210  6.909 10.709  0.000
20.850 24.610 16.558 16.527 10.946 12.552  4.935  0.000
15.604 18.394  8.791 14.366 11.402 16.188  8.316  9.609  0.000
13.949 13.565  5.751 17.301 10.771 15.314 10.725 12.661  4.876  0.000
24.234 12.620 19.140 31.786 17.166 14.790 20.733 21.202 20.224 16.665  3.000
 9.731 17.355  9.936 16.942  3.841  8.797 10.327 14.566 13.021 11.226 17.023 0.000
16.856  5.985 12.706 27.454 15.011 16.146 19.829 22.156 17.744 13.079  9.541 12.996  0.000
 ```

See `*.input` and `*.sh` files for  more examples, and [README_example_usage.txt](README_example_usage.txt) for worked example of building database and query, and running the query against the database.

## Reference

If you use our software, data, or results in your research, please cite:

Stivala, A., Stuckey, P., Wirth, A.,
[Fast and accurate protein substructure searching with simulated annealing and GPUs](http://www.biomedcentral.com/1471-2105/11/446) BMC Bioinformatics 2010, 11:446

For details on how the z-scores and p-values are computed, see http://munk.cis.unimelb.edu.au/pro-origami/tableaustats.pdf, or Chapter 6 of:

Stivala, A. D. (2010). [Algorithms for the study of RNA and protein structure](http://hdl.handle.net/11343/35772). PhD thesis, Computer Science and Software Engineering, The University of Melbourne.
