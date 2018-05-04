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
buildtableauxdb.py -t dssp -3 -5 -p none 
The distance matrices were built with 
buildtableauxdb.py -t dssp -3 -5 -p none -d and the combined
tableaux and distance matrix database in ASCII format created with
the convdb2.py script, and sorted by size with the 
sort_tableux_db_pbs_script.sh script.
These scripts are included in the [scripts.tar.gz](scripts.tar.gz)
file.

## Reference

If you use our software, data, or results in your research, please cite:

Stivala, A., Stuckey, P., Wirth, A.,
[Fast and accurate protein substructure searching with simulated annealing and GPUs](http://www.biomedcentral.com/1471-2105/11/446) BMC Bioinformatics 2010, 11:446

