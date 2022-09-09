Example query building and running (on small test database in github):

    [stivala@icslogin01 spike]$ echo "/home/stivala/public_github/cuda_satabsearch/old/nvcc_src_cuda5/tableauxdistmatrixdb.small.ascii" > ubiquiin.query
    [stivala@icslogin01 spike]$ echo "T T F" >> ubiquiin.query 
    [stivala@icslogin01 spike]$ python2 ~/proorigami/pro-origami/pytableaucreate.py -3 -5 -t stride -b -f  ~/pconpy/tests/pdb_files/1ubq.pdb >> ubiquiin.query 
    [stivala@icslogin01 spike]$ 
    [stivala@icslogin01 spike]$ cat ubiquiin.query  
    /home/stivala/public_github/cuda_satabsearch/old/nvcc_src_cuda5/tableauxdistmatrixdb.small.ascii
    T T F
       1UBQ    8
    e  
    OT e  
    LE RT xa 
    PD OS RD xg 
    RT LE RT LS e  
    LE RT LE LS OT e  
    RT LS LS RD PE OS xg 
    PE RT LE RD OT PE RT e  
     0.000 
     4.249  0.000 
    11.944 10.374  1.000 
    17.953 17.702  9.833  3.000 
    11.495 13.762 11.857 10.583  0.000 
    15.612 18.778 17.210 15.380  6.506  0.000 
    14.617 17.265 16.292 20.111 13.272 11.643  3.000 
     8.878 11.115 12.271 12.482  4.567  9.961 15.699  0.000 
    [stivala@icslogin01 spike]$ 
    [stivala@icslogin01 spike]$ ~/cuda_satabsearch/nvcc_src_current/cudaSaTabsearch  -c < ubiquiin.query  > ubiquitin_results_cpu.out
    MAXDIM = 111 (edit in saparams.h and recompile to change)
    Read 1 query structures
    Loading database...
    Loaded 586 db entries (0 order > 96) in 19.839001 ms
    maxstart = 128
    Executing simulated annealing tableaux match kernel on host for query 1UBQ...
    running on host
    host execution time 665.526978 ms
    11.270467 million iterations/sec
    [stivala@icslogin01 spike]$ 
    [stivala@icslogin01 spike]$ head ubiquitin_results_cpu.out 
    # cudaSaTabsearch LTYPE = T LORDER = T LSOLN = F
    # QUERY ID = 1UBQ    
    # DBFILE = /home/stivala/public_github/cuda_satabsearch/old/nvcc_src_cuda5/tableauxdistmatrixdb.small.ascii
    d1kcul1  6 0.6 -1.27278 0.943444
    d1rroa_  4 0.421053 -1.27278 0.943444
    d1weva_  2 0.307692 -1.27278 0.943444
    d2pq3a1  4 0.533333 -1.27278 0.943444
    d1nldl1  8 0.842105 -1.27278 0.943444
    d1lfwa2  5 0.434783 -1.27278 0.943444
    d1lfwa1  7 0.466667 -1.27278 0.943444
    [stivala@icslogin01 spike]$ 
    [stivala@icslogin01 spike]$ srun -p gpu --mem=8G --time=1:0:0 --pty  bash -i
    [stivala@icsnode05 spike]$ ~/cuda_satabsearch/nvcc_src_current/cudaSaTabsearch   < ubiquiin.query > ubiquitin_results_gpu.out 
    MAXDIM = 111 (edit in saparams.h and recompile to change)
    Read 1 query structures
    Loading database...
    Loaded 586 db entries (0 order > 96) in 23.226999 ms
    found 1 CUDA devices
    found modern architecture (compute capability 8.0) device 0: NVIDIA A100-PCIE-40GB
    using device 0: NVIDIA A100-PCIE-40GB
    totalGlobalMem              = 39.5861 GB
    sharedMemPerBlock           = 48 KB
    warpSize                    = 32
    maxThreadsPerBlock          = 1024
    clockRate                   = 1410 MHz
    totalConstMem               = 64 KB
    multiProcessorCount         = 108
    maxThreadsPerMultiProcessor = 2048
    sharedMeMPerMultiprocessor  = 32 KB
    maxBlocksPerMultiProcessor  = 0
    maxstart = 128
    Execution configuration: Grid = (128,1,1) Block = (128,1,1)
    using shared memory for small db structs: YES
    Copying database to device...
    Initialized device RNG with 16384 states (768 KB) in 2.842000 ms
    d_tableaux.pitch == 512 xsize == 96 ysize == 96
    d_distmatrices.pitch == 512 xsize == 384 ysize == 96
    srcPtr.pitch == 96
    distmatrices srcPtr.pitch == 384
    Copied 586 entries to GPU in 22.931000 ms
    XXX c_qn_addr = 0x7fdd815f8500 , qn = 8
    Copying query to constant memory took 0.054000 ms
    Executing simulated annealing tableaux match kernel ( shared memory) on GPU for qid 1UBQ...
    xxx_qn=8
    GPU execution time 5.931000 ms
    1264.677071 million iterations/sec
    [stivala@icsnode05 spike]$ 
    [stivala@icsnode05 spike]$ exit
    [stivala@icslogin01 spike]$ head ubiquitin_results_gpu.out 
    # cudaSaTabsearch LTYPE = T LORDER = T LSOLN = F
    # QUERY ID = 1UBQ    
    # DBFILE = /home/stivala/public_github/cuda_satabsearch/old/nvcc_src_cuda5/tableauxdistmatrixdb.small.ascii
    d1kcul1  11 1.1 0.903563 0.161558
    d1rroa_  5 0.526316 -1.27278 0.943444
    d1weva_  2 0.307692 -1.27278 0.943444
    d2pq3a1  3 0.4 -1.27278 0.943444
    d1nldl1  11 1.15789 0.903563 0.161558
    d1lfwa2  9 0.782609 -1.27278 0.943444
    d1lfwa1  9 0.6 -1.27278 0.943444
    [stivala@icslogin01 spike]$ 
    [stivala@icslogin01 spike]$ cat ubiquitin_results_gpu.out  | grep -v '^#' | sort -k 2,2nr |head
    d1c3ta_  42 5.25 9.60895 2.49513e-06
    d2faza1  41 5.46667 9.60895 2.49513e-06
    d1uela_  36 4.5 7.4326 4.06742e-05
    d3seba2  25 2.77778 3.07991 0.0107511
    d1qnga_  24 2.28571 3.07991 0.0107511
    d1xo7a_  24 2.4 3.07991 0.0107511
    d1ahsa_  20 2 3.07991 0.0107511
    d1ihga2  19 1.80952 0.903563 0.161558
    d1rm6c2  19 2.23529 3.07991 0.0107511
    d1tlya_  19 1.58333 0.903563 0.161558
    [stivala@icslogin01 spike]$ 

shows that the top-scoring structures (and with small p-values) are ubiquitin-like e.g d1c3ta_


This shows how to build the database. Starting with the ASTRAL SCOP 40% sequence identity subset (v2.07) downloaded from https://scop.berkeley.edu/astral/pdbstyle/ver=2.08 which when uncompressed is a hierarchy of 14 323 PDB files (about 4.2GB), we end up with the 22 MB tableaux + distance matrix db file with 14 297 entries (some have no SSEs so are skipped), in ASCII format, sorted by size, for use with cudaSaTabsearch:

    [stivala@icslogin01 ~]$ cat /scratch/stivala/astral/README 
    pdbstyle-sel-gs-bib-40-2.07.tgz Downloaded from 
    https://scop.berkeley.edu/downloads/pdbstyle/pdbstyle-sel-gs-bib-40-2.07.tgz
    with
    
     wget --no-check-certificate https://scop.berkeley.edu/downloads/pdbstyle/pdbstyle-sel-gs-bib-40-2.07.tgz
    
    
    See https://scop.berkeley.edu/astral/pdbstyle/ver=2.08
    
    "If you only need smaller subsets of these files, you can instead download archives containing the PDB-style files corresponding to our two most commonly requested genetic domain sequence subsets: the 40% ID filtered subset (992.1 MB)  ..."
    
    9 Sept 2022
    [stivala@icslogin01 ~]$ 
    
    
    [stivala@icsnode17 astral]$ time python2 ~/cuda_satabsearch/scripts/buildtableauxdb.py -t stride -3 -5 -p none ./pdbstyle-2.07 astral-sel-gs-bib-40-2.07-tableaux.pickle > buildtabldauxdb1.out 2> buildtabldauxdb1.err                        
    
    real    48m9.281s
    user    45m10.867s
    sys     3m4.331s
    [stivala@icsnode17 astral]$ 
    
    /home/stivala/cuda_satabsearch/scripts/buildtableauxdb.py: version is: Revision 3839M, Tue Jun 29 17:03:49 EST 2010
    /home/stivala/cuda_satabsearch/scripts/buildtableauxdb.py: options are: ['-t', 'stride', '-3', '-5', '-p', 'none', './pdbstyle-2.07', 'astral-sel-gs-bib-40-2.07-tableaux.pickle']
    processed 14323 files
    resulting in 14323 db entries
                 14297 tableaux
    with 0 duplicate key errors
    writing tableaux db to astral-sel-gs-bib-40-2.07-tableaux.pickle...
    done.
    [stivala@icsnode17 astral]$ 
    [stivala@icsnode17 astral]$ time python2 ~/cuda_satabsearch/scripts/buildtableauxdb.py -t stride -3 -5 -p none -d ./pdbstyle-2.07 astral-sel-gs-bib-40-2.07-distmatrices.pickle > buildtabldauxdb2.out 2> buildtabldauxdb2.err                 
    
    real    30m53.349s
    user    28m40.969s
    sys     2m35.274s
    [stivala@icsnode17 astral]$ 
    [stivala@icsnode17 astral]$ cat buildtabldauxdb2.out 
    /home/stivala/cuda_satabsearch/scripts/buildtableauxdb.py: version is: Revision 3839M, Tue Jun 29 17:03:49 EST 2010
    /home/stivala/cuda_satabsearch/scripts/buildtableauxdb.py: options are: ['-t', 'stride', '-3', '-5', '-p', 'none', '-d', './pdbstyle-2.07', 'astral-sel-gs-bib-40-2.07-distmatrices.pickle']
    processed 14323 files
    resulting in 14323 db entries
                 14297 SSE distance matrices
    with 0 duplicate key errors
    writing SSE distance matrix db to astral-sel-gs-bib-40-2.07-distmatrices.pickle...
    done.
    [stivala@icsnode17 astral]$ 
    [stivala@icsnode17 astral]$ time python2 ~/cuda_satabsearch/scripts/convdb2.py -s astral-sel-gs-bib-40-2.07-tableaux.pickle astral-sel-gs-bib-40-2.07-distmatrices.pickle > tableauxdistmatrixdb-sel-gs-bib-40-2.07.sorted.ascii
    loading tableaux...
    loading distance matrices...
    writing ASCII format tableaux+distmatrices...
    sorting database...
    wrote 14297 tableaux+distmatrices for 14297 entries
    
    real    0m15.047s
    user    0m14.934s
    sys     0m1.154s
    [stivala@icsnode17 astral]$ 
    

And now we can do a query in this database (I just changed the first line of the query to have the filename of the new database file just created):


    [stivala@icsnode05 astral]$ cat ubiquitin.query                                
    /scratch/stivala/astral/tableauxdistmatrixdb-sel-gs-bib-40-2.07.sorted.ascii   
    T T F
       1UBQ    8                                                                   
    e
    OT e                                                                           
    LE RT xa                                                                       
    PD OS RD xg
    RT LE RT LS e
    LE RT LE LS OT e
    RT LS LS RD PE OS xg
    PE RT LE RD OT PE RT e                                                         
     0.000                                                                         
     4.249  0.000                                                                  
    11.944 10.374  1.000
    17.953 17.702  9.833  3.000
    11.495 13.762 11.857 10.583  0.000                                             
    15.612 18.778 17.210 15.380  6.506  0.000                                      
    14.617 17.265 16.292 20.111 13.272 11.643  3.000                               
     8.878 11.115 12.271 12.482  4.567  9.961 15.699  0.000                        
    
    stivala@icsnode05 astral]$ ~/cuda_satabsearch/nvcc_src_current/cudaSaTabsearch   < ubiquitin.query > ubiquitin_results_gpu.out                               
    
    MAXDIM = 111 (edit in saparams.h and recompile to change)
    Read 1 query structures
    Loading database...                                                                                                                                           
    Tableau d4c2ma_ order 115 is too large (max is 111)
    WARNING: excluded database structure d4c2ma_ as it is too large
    Tableau d4g7hd_ order 131 is too large (max is 111)
    WARNING: excluded database structure d4g7hd_ as it is too large
    WARNING: skipped 2 tableaux of order > 111
    Loaded 14295 db entries (4 order > 96) in 774.919983 ms                                                                                                       
    found 1 CUDA devices
    found modern architecture (compute capability 8.0) device 0: NVIDIA A100-PCIE-40GB                                                                            
    using device 0: NVIDIA A100-PCIE-40GB
    totalGlobalMem              = 39.5861 GB                                                                                                                      
    sharedMemPerBlock           = 48 KB
    warpSize                    = 32                                                                                                                              
    maxThreadsPerBlock          = 1024
    clockRate                   = 1410 MHz                                                                                                                        
    totalConstMem               = 64 KB                                                                                                                           
    multiProcessorCount         = 108
    maxThreadsPerMultiProcessor = 2048                                                                                                                            
    sharedMeMPerMultiprocessor  = 32 KB                                                                                                                           
    maxBlocksPerMultiProcessor  = 0
    maxstart = 128                                                                                                                                                
    Execution configuration: Grid = (128,1,1) Block = (128,1,1)
    using shared memory for small db structs: YES
    Copying database to device...
    Initialized device RNG with 16384 states (768 KB) in 2.826000 ms
    d_tableaux.pitch == 512 xsize == 96 ysize == 96
    d_distmatrices.pitch == 512 xsize == 384 ysize == 96                                                                                                          
    srcPtr.pitch == 96                                                                                                                                            
    distmatrices srcPtr.pitch == 384                                                                                                                              
    Copied 14291 entries to GPU in 96.946999 ms                                                                                                                   
    XXX c_qn_addr = 0x7f20df5f8500 , qn = 8
    Copying query to constant memory took 0.054000 ms
    Executing simulated annealing tableaux match kernel ( shared memory) on GPU for qid 1UBQ...
    xxx_qn=8
    GPU execution time 96.112999 ms
    1903.226431 million iterations/sec
    Copying large structure database to device...
    d_tableaux.pitch == 512 xsize == 111 ysize == 111
    d_distmatrices.pitch == 512 xsize == 444 ysize == 111
    srcPtr.pitch == 111
    distmatrices srcPtr.pitch == 444
    Copied 4 large entries to GPU in 1.330000 ms
    XXX c_qn_addr = 0x7f20df341400 , qn = 8
    Copying query to constant memory took 0.043000 ms
    Executing simulated annealing tableaux match kernel (no shared memory) on GPU for qid 1UBQ...
    xxx_qn_noshared=8
    GPU (no shared memory) execution time 5.231000 ms
    9.787804 million iterations/sec
    [stivala@icsnode05 astral]$ 
    [stivala@icsnode05 astral]$ cat ubiquitin_results_gpu.out  | grep -v '^#' | sort -k 2,2nr |head
    d3m62b_  54 6.75 11.7853 1.53059e-07
    d2bwfa_  53 6.625 11.7853 1.53059e-07
    d2faza1  52 6.5 11.7853 1.53059e-07
    d3a9ja_  52 6.5 11.7853 1.53059e-07
    d3zdmc_  51 6.375 11.7853 1.53059e-07
    d4eewa_  48 6 11.7853 1.53059e-07
    d3kuza1  47 5.22222 9.60895 2.49513e-06
    d2wyqa_  46 5.75 9.60895 2.49513e-06
    d4icva_  44 5.5 9.60895 2.49513e-06
    d1yqba1  42 5.25 9.60895 2.49513e-06
    
