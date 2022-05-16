# polyploid-potato-assembly
Code for assembly approach presented in "Haplotype-resolved assembly of a tetraploid potato genome using long reads and low-depth offspring data"

### 1. Dosage estimation of the nodes in the hifiasm assembly graph:
  run `snakemake` in the `coverage-analysis` directory. Requires `minimap2` and `samtools`.

### 2. K-mer analysis:
####   Installation:
   Note: Requires an installation of the `jellyfish` package. If not installed yet, you can install it via `conda install jellyfish`
   
   `git clone git@github.com:rebeccaserramari/polyploid-potato-assembly.git`
   
   `cd polyploid-potato-assembly`
   
   `mkdir build; cd build; cmake ..; make`
####  Running k-mer counting procedure
1. Find k-mers of length \<len\> that are uniquely present in \<targetfile\> and not in \<samplefile\>:

    `./polyassembly_findkmers find_kmers -r <targetfile> -s <samplefile> -k <kmerfile> -l <len>`
  
2. Count k-mers from previously computed <kmerfile> in short reads of progeny samples:
  
    run `snakemake` within the `kmer-counting` directory.
  
### 3. Run the clustered phasing procedure

 
