# polyploid-potato-assembly
Code for assembly approach presented in "Haplotype-resolved assembly of a tetraploid potato genome using long reads and low-depth offspring data"

### 1. Dosage estimation of the nodes in the hifiasm assembly graph:
  run `snakemake` in the `coverage-analysis` directory. Requires `minimap2` and `samtools`.

### 2. K-mer analysis:
####   Installation:
   Note: Requires an installation of the `jellyfish` package. If not installed yet, you can install it via `conda install jellyfish`
   
   `git clone git@github.com:rebeccaserramari/polyploid-potato-assembly.git`
   
   `cd polyploid-potato-assembly/kmer-counting`
   
   `mkdir build; cd build; cmake ..; make`
####  Running k-mer counting procedure
1. To run the full procedure, including finding unique k-mers in \<targetfile\>, counting the found unique k-mers in a set of sequences samples, and merging the resulting files:

    run `snakemake` within the `kmer-counting` directory.
    
Make sure to update the config files accordingly!    
2. To run the first step individually, i.e. find k-mers of length \<len\> that are uniquely present in \<targetfile\> and not in \<comparisonfile\>:

    `./polyassembly_findkmers find_kmers -r <targetfile> -s <comparisonfile> -k <kmerfile> -l <len>`

The resulting k-mers are stored in \<kmerfile\>.
  
3. To run the second step individually, i.e. count the unique k-mers in \<samplefile\>:
  
    `/polyassembly_findkmers count_kmers -s <samplefile> -k <kmerfile> -c <output> -l <len>`
  
### 3. Phased clustering

To run the full clustering procedure:

  run `snakemake` in the `cluster-phasing` directory.
 
