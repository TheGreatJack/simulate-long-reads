# Nuclear + Mitochondrial lonng read simulator


This is based on another [repository](https://github.com/BCCDC-PHL/simulate-long-reads). And is mainly a wrapper of a [long read simulator](https://github.com/rrwick/Badread). **This is just a prototype and errors are expected.** 

The main idea is that you take the nuclear genome and mitochondrial genome and you run Badread on each of the genomes. With varying degrees of coverage for each genome. Then you combine the resulting fastq files. The rest is mapping the reads back into the genome and getting some statistics.

To run this wrapper isn't necessary to download or clone the repo. Nextflow automatically downloads the repo for you, it also installs the needed conda environment for the run. Nextflow may be downloaded via conda. 

To run simply use:

```
nextflow run TheGreatJack/simulate-long-reads -r main -with-conda -profile conda --cache ~/miniconda3/envs --nuc_assembly $path_to_nuclear_fasta --mito_assembly $path_to_mitochondrial_fasta --nuc_depth $coverage_of_genome --mito_depth $coverage_of_mitochondria --outdir $out_dir

```

The main result files will be in the folder name of nuclear genome. 