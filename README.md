# DARRC: Dynamic Alignment-free and Reference-free Read Compression

This repository contains the source code of DARRC, a new alignment-free and reference-free compression method. It addresses the problem of pan-genome compression by encoding the sequences of a pan-genome as a guided de Bruijn Graph. The novelty of this method is its ability to incrementally update DARRC archives with new genome sequences without full decompression of the archive. DARRC can compress both single-end and paired-end read sequences of any length using all symbols of the IUPAC nucleotide code. On a large *P. aeruginosa* dataset, our method performs as little as 0.204 bits per base in single-end mode.

## Dependencies

In order to compile and run DARRC, you need a machine running a 64 bits Linux or Mac OS operating system, POSIX compliant and accepting SSE4 instructions.

DARRC depends on four external libraries: Judy Arrays (http://judy.sourceforge.net), Jemalloc (https://github.com/jemalloc/jemalloc), 7zip (http://www.7-zip.org/) and the Bloom Filter Trie (https://github.com/GuillaumeHolley/BloomFilterTrie).

All of them can be downloaded and installed by following the instructions on their respective websites. It is however most likely that at least a few of them are available via a package manager for your operating system.

If you operating system is Ubuntu or Debian, you can install the Judy Arrays, Jemalloc and 7zip with the following command:
```
sudo apt-get install libjudydebian1 libjudy-dev libjemalloc1 libjemalloc-dev p7zip
```

If you operating system is Mac OS, Jemalloc and 7zip can be easily downloaded and installed via Homebrew:
```
brew install jemalloc
brew install p7zip
```

## Compilation and installation

DARRC compiles with GNU GCC and G++ (compatibility with Clang is in preparation). DARRC successfully compiles and runs on Ubuntu 14.04, 15.04 and 16.10.

**Important**: DARRC creates and opens simultaneously around a thousand of temporary files, which can be above the default limit of your OS, especially if this OS is macOS. Make sure before running DARRC that you have changed the default limit to an upper bound, using the following command:
```
ulimit -n 2048
```

### Linux

On Linux, you can verify the presence of gcc and g++ on your system with:
```
gcc -v
g++ -v
```

If not present (unlikely), they can be installed for Ubuntu or Debian with:
```
sudo apt-get install build-essential
```

Compiling DARRC should then be as simple as:
```
cd <DARRC_directory>
./configure
make
make install
```

You can also install DARRC in a specific directory (for example because you are not root on the machine you are using) with:
```
cd <DARRC_directory>
./configure --prefix=<your_directory>
make
make install
```

Make sure that your environment variables are all set with <your_directory>.

### Mac OS

**Disclaimer**: DARRC compiles but was not tested with Mac OS.

For Mac OS, you will need the "real" GCC and G++ (not the Clang interface for GCC or G++, default on macOS). Both can be installed via Homebrew:
```
brew install gcc-x
brew install g++-x
```

in which *x* is the latest major version of GCC and G++.

Compiling DARRC should then be as simple as:
```
cd <DARRC_directory>
./configure CC=gcc-x
make
make install
```

in which *x* is the version of GCC and G++ that you installed via Homebrew or other.

To install DARRC in a specific directory, see Linux compilation and installation.

## Binary usage:
```
darrc <action parameter> <action-specific parameters> <general parameters>

<action parameter>:

	-c	[--compress]		compression
	-u	[--update]		update
	-d	[--decompress]		decompression

	-v	[--version]		print version info
	-h	[--help]		print help info

<compression parameters>:

	-k	[--kmer]	arg	length of k-mers, must be either 18, 27, 36 (default), 45, 54 or 63
	-o	[--overlap]	arg	length of k-mers overlap, must be between 8 and 15 (default: 11)

<compression and update parameters>:

	-min	[--minimizer]	arg	length of minimizers, must be between 8 and 10 (default: 9)
	-mis	[--mismatch]	arg	number of mismatches allowed during merging (default: 5)
	-1	[--mate1]	arg	input FASTA/Q file: single-end reads or first mate of paired-end reads
	-2	[--mate2]	arg	input FASTA/Q file: second mate of paired-end reads
	-l1	[--listmate1]	arg	list of input FASTA/Q files: single-end reads or first mate of paired-end reads
	-l2	[--listmate2]	arg	list of input FASTA/Q files: second mate of paired-end reads

<decompression parameters>: None

<general parameters>:

	-g	[--graph]	arg	graph filename prefix (output for compression, input otherwise)
	-m	[--meta]	arg	meta filename prefix (output for compression/update, input otherwise)
	-lm	[--listmeta]	arg	list of meta filename prefixes (output for compression/update, input otherwise)
	-dir	[--directory]	arg	compression and update:	directory for temporary files
					decompression: directory for temporary and output files
```

## Usage examples


### Compression

Compressing single-end reads in *strain_u.fastq*:
```
darrc -c -1 strain_u.fastq -g species -m strain_u
```

Compressing paired-end reads (first mates in file *mate_1.fastq* and second mates in file *mate_2.fastq*):
```
darrc -c -1 mate_1.fastq -2 mate_2.fastq -g species -m strain_u
```

Both commands will produce a graph file named *species.graph.darrc* and a meta file named *strain_u.meta.darrc*.
List of files to compress can be also used with the *-l1* and/or *-l2* parameters, similarly as the *-1* and *-2* parameters, combined with the *-lm* parameter instead of *-m*. Note you can input single-end read files and paired-end read files in the list of files to compress. Example:
```
darrc -c -l1 mates_1.txt -l2 mates_2.txt -g species -m strains.txt
```

*mates_1.txt*
```
strain_u_1.fastq
strain_v.fastq
strain_w_1.fastq
```

*mates_2.txt*
```
strain_u_2.fastq

strain_w_2.fastq
```

*strains.txt*
```
strain_u
strain_v
strain_w
```

In *mates_2.txt*, the second line is empty because *strain_v.fastq* in *mates_1.txt* corresponds to single-end reads. 

### Update

Updating a compressed archive works similarly as the compression except that the parameter *-g* is used for the (input) graph to update. Assuming a previously created graph file *species.graph.darrc*, updating it with single-end reads is:
```
darrc -u -1 strain_x.fastq -g species -m strain_x
```

Updating with paired-end reads (first mates in file *mate_1.fastq* and second mates in file *mate_2.fastq*):
```
darrc -u -1 mate_1.fastq -2 mate_2.fastq -g species -m strain_x
```

Both commands will update *species.graph.darrc* and create a new meta file *strain_x.meta.darrc*. List of files to compress (*-l1* and/or *-l2* + *lm* parameters) can be used as well.

### Decompression

Decompressing reads of file *strain_v.fastq* compressed in graph file *species.graph.darrc* and meta file *strain_v.meta.darrc* is:
```
darrc -d -g species -m strain_v
```

This will produce an output file *strain_v.reads*. If used for a paired-end meta file *strain_u.meta.darrc*, it will output two files *strain_u.reads_1* and *strain_u.reads_2*.

## Contact

For any question, feedback or problem, please contact me at gholley{aT}cebitec{d0t}uni-bielefeld{d0t}de
