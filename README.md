# pal-MEM
Identifying inverted repeats in large genomic read files using efficient computation of Maximal Exact Matches (MEMs)

## DESCRIPTION

pal-MEM finds inverted terminal repeats (ITRs) in large genomes, including metagenomes. It is modified from the program [E-MEM](https://github.com/lucian-ilie/E-MEM) by N. Khiste, L. Ilie [E-MEM: efficient computation of maximal exact matched for very large genomes](http://bioinformatics.oxfordjournals.org/content/31/4/509.short)

## INSTALLATION

### Download binaries
[Available here](https://github.com/blue-moon22/pal-MEM/releases) for Mac OS (built on Catalina v10.15.6) and Linux

### Or install from source
#### Requirements
- pal-MEM requires a 64-bit system
- [Boost](https://www.boost.org/)
- [cmake](https://cmake.org/download/)

First clone the repository and go into the pal-MEM directory, then follow steps below depending on your operating system.
```
git clone https://github.com/blue-moon22/pal-MEM.git
cd pal-MEM
mkdir build
cd build
```

#### For Mac OS
- By default, Mac OS uses clang as its compiler for C++. If not already installed, you must install cmake, the GCC compiler and boost. The easiest way is to install using [Homebrew](https://brew.sh/). Once brew is installed:
```
brew install cmake boost gcc
CXX=/usr/local/Cellar/gcc/<gcc_version>/bin/g++-<version> cmake ..
make
```

#### For Linux
```
cmake ..
make
```

#### For Linux (where boost is not installed with root access)
```
boost_location=<wherever boost installed containing lib and include directories>
LD_LIBRARY_PATH=${boost_location}/lib:$LD_LIBRARY_PATH
LIBRARY_PATH=${boost_location}/lib:$LIBRARY_PATH
INCLUDE=${boost_location}/include/boost:$INCLUDE
C_INCLUDE_PATH=${boost_location}/include/boost:$C_INCLUDE_PATH
INCLUDE_PATH=${boost_location}/include/boost:$INCLUDE_PATH

BOOST_INCLUDEDIR=${boost_location}/include BOOST_LIBRARYDIR=${boost_location}/lib cmake ..
make
```
- If using HPC when running pal-mem, ensure paths are set above.


### USAGE
```
pal-mem  -f1 <paired-end fasta file 1>  -f2 <paired-end fasta file 2>  -o <output prefix>  [options]
```
OR
```
pal-mem  -fu <single fasta file>  -o <output prefix>  [options]
```

Type *pal-mem -h* for a list of options.

### OUTPUT

pal-MEM outputs a tab-delimited file and two fasta files (if a single fasta file is specified) or four fasta files (if paired fasta files are specified). The tab-delimited file with suffix `_IR.tab` contains the original sequence names of reads containing the inverted repeats (IRs) with the first and second columns representing the first and second read of the pair, respectively. The sequence name is reported with "LCoord" and "RCoord" values representing the coordinates of the ITR.

    Seq2691_ERR589346.2764_FCC4C01ACXX:6:1101:14662:3637#AAGTCTCT/1_f1_LCoord_38_RCoord_64	Seq422_ERR589346.468_FCC4C01ACXX:6:1101:11875:2371#AAGTCTCT/1_f1_LCoord_30_RCoord_56
    Seq4972_ERR589346.5141_FCC4C01ACXX:6:1101:21190:4913#AAGTCTCT/1_f1_LCoord_45_RCoord_60	Seq3747_ERR589346.3872_FCC4C01ACXX:6:1101:3886:4252#AAGTCTCT/1_f1_LCoord_51_RCoord_66
    ...

The fasta file(s) with suffix `_IR.fasta` (if single file) or `_IR_1.fasta` and `_IR_2.fasta` (if paired files) contains reads with inverted repeats (IRs).

    >Seq2691_ERR589346.2764_FCC4C01ACXX:6:1101:14662:3637#AAGTCTCT/1
    AGCCTCTGGTAGGGCTAATTGCACCAAATGTCCCAATGCCCAAGTAACGGCATAATCATTACCTTGCAAATACCCATCTTGTTTTTCGGTTGCCCCCAC
    >Seq422_ERR589346.468_FCC4C01ACXX:6:1101:11875:2371#AAGTCTCT/1
    GCAAGTGAAAAACAAGATGGATATTTGTTAGGTAATGATTATGCCGTTACTTGGGCTCTAGGGCATTTGGTGCAATTAGCCCTCCCAGAGGCTTATGGTT

The fasta file(s) with suffix `_discord_non_IR.fasta` (if single file) or `_discord_non_IR_1.fasta` and `_discord_non_IR_2.fasta` (if paired files) contains reads that do not contain IRs but are paired to reads that do contain IRs. For example, if the reads paired to the above do not contain IRs, then the second file `_discord_non_IR_2.fasta` would contain:

    >Seq2691_ERR589346.2764_FCC4C01ACXX:6:1101:14662:3637#AAGTCTCT/2_f2
    CAATCCTCAAAAGGACATTACTGATAAAGTTTCTTCTACCAAACAAAAAGCTGAAACTTCTAAAGCCAAAGAAGAAAAACAACCTCAAAAGCAATCAGAA
    >Seq422_ERR589346.468_FCC4C01ACXX:6:1101:11875:2371#AAGTCTCT/2_f2
    CCTTTGAAAAGGTTTTTGACAGTGCAAATATTCATAGATATAGCGGAAGATCAGTTCTCCTTCACGACCTGCATCGGTGGCTACAATAATGGAACTACAC

## OPTIONS

The program can be run in both serial and parallel mode. The parallel mode has an advantage in terms of time with respect to serial mode. The options for pal-mem are:

-l set the minimum length of a match. Default: 24

-m set the maximum length of a match. Default: 100

-k set the k-mer length. Default: 15

-t number of threads. Default: 1. The option -t is used for running the program in parallel mode. The default value is set to 1, which means serial mode. This option with value > 1 will reduce overall running time of the program.

-h show possible options

## EXAMPLE
To get ITRs with a minimum length of 30 and maximum length of 50 from paired-end metagenomic fasta files with a k-mer length of 18 on a machine with 8 threads:
```
pal-mem -f1 example_1.fasta -f2 example_1.fasta -o example -l 30 -m 50 -k 18 -t 8
```
