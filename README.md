# pal-MEM
Identifying inverted repeats in large genomic read files using efficient computation of Maximal Exact Matches (MEMs)

## DESCRIPTION

pal-MEM finds inverted terminal repeats (ITRs) in large genomes, including metagenomes. It is modified from the program [E-MEM](https://github.com/lucian-ilie/E-MEM) by N. Khiste, L. Ilie [E-MEM: efficient computation of maximal exact matched for very large genomes](http://bioinformatics.oxfordjournals.org/content/31/4/509.short)

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

pal-MEM outputs a tab-delimited file and two fasta files. The tab-delimited file contains the original sequence names of reads containing the ITR pair with the first and second columns representing the first and second read of the pair, respectively.

    Seq1409_ERR589353.1592_FCC4C01ACXX:7:1101:11563:2918#GCGGAACT/1_f1	Seq35410689_ERR589541.11558754_FCD14R7ACXX:2:1308:15750:53071#CTCGTCCG/1_f1
    Seq3649_ERR589353.4045_FCC4C01ACXX:7:1101:6219:4203#GCGGAACT/1_f1	Seq29675265_ERR589541.5713831_FCD14R7ACXX:2:1204:4038:103263#CTCGTCCG/1_f1
    Seq5633_ERR589353.6277_FCC4C01ACXX:7:1101:7124:5063#GCGGAACT/1_f1	Seq17102433_ERR589353.17765993_FCC4C01ACXX:7:2205:21103:80820#GCGGAACT/1_f1
    ...

One fasta file contains reads with inverted terminal repeats (ITRs). The sequence name is reported with "LCoord" and "RCoord" values representing the coordinates of the ITR. For example, the ITR pair of length 41 is situated in the first and second read from 28th to 68th and from 31st to 71st nucleotide, respectively.

    >Seq1409_ERR589353.1592_FCC4C01ACXX:7:1101:11563:2918#GCGGAACT/1_f1_LCoord_28_RCoord_68
    AAGGTTTGATCCTTCTTATCGTCATTATCGAAGGTCTTAGGTCCTGCGAACTCATCCGTTACAACTTCATAGCCCTTGTCTGTCAATTCTTTCAGACGGG
    >Seq35410689_ERR589541.11558754_FCD14R7ACXX:2:1308:15750:53071#CTCGTCCG/1_f1_LCoord_31_RCoord_71
    CGTCTGAAAGAATTGACAGACAAGGGTTACGAAGTTGTAACGGATGAGTTCGCAGGACCTAAGACCTTCGACAATGATGATAAGAAGGATCAAACCTTCA

The other fasta file contains the rest of the reads that do not contain ITRs with their original headers.

## OPTIONS

The program can be run in both serial and parallel mode. The parallel mode has an advantage in terms of time with respect to serial mode. The options for pal-mem are:

-l set the minimum length of a match. Default: 24

-k set the k-mer length. Default: 15

-d set the split size. Default: 1. The option -d is used for splitting the sequences into two or more parts. By default this value is set to 1, which means no splitting. This option with value >1 will reduce the overall memory requirement of the program. For large genomic files (such as metagenomic files over 4 Gigabytes), it is recommended to increase d when RAM is limited, otherwise the program will fail (see example below).

-t number of threads. Default: 1. The option -t is used for running the program in parallel mode. The default value is set to 1, which means serial mode. This option with value > 1 will reduce overall running time of the program.

-h show possible options

## EXAMPLE
To get ITRs with a minimum length of 30 from paired-end metagenomic fasta files with a k-mer length of 18, a split size of 20 on a machine with 8 threads:
```
pal-mem -f1 example_1.fasta -f2 example_1.fasta -o example -l 30 -k 18 -d 20 -t 8
```

## INSTALLATION

### Requirements:
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

### For Mac OS X
- By default, Mac OS X uses clang as its compiler for C++. If not already installed, you must install the GCC compiler and locate the executable, which should be */usr/local/bin/g++-(version)*. Use this as the CXX option for cmake.
```
CXX=/usr/local/bin/g++-9 cmake ..
make
```

### For Linux
```
cmake ..
make
```

### For Linux (where boost is not installed with root access)
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
