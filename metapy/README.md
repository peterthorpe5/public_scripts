[![THAPBI-pycits TravisCI build status](https://api.travis-ci.org/widdowquinn/THAPBI-pycits.svg?branch=master)](https://travis-ci.org/widdowquinn/THAPBI-pycits/branches)
[![codecov.io coverage status](https://img.shields.io/codecov/c/github/widdowquinn/THAPBI-pycits.svg)](https://codecov.io/gh/widdowquinn/THAPBI-pycits)
[![Code Health](https://landscape.io/github/widdowquinn/THAPBI-pycits/master/landscape.svg?style=flat)](https://landscape.io/github/widdowquinn/THAPBI-pycits/master)

* [Travis-CI page](https://travis-ci.org/widdowquinn/THAPBI-pycits/branches)
* [`codecov.io` page](https://codecov.io/gh/widdowquinn/THAPBI-pycits)
* [`landscape.io` page](https://landscape.io/github/widdowquinn/THAPBI-pycits/1)

# README.py - `THAPBI-pycits`
This repository is for development of ITS1-based diagnostic/profiling tools for the THAPBI Phyto-Threats project, funded by BBSRC.

# USAGE:

SEE ``USAGE.md`` on how to use the tool


# DEVELOPER NOTES

## Python style conventions

In this repository, we're trying to keep to the Python [PEP8 style convention](https://www.python.org/dev/peps/pep-0008/), the [PEP257 docstring conventions](https://www.python.org/dev/peps/pep-0257/), and the [Zen of Python](https://www.python.org/dev/peps/pep-0020/). To help in this, a pre-commit hook script is provided in the `git_hooks` subdirectory that, if deployed in the `Git` repository, checks Python code for PEP8 correctness before permitting a `git commit` command to go to completion.

If the `pep8` module is not already present, it can be installed using `pip install pep8`

Whether you choose to use this or not, the `THAPBI-pycits` repository is registered with `landscape.io`, and the "health" of the code is assessed and reported for every repository push.

* [`landscape.io` page](https://landscape.io/github/widdowquinn/THAPBI-pycits)

### Installing the `git hook`

To install the pre-commit hook:

1. clone the repository with `git clone https://github.com/widdowquinn/THAPBI` (you may already have done this)
2. change directory to the root of the repository with `cd THAPBI-pycits`
3. copy the pre-commit script to the `.git/hooks` directory with `cp git_hooks/pre-commit .git/hooks/`

## Using a virtual environment with the repository

In the root directory of the repository:

```
$ virtualenv -p python3.5 venv-THAPBI-pycits
$ source venv-THAPBI-pycits/bin/activate
<activity>
$ deactivate
```

# INSTALLATION

## Dependencies: Python modules

All Python module dependencies are described in `requirements.txt` and can be installed using

```
pip install -r requirements.txt
```

There may be issues with `biom-format` and `biopython` installations due to ordering of module installation. If this is the case for you, then it might be solved by installing `numpy` at the command-line first, with:

```
pip install numpy
pip install -r requirements.txt
```

## Dependencies: Third-party applications


### `pear`

* [home page](http://sco.h-its.org/exelixis/web/software/pear/)

`pear` is a paired-end read merger, used by the pipeline to merge ITS paired-end reads into a single ITS sequence. It is available from the [`pear` home page](http://sco.h-its.org/exelixis/web/software/pear/) as a precompiled executable that can be placed in your `$PATH`, and it can be installed on the Mac with [Homebrew](http://brew.sh/) and [homebrew-science](https://github.com/Homebrew/homebrew-science), using:
If you are having trouble on some Linux system due to bzip+zlib issues, a precompiled version is availble: https://github.com/xflouris/PEAR/blob/master/bin/pear-0.9.5-bin-64
```
brew install pear
```
please rename the binary to pear:
```
wget http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.10-bin-64.tar.gz
tar -zxvf pear-0.9.10-bin-64.tar.gz
cp pear-0.9.10-bin-64/pear-0.9.10-bin-64 pear-0.9.10-bin-64/pear
put this in your PATH
```

### `Trimmomatic`

* [home page](http://www.usadellab.org/cms/?page=trimmomatic)

`Trimmomatic` is used to trim and quality-control the input reads. `pycits` expects `Trimmomatic` to be available at the command-line as `trimmomatic`. You can check if the tool is installed this way with the command:

```
which trimmomatic
```

To obtain `Trimmomatic` with this installation type on Linux systems, you can use:

```
apt-get install trimmomatic
```

and on the Mac (with [Homebrew](http://brew.sh/) and [homebrew-science](https://github.com/Homebrew/homebrew-science)):

```
brew install trimmomatic
```

If you have downloaded the Java `.jar.` file from [`trimmomatic`'s home page](http://www.usadellab.org/cms/?page=trimmomatic), you can wrap the `.jar` file with a Bash script called `trimmomatic` in your `$PATH`, such as

```
#!/bin/bash
exec java -jar $TRIMMOMATIC "$@"
```

where `$TRIMMOMATIC` is the path to your `trimmomatic .jar` file.


### `muscle`

* [home page](http://www.drive5.com/muscle/muscle.html)

`muscle` MUSCLE is a program for creating multiple alignments of amino acid or nucleotide sequences [`muscle` download page]http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux32.tar.gz as a precompiled executable that can be placed in your `$PATH` when decompressed. Please rename to muscle

```
wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux32.tar.gz
tar -zxvf muscle3.8.31_i86linux32.tar.gz
cp muscle3.8.31_i86linux32 muscle
put this in your PATH
```


### `fastqc`

* [home page](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

`fastqc` fastqc A quality control tool for high throughput sequence data. [`fastqc` download page]http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux32.tar.gz as a precompiled executable that can be placed in your `$PATH` when decompressed.

```
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip
cd FastQC
chmod 755 fastqc
put this in your PATH
```

### `spades`

* [home page](http://bioinf.spbau.ru/spades)

`spades` spades is an assembly program which we use for error correction. [`spades` download page]http://spades.bioinf.spbau.ru/release3.9.1/SPAdes-3.9.1-Linux.tar.gz as a precompiled executable that can be placed in your `$PATH` when decompressed.

```
wget http://spades.bioinf.spbau.ru/release3.9.1/SPAdes-3.9.1-Linux.tar.gz
tar -zxvf SPAdes-3.9.1-Linux.tar.gz
cd ./SPAdes-3.9.1-Linux/bin/
put this in your PATH
```

### `flash`

* [home page](https://ccb.jhu.edu/software/FLASH/)

`flash` flash is a pair end read assembly program. [`flash` download page]https://sourceforge.net/projects/flashpage/files/FLASH-1.2.11.tar.gz as a precompiled executable that can be placed in your `$PATH` when decompressed.

```
wget https://sourceforge.net/projects/flashpage/files/FLASH-1.2.11.tar.gz
tar -zxvf FLASH-1.2.11.tar.gz
cd FLASH-1.2.11
put this in your PATH
```

### `swarm`

* [home page](https://ccb.jhu.edu/software/FLASH/)

`swarm` swarm is a clustering program. [`swarm` download page]https://github.com/torognes/swarm 

```
git clone https://github.com/torognes/swarm.git
cd swarm/src
make
put PATH_TO/swarm/bin in your PATH
```

### `blastclust`

* [home page](https://ccb.jhu.edu/software/FLASH/)

`blastclust` blastclust is a clustering program. This is not essential to download. [`blastclust` download page]ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/2.2.26/blast-2.2.26-x64-linux.tar.gz as a precompiled executable that can be placed in your `$PATH` when decompressed.

```
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/2.2.26/blast-2.2.26-x64-linux.tar.gz
tar -zxvf blast-2.2.26-x64-linux.tar.gz
put this in your PATH
```
  
### `cd-hit`

* [home page](https://github.com/weizhongli/cdhit)

`cd-hit` cd-hit is a clustering program. [`cd-hit` download page]https://github.com/weizhongli/cdhit

```
git clone git@github.com:weizhongli/cdhit.git
cd cdhit
make
put PATH_TO/cdhit in your PATH
```

### `bowtie_2.2.5`

* [home page](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

`bowtie_2.2.5` bowtie_2.2.5 is a read mapping program. [`bowtie_2.2.5` download page]https://depot.galaxyproject.org/package/linux/x86_64/bowtie2/bowtie2-2.2.5-Linux-x86_64.tar.gz
We have use the pre compiled binary to reduce difficulty in getting this to work.

```
mkdir bowtie_2.2.5
cd bowtie_2.2.5
wget https://depot.galaxyproject.org/package/linux/x86_64/bowtie2/bowtie2-2.2.5-Linux-x86_64.tar.gz
tar -zxvf bowtie2-2.2.5-Linux-x86_64.tar.gz
put this in your path
export PATH=$HOME/bowtie_2.2.5/bin/:$PATH
```

### `samtools_1.2`

* [home page](http://www.htslib.org/doc/samtools.html)

`samtools_1.2` samtools_1.2 is a program to do lots of thing with sam/ bam files. [`samtools_1.2` download page]https://depot.galaxyproject.org/package/linux/x86_64/samtools/samtools-1.2-Linux-x86_64.tgz
We have use the pre compiled binary to reduce difficulty in getting this to work.

```
mkdir samtools_1.2
cd samtools_1.2
wget https://depot.galaxyproject.org/package/linux/x86_64/samtools/samtools-1.2-Linux-x86_64.tgz
tar -zxvf samtools-1.2-Linux-x86_64.tgz
put this in your path
export PATH=$HOME/samtools_1.2/bin/:$PATH
```

### `vsearch`

* [home page](https://github.com/torognes/vsearch)

`vsearch` vsearch Versatile open-source tool for metagenomics. [`vsearch` download page]https://github.com/torognes/vsearch

```
wget https://github.com/torognes/vsearch/releases/download/v2.4.0/vsearch-2.4.0-linux-x86_64.tar.gz
tar xzf vsearch-2.4.0-linux-x86_64.tar.gz
```


### More information

* Git hooks (`git`): [https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks](https://git-scm.com/book/en/v2/Customizing-Git-Git-Hooks)
* Git hooks (tutorial): [http://githooks.com/](http://githooks.com/)
* PEP8: [https://www.python.org/dev/peps/pep-0008/](https://www.python.org/dev/peps/pep-0008/)
* PEP257: [https://www.python.org/dev/peps/pep-0257/](https://www.python.org/dev/peps/pep-0257/)
* Zen of Python (PEP20): [https://www.python.org/dev/peps/pep-0020/](https://www.python.org/dev/peps/pep-0020/)


