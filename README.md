[![Travis CI](https://travis-ci.org/wtsi-hpag/HiLine.svg?branch=master)](https://travis-ci.org/github/wtsi-hpag/HiLine)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/hiline/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/hiline/badges/downloads.svg)](https://anaconda.org/bioconda/hiline)
# HiLine
HiC alignment and classification pipeline.

# Bioconda
HiLine is available on [bioconda](https://bioconda.github.io/)

```bash
conda install hiline
```

# As a Command Line Program
```bash
HiLine --help
```

# As a Python Class
```python
from HiLine import Pipeline
help(Pipeline)
```

# Requirments, Running
* Unix OS (MacOS / Linux)
* [Python](https://www.python.org/) version 3.8 or later
* At least one of
    * [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) version 2.0 or later
    * [bwa](https://github.com/lh3/bwa) version 0.7.17 or later
    * [minimap2](https://github.com/lh3/minimap2) version 2.17-r941 or later
        * [gawk](https://www.gnu.org/software/gawk/) version 5.1.0 or later (if using minimap2)
        * [perl](https://www.perl.org/) version 5.30.3 or later (if using minimap2)
* [Samtools](http://www.htslib.org/) version 1.10 or later

# Requirments, Installing
* [Python](https://www.python.org/) version 3.8 or later
* [Pip](https://pypi.org/project/pip/)
* [Setuptools](https://setuptools.readthedocs.io/en/latest/)
* [NumPy](https://numpy.org/) version 1.18.1 or later
* C/C++ compiler (tested with [Clang](https://clang.llvm.org/) version 9 and [GCC](https://gcc.gnu.org/) version 7.3)<br/>
* [Git](https://git-scm.com/)
* [Ninja](https://ninja-build.org/)
* [Hyperscan](https://github.com/intel/hyperscan) build requirements
    * [CMake](https://cmake.org/)
    * [Boost](https://www.boost.org/) version 1.61 or later
    * [Ragel](http://www.colm.net/open-source/ragel/)
    * [PCRE](https://www.pcre.org/)
    * [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/)
* [pytest](https://docs.pytest.org/en/latest/) (for testing only)

```bash
cd HiLine
pip install .
pytest
```

# Notes
* bwa-mem2 is the default mapper used by HiLine and is recommended over bwa, sans technical restrictions e.g. memory usage. 
* By default, HiLine performs a 'read-trimming' step as part of the alignment process. A restriction enzyme-digested version of the reference genome is created, with filled-in blunt (biotin) ends. Reads are first mapped to the digested reference, portions of reads that align past the end of a restriction fragment (not including the blunt ends) are trimmed away. Trimmed reads are then aligned to the original reference. Trimming results in more accurate alignments, but halves the overall speed of the alignment process. Turn off read-trimming with the '--no-trim' options.
* Minimap2 alignments are processed through a script from the [Arima Genomics mapping pipeline](https://github.com/ArimaGenomics/mapping_pipeline) to select for 5' mappings, which filters out all other chimeric supplementary alignments. Such alignments, therefore, will not be available for output. Minimap2 mode is only recommended for quick, approximate mappings.

# Third-Party Acknowledgements
HiLine uses the following third-party software and libraries<br/>
* Software
  * [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
  * [bwa](https://github.com/lh3/bwa)
  * [minimap2](https://github.com/lh3/minimap2)
  * [Samtools](http://www.htslib.org/)
* Libraries
    * C/C++
      * [Hyperscan](https://github.com/intel/hyperscan)
      * [stb_sprintf.h](https://github.com/nothings/stb/blob/master/stb_sprintf.h)
    * Python
      * [Click](https://palletsprojects.com/p/click/)
      * [Pandas](https://pandas.pydata.org/)
      * [NumPy](https://numpy.org/)
      * [Seaborn](https://seaborn.pydata.org/)
      * [Matplotlib](https://matplotlib.org/)
      * [Biopython](https://biopython.org/)
    * Perl
      * [filter_five_end.pl](https://github.com/ArimaGenomics/mapping_pipeline/blob/master/filter_five_end.pl)
