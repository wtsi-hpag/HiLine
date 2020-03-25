# HiLine
HiC alignment and classification pipeline.

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
* [BWA](https://github.com/lh3/bwa) version 0.7.17 or later
* [Samtools](http://www.htslib.org/) version 1.10 or later

# Requirments, Installing
* [Python](https://www.python.org/) version 3.8 or later
* [Pip](https://pypi.org/project/pip/)
* [Setuptools](https://setuptools.readthedocs.io/en/latest/)
* [NumPy](https://numpy.org/) version 1.18.1 or later
* C/C++ compiler (tested with [Clang](https://clang.llvm.org/) version 9 and [GCC](https://gcc.gnu.org/) version 7.3)<br/>
* [Git](https://git-scm.com/)
* [Ninja](https://ninja-build.org/)
* [Hyperscan](https://github.com/intel/hyperscan) build requirments
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

# Third-Party Acknowledgements
HiLine uses the following third-party software and libraries<br/>
* Software
  * [BWA](https://github.com/lh3/bwa)
  * [Samtools](http://www.htslib.org/)
* Libraries
    * C/C++
      * [Hyperscan](https://github.com/intel/hyperscan)
      * [stb_sprintf.h](https://github.com/nothings/stb/blob/master/stb_sprintf.h)
    * Python
      * [Click](https://palletsprojects.com/p/click/)
      * [pandas](https://pandas.pydata.org/)
      * [NumPy](https://numpy.org/)
      * [Seaborn](https://seaborn.pydata.org/)
      * [Matplotlib](https://matplotlib.org/)
      * [Biopython](https://biopython.org/)
