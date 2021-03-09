# Copyright (c) 2020 Ed Harry, Wellcome Sanger Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import logging
import os
import re
import sys
import warnings
from binascii import hexlify
from functools import update_wrapper
from os import makedirs
from os.path import isfile, join, basename, normpath
from subprocess import Popen, PIPE, CalledProcessError, check_output, STDOUT
from sys import platform
from threading import Thread

import click
import numpy as np
import pandas as pd
from Bio.Restriction import AllEnzymes, RestrictionBatch
from _HiLine import _HiLine_Main

NAME = "HiLine"
VERSION = "0.2.1"
DESCRIPTION = "A HiC alignment and classification pipeline."
LICENCE = "Copyright (c) 2020 Ed Harry, Wellcome Sanger Institute."

logger = logging.getLogger(NAME)
my_pid = os.getpid()


def create_logger_handle(stream, typeid, level):
    class LogFilter(logging.Filter):
        def __init__(self, level):
            super().__init__()
            self.__level = level

        def filter(self, record):
            return record.levelno == self.__level

    handle = logging.StreamHandler(stream=stream)
    handle.setLevel(level=level)
    handle.setFormatter(
        logging.Formatter("[%(name)s {id}] :: %(message)s".format(id=typeid))
    )
    handle.addFilter(LogFilter(level=level))

    return handle


class Pipeline(object):
    """
    HiC alignment, classification and analysis pipeline runner

    Performs HiC read alignment (using bwa mem, bwa-mem2 or minimap2) / reads in an existing alignment -> performs optional duplicate marking (using samtools markdups) -> performs HiC read classification based on an in-silico reference digest -> optionally outputs reads to different files based on classification -> returns read statistics


    Usage example
    -------------
    pipeline = Pipeline()

    pipeline.logger = <python logger>

    pipeline.reference = "<path to reference FASTA>"
    pipeline.restriction_sites = "<enzyme set>"
    pipeline.threads = <number of threads to use>
    pipeline.min_mapq = <minimum mapping quality>

    pipeline.register_alignment_sam_reads(reads="<path to reads in SAM/BAM/CRAM format>", mark_dups=<True / False>)
    pipeline.register_output_file(file_name="<path to output SAM/BAM/CRAM>", sort=<True / False>, handle=pipeline.output.<handle>)
    pipeline.save_stats_path = "<path to output stats folder>"

    pipeline.run()


    Initialisation
    --------------
    pipeline = Pipeline()


    Required attributes (all must be set prior to running the pipeline)
    -------------------------------------------------------------------
    pipeline.reference (str): Path to a reference genome file in FASTA format. May be compressed with gzip.
    pipeline.restriction_sites (str): /^<restriction>(, <restriction>)*$/ A comma-separated list of restriction site definitions. See 'restriction_sites' documentation for accepted site formats.
    pipeline.threads (int): Number of threads to use. Must be at least 3.
    pipeline.min_mapq (int): Minimum mapping quality. Reads below this threshold are classified as such and are considered as bad reads.
    pipeline.sam_input (Pipeline.SamReader / Pipeline.Align): SAM formatted input. One of the two indicated sub-classes, see their respective documentation. Also see the 'Registering SAM input' section below.


    Optional attributes
    -------------------
    pipeline.logger (logging.Logger): Message logger. Defaults to None i.e. no messages logged.
    pipeline.save_stats_path (str): Save path for statistics. Defaults to None i.e. statistics not saved to disk.
    pipeline.exception_callback (collections.abc.Callable): Callback when an exception is encountered, must have the signature: callback(Exception). Defaults to 'raise exception'.


    Registering SAM input
    ---------------------
    Helper functions that fill the 'sam_input' attribute.

    pipeline.register_sam_input: Sets a SAM/BAM/CRAM formatted alignment source (can be <stdin>). Used to process an external alignment.
    pipeline.register_alignment: Sets a read source in FASTQ format for alignment as part of the pipeline.
    pipeline.register_alignment_sam_reads: Sets a read source in SAM/BAM/CRAM format for alignment as part of the pipeline.
    pipeline.register_alignment_index_only: Just creates and saves a reference index rather than performing an alignment.

    See each function's documentation for details.


    Registering SAM output
    ----------------------
    pipeline.register_output_file: Function for setting output files. See function documentation for details.


    Running the pipeline
    --------------------
    pipeline.run: Function that runs the pipeline. Returns a Pipeline.Stats object. See function documentation for details.


    Read classification
    -------------------
    Reads are classified into 18 types:
        Good reads:
            Valid pairs:
                1. FF
                2. FR
                3. RF
                4. RR
                F/R refer to the read direction relative to the respective reference sequence of the 1st/2nd read respectively. If both reads in a pair map to the same reference sequence then the 1st read is defined as the read with the lowest mapped 5' coordinate. Otherwise for read-pairs mapped to different reference sequences read 1/2 assignment is identical to the SAM assignment.

            Invalid pairs:
                5. Self-circle (RF)
                6. Dangling-end (FR)
                7. Same-fragment-and-strand (FF or RR)
                The three types above are read-pairs that map to the same restriction fragment.
                8. Re-ligation
                Re-ligated pairs are read pairs that map to adjacent restriction fragments.

            Invalid reads:
                9. Dumped
                Dumped reads are reads classified as good but whose mate read is classified as bad, hence cannot form a read-pair.

        Bad reads:
            10. Below-min-mapq
            Reads with a mapping quality below the minimum threshold.

            11. Too-far-from-restriction-site
            Reads whose 5' mapping coordinate is too far from a restriction site, defined to be twice the read length.

            12. Invalid-reference
            Reads mapping to a sequence with an ID not found in the defined reference.

            13. Unmapped
            Reads with the 0x4 SAM flag set.

            14. Unpaired
            Reads with the 0x1 SAM flag unset.

            15. Secondary
            Reads with the 0x100 SAM flag set.

            16. QC-fail
            Reads with the 0x200 SAM flag set.

            17. Duplicate
            Reads with the 0x400 SAM flag set.

            18. Supplementary
            Reads with the 0x800 SAM flag set.
    """

    class OutputHandles(object):
        """Output handles for each read-type"""

        class Handle(object):
            """Output handle, wrapper around a file descriptor"""

            def __init__(self):
                self._fd = 0

            @property
            def fd(self):
                return self._fd

            @fd.setter
            def fd(self, fd):
                if isinstance(fd, int):
                    self._fd = fd
                elif hasattr(fd, "fileno"):
                    self._fd = fd.fileno()
                else:
                    raise TypeError(
                        "{fd} does not have a file descriptor".format(fd=fd)
                    )

            def fileno(self):
                return self.fd

        def __init__(self):
            self._unPaired = self.Handle()
            self._unMapped = self.Handle()
            self._secondary = self.Handle()
            self._qcFail = self.Handle()
            self._duplicate = self.Handle()
            self._supplementary = self.Handle()
            self._invalidReferenceName = self.Handle()
            self._belowMinMapq = self.Handle()
            self._tooFarFromRestrictionSite = self.Handle()
            self._dump = self.Handle()
            self._selfCircle = self.Handle()
            self._danglingEnd = self.Handle()
            self._sameFragmentAndStrand = self.Handle()
            self._reLigated = self.Handle()
            self._FF = self.Handle()
            self._FR = self.Handle()
            self._RF = self.Handle()
            self._RR = self.Handle()

        @property
        def unpaired(self):
            """Unpaired reads"""
            return self._unPaired

        @unpaired.setter
        def unpaired(self, fd):
            self._unPaired.fd = fd

        @property
        def unmapped(self):
            """Unmapped reads"""
            return self._unMapped

        @unmapped.setter
        def unmapped(self, fd):
            self._unMapped.fd = fd

        @property
        def secondary(self):
            """Secondary reads"""
            return self._secondary

        @secondary.setter
        def secondary(self, fd):
            self._secondary.fd = fd

        @property
        def qcfail(self):
            """QC-failed reads"""
            return self._qcFail

        @qcfail.setter
        def qcfail(self, fd):
            self._qcFail.fd = fd

        @property
        def duplicate(self):
            """Duplicate reads"""
            return self._duplicate

        @duplicate.setter
        def duplicate(self, fd):
            self._duplicate.fd = fd

        @property
        def supplementary(self):
            """Supplementary reads"""
            return self._supplementary

        @supplementary.setter
        def supplementary(self, fd):
            self._supplementary.fd = fd

        @property
        def invalidreferencename(self):
            """Reads with an invalid reference"""
            return self._invalidReferenceName

        @invalidreferencename.setter
        def invalidreferencename(self, fd):
            self._invalidReferenceName.fd = fd

        @property
        def belowminmapq(self):
            """Reads below minimum mapping quality"""
            return self._belowMinMapq

        @belowminmapq.setter
        def belowminmapq(self, fd):
            self._belowMinMapq.fd = fd

        @property
        def toofarfromrestrictionsite(self):
            """Reads too far from a restriction site"""
            return self._tooFarFromRestrictionSite

        @toofarfromrestrictionsite.setter
        def toofarfromrestrictionsite(self, fd):
            self._tooFarFromRestrictionSite.fd = fd

        @property
        def dump(self):
            """Dumped reads (valid_pairs reads with an invalid partner)"""
            return self._dump

        @dump.setter
        def dump(self, fd):
            self._dump.fd = fd

        @property
        def selfcircle(self):
            """Self-circular read pairs"""
            return self._selfCircle

        @selfcircle.setter
        def selfcircle(self, fd):
            self._selfCircle.fd = fd

        @property
        def danglingend(self):
            """Dangling-end read pairs"""
            return self._danglingEnd

        @danglingend.setter
        def danglingend(self, fd):
            self._danglingEnd.fd = fd

        @property
        def samefragmentandstrand(self):
            """Read pairs on the same restriction fragment and strand"""
            return self._sameFragmentAndStrand

        @samefragmentandstrand.setter
        def samefragmentandstrand(self, fd):
            self._sameFragmentAndStrand.fd = fd

        @property
        def religated(self):
            """Re-ligated read pairs"""
            return self._reLigated

        @religated.setter
        def religated(self, fd):
            self._reLigated.fd = fd

        @property
        def ff(self):
            """Valid ff read pairs"""
            return self._FF

        @ff.setter
        def ff(self, fd):
            self._FF.fd = fd

        @property
        def fr(self):
            """Valid fr read pairs"""
            return self._FR

        @fr.setter
        def fr(self, fd):
            self._FR.fd = fd

        @property
        def rf(self):
            """Valid rf read pairs"""
            return self._RF

        @rf.setter
        def rf(self, fd):
            self._RF.fd = fd

        @property
        def rr(self):
            """Valid rr read pairs"""
            return self._RR

        @rr.setter
        def rr(self, fd):
            self._RR.fd = fd

        @property
        def descriptors(self):
            return {
                self.unpaired: "unpaired reads",
                self.unmapped: "unmapped reads",
                self.secondary: "secondary reads",
                self.qcfail: "qc-failed reads",
                self.duplicate: "duplicate reads",
                self.supplementary: "supplementary reads",
                self.invalidreferencename: "reads with an invalid reference",
                self.belowminmapq: "reads below minimum mapping quality",
                self.toofarfromrestrictionsite: "reads too far from a restriction site",
                self.dump: "dumped reads (good reads with a bad mate read)",
                self.selfcircle: "self-circular read pairs",
                self.danglingend: "dangling-end read pairs",
                self.samefragmentandstrand: "read pairs on the same restriction fragment and strand",
                self.religated: "re-ligated read pairs",
                self.ff: "valid ff read pairs",
                self.fr: "valid fr read pairs",
                self.rf: "valid rf read pairs",
                self.rr: "valid rr read pairs",
            }

        @property
        def valid_pairs(self):
            """All valid_pairs pair-types"""
            return self.ff, self.fr, self.rf, self.rr

        @property
        def invalid_pairs(self):
            """All invalid pair-types"""
            return (
                self.selfcircle,
                self.danglingend,
                self.samefragmentandstrand,
                self.religated,
            )

        @property
        def invalid_reads(self):
            """All invalid read-types"""
            return (
                self.selfcircle,
                self.danglingend,
                self.samefragmentandstrand,
                self.religated,
                self.dump,
            )

        @property
        def good_reads(self):
            """All good reads"""
            return (
                self.ff,
                self.fr,
                self.rf,
                self.rr,
                self.selfcircle,
                self.danglingend,
                self.samefragmentandstrand,
                self.religated,
                self.dump,
            )

        @property
        def bad_reads(self):
            """All bad read-types"""
            return (
                self.unpaired,
                self.unmapped,
                self.secondary,
                self.qcfail,
                self.duplicate,
                self.supplementary,
                self.invalidreferencename,
                self.belowminmapq,
                self.toofarfromrestrictionsite,
            )

        @property
        def all_reads(self):
            """All read types"""
            return (
                self.ff,
                self.fr,
                self.rf,
                self.rr,
                self.selfcircle,
                self.danglingend,
                self.samefragmentandstrand,
                self.religated,
                self.dump,
                self.unpaired,
                self.unmapped,
                self.secondary,
                self.qcfail,
                self.duplicate,
                self.supplementary,
                self.invalidreferencename,
                self.belowminmapq,
                self.toofarfromrestrictionsite,
            )

        def description(self, handle):
            return self.descriptors[handle]

        class Exception(Exception):
            pass

        def validate(self):
            for fd, description in self.descriptors.items():
                if not (
                    (isinstance(fd, int) and fd >= 0)
                    or (hasattr(fd, "fileno") and fd.fileno() >= 0)
                ):
                    raise self.Exception(
                        "invalid file descriptor: {fd}, '{des}'".format(
                            fd=fd, des=description
                        )
                    )

    def __init__(self):
        self._reference = None
        self._restrictionSites = None
        self._threads = None
        self._minMapq = None
        self._samInput = None

        self._samInputWrapper = None
        self._output = self.OutputHandles()
        self._outputFiles = {}
        self._outputCommands = []

        self._fastaHandle = None

        self._logger = self.LoggerWrapper()
        self._loggerHandle = None

        self._runParms = None

        def _exception_callback_default(exception):
            raise exception

        self._exceptionCallback = _exception_callback_default

        self._saveStatsPath = None

        self._commandLine = None

    class Exception(Exception):
        pass

    @property
    def output(self):
        return self._output

    @property
    def output_files(self):
        return self._outputFiles

    @property
    def output_commands(self):
        return self._outputCommands

    @property
    def logger(self):
        """logging.Logger, message logger for pipeline processes"""
        return self._logger

    @logger.setter
    def logger(self, logger):
        self._logger.logger = logger

    @property
    def save_stats_path(self):
        return self._saveStatsPath

    @save_stats_path.setter
    def save_stats_path(self, path):
        if self._saveStatsPath is not None and self._saveStatsPath != path:
            self.log_error(
                exception=self.Exception(
                    "Save stats path set twice: '{orig}' and '{new}'".format(
                        orig=self._saveStatsPath, new=path
                    )
                )
            )
        self._saveStatsPath = path

    @property
    def exception_callback(self):
        return self._exceptionCallback

    @exception_callback.setter
    def exception_callback(self, callback):
        if not callable(callback):
            self.log_error(
                exception=self.Exception(
                    "Exception callback '{callback}' is not callable".format(
                        callback=callback
                    )
                )
            )
        self._exceptionCallback = callback

    @property
    def command_line(self):
        return (
            self._commandLine
            if self._commandLine is not None
            else "Called as python library, num_threads: {threads}, min_mapq: {mapq}, reference: {ref}, restriction_sites: {sites}, {sam_in}, {sam_out}{save_stats}".format(
                threads=self.threads,
                mapq=self.min_mapq,
                ref=self.reference,
                sites=self.restriction_sites,
                sam_in=self.sam_input,
                sam_out=", ".join(
                    "{handle}: {out_file}".format(
                        handle=self.output.description(handle=handle), out_file=out_file
                    )
                    for out_file, handles in self.output_files.items()
                    for handle in handles
                ),
                save_stats=""
                if self.save_stats_path is None
                else ", save_stats_path: {path}".format(path=self.save_stats_path),
            )
        )

    @command_line.setter
    def command_line(self, cl):
        self._commandLine = cl

    @property
    def reference(self):
        """Path to reference genome in (gzipped) FASTA format."""
        return self._reference

    @reference.setter
    def reference(self, reference):
        if self._reference is not None and self._reference != reference:
            self.log_error(
                exception=self.Exception(
                    "Reference reference set twice: '{orig}' and '{new}'".format(
                        orig=self._reference, new=reference
                    )
                )
            )
        self._reference = reference

    @property
    def restriction_sites(self):
        """
        Restriction site(s) specification string, a comma-separated list of restriction site definitions.
        /^<restriction>(, <restriction>)*$/

        where <restriction> is either: the name of a HiC kit provider, the name of a restriction enzyme, or the string-definition of the restriction site.
        <restriction> = /(<kit>|<enzyme>|<site>)/

        <kit> is the name of a HiC kit provider, currently supported are:
            Arima Genomics version 2: "Arima_v2"
            Arima Genomics: "Arima"
            Dovetail Genomics: "Dovetail"
            Phase Genomics: "Phase"
            Qiagen: "Qiagen"
        <kit> = /((Arima_v2)|(Arima)|(Dovetail)|(Phase)|(Qiagen))/, note: names are case-insensitive.

        <enzyme> is the name of a restriction enzyme e.g. "DpnII", currently all known enzymes defined in the 'Biopython.Restriction' package are supported.
        Note: enzyme names are case-sensitive.

        <site> is the IUPAC string definition of a restriction site e.g. "^GATC".
        <site> = /([ACGTWSMKRYBDHVN]*\^[ACGTWSMKRYBDHVN]*_?[ACGTWSMKRYBDHVN]*)/
        A caret character "^" is required to appear once and only once, and defines the cut location.
        If the site is not palindromic, an underscore character "_" is required to appear once and only once, and defines the cut location on the reverse strand. Otherwise it is not required.
        """
        return self._restrictionSites

    @restriction_sites.setter
    def restriction_sites(self, restriction_sites):
        if (
            self._restrictionSites is not None
            and self._restrictionSites != restriction_sites
        ):
            self.log_error(
                exception=self.Exception(
                    "Restriction sites set twice: '{orig}' and '{new}'".format(
                        orig=self._restrictionSites, new=restriction_sites
                    )
                )
            )
        self._restrictionSites = restriction_sites

    @property
    def threads(self):
        return self._threads

    @threads.setter
    def threads(self, threads):
        if self._threads is not None and self._threads != threads:
            self.log_error(
                exception=self.Exception(
                    "Number of threads set twice: '{orig}' and '{new}'".format(
                        orig=self._threads, new=threads
                    )
                )
            )
        self._threads = threads

    @property
    def min_mapq(self):
        return self._minMapq

    @min_mapq.setter
    def min_mapq(self, min_mapq):
        if self._minMapq is not None and self._minMapq != min_mapq:
            self.log_error(
                exception=self.Exception(
                    "Minimum mapping quality set twice: '{orig}' and '{new}'".format(
                        orig=self._minMapq, new=min_mapq
                    )
                )
            )
        self._minMapq = min_mapq

    @property
    def sam_input(self):
        """Either a Pipeline.SamReader or a Pipeline.Align object"""
        return self._samInput

    @sam_input.setter
    def sam_input(self, sam):
        if self._samInput is not None and self._samInput != sam:
            self.log_error(
                exception=self.Exception("More than one SAM input source set")
            )
        self._samInput = sam

    class FastaHandle(object):
        def __init__(self, path):
            class ZipCat(object):
                def __init__(self, path):
                    self.path = path

                @property
                def cat(self):
                    if self.gz(path=self.path):
                        cat = "zcat" if "linux" in platform else "gzcat"
                    elif self.bz2(path=self.path):
                        cat = "bzcat"
                    else:
                        cat = "cat"

                    return "{cat} {file}".format(cat=cat, file=self.path)

                @staticmethod
                def test(path, magic):
                    with open(path, "rb") as file:
                        return hexlify(file.read(len(magic) // 2)) == magic

                @classmethod
                def gz(cls, path):
                    return cls.test(path=path, magic=b"1f8b")

                @classmethod
                def bz2(cls, path):
                    return cls.test(path=path, magic=b"425a68")

            self.proc = Popen(ZipCat(path).cat.split(), stdout=PIPE)
            self.name = path

        def __enter__(self):
            self.proc.__enter__()
            return Pipeline.HandleWrapper(handle=self.proc.stdout, name=self.name)

        def __exit__(self, type, value, traceback):
            self.proc.__exit__(type, value, traceback)

    class HandleWrapper(object):
        def __init__(self, handle, name):
            self.handle = handle
            self.name = name

        def fileno(self):
            return self.handle.fileno()

    class SamPGChain(object):
        class Exception(Exception):
            pass

        class PGLine(object):
            def __init__(self, line):
                self._lineDict = {}
                for it in line.decode("utf-8", errors="ignore")[:-1].split("\t")[1:]:
                    self._lineDict.setdefault(it[:2], []).append(it[3:])

            @property
            def id(self):
                return self._lineDict["ID"][0]

            @id.setter
            def id(self, id):
                self._lineDict["ID"][0] = id

            @property
            def pp(self):
                return self._lineDict.setdefault("PP", [None])[0]

            @pp.setter
            def pp(self, id):
                self._lineDict.setdefault("PP", [None])[0] = id

            def __str__(self):
                dic = self._lineDict.copy()
                line = ["@PG", "ID:" + dic.pop("ID")[0]]
                if "PN" in dic:
                    line += ["PN:" + dic.pop("PN")[0]]
                if "PP" in dic:
                    if dic["PP"][0] is not None:
                        line += ["PP:" + dic.pop("PP")[0]]
                    else:
                        del dic["PP"]
                for key, item in dic.items():
                    for it in item:
                        line += ["{id}:{val}".format(id=key, val=it)]
                return "\t".join(line) + "\n"

        def __init__(self, lines=[]):
            lines = [self.PGLine(line) for line in lines]
            pps = [line.pp for line in lines]

            end_idx = []
            for i, line in enumerate(lines[::-1]):
                if line.id not in pps:
                    end_idx.append(len(lines) - i - 1)

            if len(end_idx) > 0 or len(lines) == 0:
                for end, idx in enumerate(end_idx):
                    lines[idx], lines[-1 - end] = lines[-1 - end], lines[idx]

                self._chain = lines
                self._num_chains = max(len(end_idx), 1)
            else:
                raise self.Exception("Invalid SAM PG chain detected")

        @property
        def ids(self):
            return [line.id for line in self._chain]

        @property
        def pps(self):
            return [line.pp for line in self._chain]

        @property
        def num_chains(self):
            return self._num_chains

        def append(self, line_in):
            for _ in range(self._num_chains):

                line = self.PGLine(line_in)

                if len(self._chain) > 0:
                    id = line.id
                    ids = self.ids
                    if id in ids:
                        match = re.match(r"^(?P<name>.+)\.(?P<num>\d+)", id)
                        if match:
                            id = match.group("name")
                            n = int(match.group("num")) + 1
                        else:
                            n = 1
                        test = ids[0]
                        while test in ids:
                            test = id + ".{n}".format(n=n)
                            n += 1
                        line.id = test

                    line.pp = self._chain[-self._num_chains].id

                self._chain.append(line)

        def __iter__(self):
            for line in self._chain:
                yield str(line).encode("utf-8")

    def __enter__(self):
        try:
            try:
                match = re.match(
                    r"^samtools (?P<major>\d+)\.(?P<minor>\d+)",
                    check_output("samtools --version".split(), stderr=STDOUT).decode(
                        "utf-8"
                    ),
                )
                if match is None:
                    raise self.Exception("Could not determine 'samtools' version")

                major = int(match.group("major"))
                minor = int(match.group("minor"))

                if not (major > 1 or (major == 1 and minor >= 10)):
                    raise self.Exception(
                        "'samtools' version {major}.{minor} found, version 1.10 or greater required".format(
                            major=major, minor=minor
                        )
                    )

                samtools_version = "{major}.{minor}".format(major=major, minor=minor)

            except FileNotFoundError:
                raise self.Exception("'samtools' not found on $PATH")
            except CalledProcessError as ex:
                raise self.Exception(str(ex))

            if self.reference is None:
                raise self.Exception("Reference reference file has not been set")
            if self.restriction_sites is None:
                raise self.Exception("Restriction sites have not been set")
            if self.threads is None:
                raise self.Exception("Number of threads has not been set")
            if self.min_mapq is None:
                raise self.Exception("Minimum mapping quality has not been set")
            if self.sam_input is None:
                raise self.Exception("SAM input has not been set")

            restriction_sites = self.RestrictionSites(logger=self.logger).create_sites(
                self.restriction_sites
            )

            self._loggerHandle = self.LoggerHandle(logger=self.logger)
            if not isinstance(self.sam_input, self.IndexOnly):
                for file, handles in self.output_files.items():
                    logger.info(
                        "Writing\n\t{desc}\n\t\tto file '{file}'{extra}".format(
                            file=file,
                            desc="\n\t".join(
                                [
                                    "'{des}'".format(
                                        des=self.output.description(handle)
                                    )
                                    for handle in handles
                                ]
                            ),
                            extra=" (sorted and indexed)" if file.sort else "",
                        )
                    )

                    logger_handle = self._loggerHandle.add_logger(
                        log_func=self.logger.info, id="writing {file}".format(file=file)
                    )
                    command1 = Popen(
                        (
                            "samtools sort -@ {threads} -l 0 -".format(
                                threads=self.threads
                            )
                            if file.sort
                            else "samtools view -b --output-fmt-option level=0 -"
                        ).split(),
                        stdin=PIPE,
                        stdout=PIPE,
                        stderr=logger_handle,
                    )
                    command2 = Popen(
                        "samtools view -@ {threads} -h -o {file} --output-fmt-option level={level}{extra}{index} -".format(
                            file=file if file.name != "<stdout>" else "-",
                            level="9" if file.name != "<stdout>" else "0",
                            extra=" --reference {ref} --output-fmt-option use_bzip2=1 --output-fmt-option use_lzma=1".format(
                                ref=self.reference
                            )
                            if file.name.endswith(".cram")
                            else "",
                            threads=self.threads,
                            index=" --write-index"
                            if file.sort and file.name != "<stdout>"
                            else "",
                        ).split(),
                        stdin=command1.stdout,
                        stderr=logger_handle,
                    )
                    command1.__enter__()
                    command2.__enter__()
                    self.output_commands.append(command1)
                    for handle in handles:
                        handle.fd = command1.stdin

                self.output.validate()

            class HandleFileGuard(object):
                def __init__(self, handle, name, file=None):
                    self.wrapper = Pipeline.HandleWrapper(handle=handle, name=name)
                    self.file = file

                def __enter__(self):
                    return self.wrapper

                def __exit__(self, exc_type, exc_val, exc_tb):
                    if self.file is not None:
                        self.file.close()

            self._pipeThreads = []
            input_file_handle = None
            num_programs = 3

            if isinstance(self.sam_input, self.SamReader):
                num_programs += 1

                name = self.sam_input.name
                if name == "<stdin>":
                    input = sys.stdin
                else:
                    input = open(name, "rb")

                input_file_handle = input
                input = Popen(
                    "samtools collate -@ {threads} -Ou - .{pid}.in.read.collate".format(
                        threads=self.threads, pid=my_pid
                    ).split(),
                    stdin=input,
                    stdout=PIPE,
                    stderr=self._loggerHandle.add_logger(
                        log_func=self.logger.info,
                        id="collating {file}".format(file=name),
                    ),
                ).stdout

            elif isinstance(self.sam_input, self.Align):
                try:
                    with Popen(
                        "_HiLine_Aligner --help".split(), stdout=PIPE, stderr=STDOUT
                    ) as process:
                        output = "".join(
                            [
                                line.decode("utf-8")
                                for line in process.stdout.readlines()
                            ]
                        )

                    match = re.search(
                        r"_HiLine_Aligner (?P<version>(\d\.)+\d)\.",
                        output,
                    )

                    if match is None:
                        raise self.Exception(
                            "Could not determine '_HiLine_Aligner' version"
                        )

                    if match.group("version") != VERSION:
                        raise self.Exception(
                            "_HiLine_Aligner version {al_ver} != HiLine version {ver}".format(
                                al_ver=match.group("version"), ver=VERSION
                            )
                        )

                except FileNotFoundError:
                    raise self.Exception("'_HiLine_Aligner' not found on $PATH")
                except CalledProcessError as ex:
                    raise self.Exception(str(ex))

                name = "Alignment"
                aligner_cmd = "_HiLine_Aligner {version} {trim} {sites} {tags} -t {threads} {index} {ref}".format(
                    sites=" ".join(
                        "--site {s} {f} {r}".format(
                            s=site.pattern, f=site.loc, r=site.rev_loc
                        )
                        for site in restriction_sites
                    ),
                    tags=" ".join(
                        "--tag {t}".format(t=t)
                        for t in (
                            self.sam_input.tags
                            if isinstance(self.sam_input, self.SamReads)
                            else []
                        )
                    ),
                    threads=self.threads,
                    ref=self.reference,
                    version="--bwa1"
                    if self.sam_input.version == 1
                    else ("--bwa2" if self.sam_input.version == 2 else "--minimap2"),
                    trim="--trim" if self.sam_input.trim else "--no-trim",
                    index="--no-align"
                    if isinstance(self.sam_input, self.IndexOnly)
                    else "",
                )
                num_align_prog = 5 if self.sam_input.trim else 1

                if isinstance(self.sam_input, self.IndexOnly):
                    with Popen(
                        (aligner_cmd + " -").split(),
                        stderr=self._loggerHandle.add_logger(
                            log_func=self.logger.info, id="Index Reference"
                        ),
                    ) as process:
                        process.communicate()
                        self._samInputWrapper = None
                        self._runParms = self.IndexRun()
                        return self._runParms

                if isinstance(self.sam_input, self.SamReads):
                    num_programs += 2 + num_align_prog

                    file_name = self.sam_input.reads1
                    if file_name == "<stdin>":
                        input = sys.stdin
                    else:
                        input = open(file_name, "rb")

                    input_file_handle = input
                    command_log_handle = self._loggerHandle.add_logger(
                        log_func=self.logger.info,
                        id="{file} -> fastq".format(file=file_name),
                    )
                    input = Popen(
                        "samtools collate -@ {threads} -Ouf - .{pid}.reads.collate".format(
                            threads=self.threads, pid=my_pid
                        ).split(),
                        stdin=input,
                        stdout=PIPE,
                        stderr=command_log_handle,
                    ).stdout

                    read, write = os.pipe()
                    sam_pg_rg_capture_input = input
                    pg_lines = []
                    rg_lines = []
                    global headermode
                    headermode = True
                    pg_rg_capture_write = write

                    def sam_pg_rg_capture():
                        global headermode
                        seen_header = False
                        with os.fdopen(pg_rg_capture_write, "wb") as f_out:
                            for line in sam_pg_rg_capture_input:
                                if headermode:
                                    if (
                                        decoded_line := line.decode(
                                            "utf-8", errors="backslashreplace"
                                        )
                                    ).startswith("@"):
                                        seen_header = True
                                        if decoded_line.startswith("@PG"):
                                            pg_lines.append(line)
                                        elif decoded_line.startswith("@RG"):
                                            rg_lines.append(line)
                                    elif seen_header:
                                        headermode = False

                                f_out.write(line)
                            headermode = False

                    thread = Thread(target=sam_pg_rg_capture)
                    thread.start()
                    self._pipeThreads.append(thread)
                    input = read

                    samtools_fastq_cmd = "samtools fastq -t{tags} -@ {threads} -0 /dev/null -F 0x900 -".format(
                        threads=self.threads,
                        tags=""
                        if len(self.sam_input.tags) == 0
                        else "T {tg}".format(
                            tg=",".join(t for t in self.sam_input.tags)
                        ),
                    )
                    input = Popen(
                        samtools_fastq_cmd.split(),
                        stdin=input,
                        stdout=PIPE,
                        stderr=command_log_handle,
                    ).stdout
                    input = Popen(
                        (aligner_cmd + " -").split(),
                        stdin=input,
                        stdout=PIPE,
                        stderr=self._loggerHandle.add_logger(
                            log_func=self.logger.info, id="Align"
                        ),
                    ).stdout

                    read, write = os.pipe()
                    sam_pg_rg_append_input = input
                    pg_rg_append_write = write

                    def sam_pg_rg_append():
                        global headermode
                        local_header_mode = True
                        seen_header = False
                        header_buffer = []
                        local_pg_lines = []
                        with os.fdopen(pg_rg_append_write, "wb") as f_out:
                            for line in sam_pg_rg_append_input:
                                if local_header_mode:
                                    if (
                                        decoded_line := line.decode(
                                            "utf-8", errors="backslashreplace"
                                        )
                                    ).startswith("@"):
                                        seen_header = True
                                        if decoded_line.startswith("@PG"):
                                            local_pg_lines.append(line)
                                        else:
                                            header_buffer.append(line)
                                    elif seen_header:
                                        while headermode:
                                            pass
                                        local_header_mode = False

                                        for headerLine in header_buffer:
                                            f_out.write(headerLine)

                                        for rgLine in rg_lines:
                                            f_out.write(rgLine)

                                        chain = self.SamPGChain(pg_lines)
                                        chain.append(
                                            "@PG\tID:samtools\tPN:samtools\tVN:{version}\tCL:{cmd}\n".format(
                                                version=samtools_version,
                                                cmd=samtools_fastq_cmd,
                                            ).encode(
                                                "utf-8"
                                            )
                                        )

                                        for pgLine in local_pg_lines:
                                            chain.append(pgLine)

                                        for pgLine in chain:
                                            f_out.write(pgLine)

                                        f_out.write(line)
                                    else:
                                        header_buffer.append(line)
                                else:
                                    f_out.write(line)

                    thread = Thread(target=sam_pg_rg_append)
                    thread.start()
                    self._pipeThreads.append(thread)
                    input = read

                else:
                    num_programs += num_align_prog

                    if self.sam_input.reads2 is None:
                        if self.sam_input.reads1 == "<stdin>":
                            input_file_handle = sys.stdin
                            input = Popen(
                                (aligner_cmd + " -").split(),
                                stdin=input_file_handle,
                                stdout=PIPE,
                                stderr=self._loggerHandle.add_logger(
                                    log_func=self.logger.info, id="Align"
                                ),
                            ).stdout
                        else:
                            input = Popen(
                                "{cmd} {reads}".format(
                                    cmd=aligner_cmd,
                                    reads=self.sam_input.reads1,
                                ).split(),
                                stdout=PIPE,
                                stderr=self._loggerHandle.add_logger(
                                    log_func=self.logger.info, id="Align"
                                ),
                            ).stdout
                    else:
                        if self.sam_input.reads1 == "<stdin>":
                            input_file_handle = sys.stdin
                            self.sam_input.reads1 = "-"
                        if self.sam_input.reads2 == "<stdin>":
                            input_file_handle = sys.stdin
                            self.sam_input.reads2 = "-"

                        if input_file_handle is None:
                            input = Popen(
                                "{cmd} {reads1} {reads2}".format(
                                    cmd=aligner_cmd,
                                    reads1=self.sam_input.reads1,
                                    reads2=self.sam_input.reads2,
                                ).split(),
                                stdout=PIPE,
                                stderr=self._loggerHandle.add_logger(
                                    log_func=self.logger.info, id="Align"
                                ),
                            ).stdout
                        else:
                            input = Popen(
                                "{cmd} {reads1} {reads2}".format(
                                    cmd=aligner_cmd,
                                    reads1=self.sam_input.reads1,
                                    reads2=self.sam_input.reads2,
                                ).split(),
                                stdin=input_file_handle,
                                stdout=PIPE,
                                stderr=self._loggerHandle.add_logger(
                                    log_func=self.logger.info, id="Align"
                                ),
                            ).stdout

            else:
                self.log_error(exception=self.Exception("unknown input type"))

            if self.sam_input.markDups:
                num_programs += 5

                command_log_handle = self._loggerHandle.add_logger(
                    log_func=self.logger.info,
                    id="marking duplicates in {file}".format(file=name),
                )

                command = Popen(
                    "samtools fixmate -@ {threads} -m --output-fmt-option level=0 - -".format(
                        threads=self.threads
                    ).split(),
                    stdin=input,
                    stdout=PIPE,
                    stderr=command_log_handle,
                )

                command = Popen(
                    "samtools sort -@ {threads} -l 0 -".format(
                        threads=self.threads
                    ).split(),
                    stdin=command.stdout,
                    stdout=PIPE,
                    stderr=command_log_handle,
                )

                command = Popen(
                    "samtools markdup -@ {threads} -S -d 100 -m s -T .{pid}.in.markdup --output-fmt-option level=0 - -".format(
                        threads=self.threads, pid=my_pid
                    ).split(),
                    stdin=command.stdout,
                    stdout=PIPE,
                    stderr=command_log_handle,
                )

                command = Popen(
                    "samtools collate -@ {threads} -Ou - .{pid}.in.markdup.collate".format(
                        threads=self.threads, pid=my_pid
                    ).split(),
                    stdin=command.stdout,
                    stdout=PIPE,
                    stderr=command_log_handle,
                )

                command = Popen(
                    "samtools view -h -".split(),
                    stdin=command.stdout,
                    stdout=PIPE,
                    stderr=command_log_handle,
                )

            else:
                num_programs += 1

                command = Popen(
                    "samtools view -h -".split(),
                    stdin=input,
                    stdout=PIPE,
                    stderr=self._loggerHandle.add_logger(
                        log_func=self.logger.info, id="reading {file}".format(file=name)
                    ),
                )

            command.__enter__()
            read, write = os.pipe()

            main_name = NAME + "_main"

            def add_pipeline_pg_line():
                local_header_mode = True
                seen_header = False
                header_buffer = []
                pg_lines = []
                with os.fdopen(write, "wb") as f_out:
                    for line in command.stdout:
                        if local_header_mode:
                            if (
                                decoded_line := line.decode(
                                    "utf-8", errors="backslashreplace"
                                )
                            ).startswith("@"):
                                seen_header = True
                                if decoded_line.startswith("@PG"):
                                    pg_lines.append(line)
                                else:
                                    header_buffer.append(line)
                            elif seen_header:
                                local_header_mode = False

                                for headerLine in header_buffer:
                                    f_out.write(headerLine)

                                pipeline_chain = []
                                for _ in range(num_programs - 3):
                                    chain = self.SamPGChain(pg_lines)
                                    pg_lines = [line for line in chain]
                                    pipeline_chain.insert(
                                        0, pg_lines[-chain.num_chains]
                                    )
                                    pg_lines = pg_lines[: -chain.num_chains]
                                chain = self.SamPGChain(pg_lines)

                                chain.append(
                                    "@PG\tID:{id}\tPN:{pn}\tCL:{cl}\tDS:{ds}\tVN:{vn}\n".format(
                                        id=NAME,
                                        pn=NAME,
                                        cl=self.command_line,
                                        vn=VERSION,
                                        ds="{des} The following {n} programs in this PG chain were called automatically as part of this pipeline.".format(
                                            des=DESCRIPTION, n=num_programs
                                        ),
                                    ).encode(
                                        "utf-8"
                                    )
                                )

                                for pgLine in pipeline_chain:
                                    chain.append(pgLine)

                                chain.append(
                                    "@PG\tID:{id}\tPN:{pn}\tVN:{vn}\n".format(
                                        id=main_name,
                                        pn=main_name + "_",
                                        vn=VERSION,
                                    ).encode("utf-8")
                                )

                                for pgLine in chain:
                                    f_out.write(pgLine)

                                f_out.write(line)
                            else:
                                header_buffer.append(line)
                        else:
                            f_out.write(line)

            thread = Thread(target=add_pipeline_pg_line)
            thread.start()
            self._pipeThreads.append(thread)

            self._samInputWrapper = HandleFileGuard(
                handle=os.fdopen(read, "rb"), name=name, file=input_file_handle
            )
            self._fastaHandle = self.FastaHandle(path=self.reference)
            lh = self._loggerHandle.__enter__()

            self._runParms = self.RunParams(
                name=NAME + "_main",
                version=VERSION,
                fasta_input=self._fastaHandle.__enter__(),
                sam_input=self._samInputWrapper.__enter__(),
                sam_output=self.output,
                restriction_sites=restriction_sites,
                n_threads=self.threads,
                min_map_q=self.min_mapq,
                info=lh.info,
                error=lh.error,
            )

        except (
            self.Exception,
            self.OutputHandles.Exception,
            self.RestrictionSites.Exception,
        ) as ex:
            self.log_error(exception=ex)
            self._runParms = None

        return self._runParms

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._runParms is not None:
            for command in self.output_commands + [
                self._fastaHandle,
                self._samInputWrapper,
                self._loggerHandle,
            ]:
                if command is not None:
                    command.__exit__(exc_type, exc_val, exc_tb)
            for thread in self._pipeThreads:
                thread.join()

    def register_output_file(self, file_name, handle, sort=False):
        """
        Assign a read-type to an output file.
        This function may be called multiple times with the same file_name parameter, assigning different read-types to the same file.
        However, the same read-type may not be assigned to multiple files.

        :param file_name: (str) Path to the SAM/BAM/CRAM output file, output type will be determined from the file extension. May be "<stdout>".
        :param handle: (Pipeline.OutputHandles.Handle) One of the output handles available in self.output.
        :param sort: (bool) Whether or not to sort and index the output file.


        Examples
        --------
        pipeline.register_output_file(file_name="self_circle.bam", handle=pipeline.output.selfcircle, sort=False)

        for handle in pipeline.output.valid_pairs:
            pipeline.register_output_file(file_name="results/valid_pairs.cram", sort=True, handle=handle)
        """

        class FileWrapper(object):
            def __init__(self, name, sort, parent):
                self.name = name
                self.sort = sort
                self.parent = parent

                if self.sort and self.name.endswith(".sam"):
                    self.parent.logger.info(
                        "WARNING: cannot sort SAM file {file}".format(file=self)
                    )
                    self.sort = False

            def __repr__(self):
                return self.name

            def __eq__(self, other):
                eq = self.name == other.name

                if eq and other.sort != self.sort:
                    self.parent.logger.info(
                        "WARNING: output file '{file}' chosen multiple times with differing sort options. Will set "
                        "file to be sorted".format(file=self)
                    )
                    self.sort = other.sort = True
                return eq

            def __hash__(self):
                return hash(self.name)

        for file, handles in self.output_files.items():
            if handle in handles and file_name != file.name:
                self.log_error(
                    exception=self.Exception(
                        "attempting to register output type '{type}' with more than one file: '{f1}' and '{f2}'".format(
                            type=self.output.description(handle=handle),
                            f1=file,
                            f2=file_name,
                        )
                    )
                )

        self.output_files.setdefault(
            FileWrapper(name=file_name, sort=sort, parent=self), []
        ).append(handle)

    class RestrictionSites(object):
        class Exception(Exception):
            pass

        @staticmethod
        def enzymes_by_company(comp):
            # https://www.biorxiv.org/content/10.1101/659623v1.full.pdf,
            # https://www.qiagen.com/jp/resources/download.aspx?id=813c3f19-b24e-426e-9fc6-23f48defd828&lang=en
            all_enz = {
                "ARIMA_V2": ("Arima Genomics version 2", ("Arima", "DdeI", "MseI")),
                "ARIMA": ("Arima Genomics", ("HinfI", "DpnII")),
                "DOVETAIL": ("Dovetail Genomics", ("DpnII",)),
                "PHASE": ("Phase Genomics", ("Sau3AI",)),
                "QIAGEN": ("Qiagen", ("DpnII",)),
            }
            return all_enz.setdefault(comp.upper(), None)

        @staticmethod
        def rev_comp(string):
            return string[::-1].translate(
                str.maketrans("ACGTWSMKRYBDHVN^_", "TGCAWSKMYRVHDBN_^")
            )

        @staticmethod
        def replace_degen(string):
            return string.translate(
                str.maketrans(
                    dict(
                        [
                            ("_", ""),
                            ("^", ""),
                            ("N", "[ACGT]"),
                            ("V", "[ACG]"),
                            ("H", "[ACT]"),
                            ("D", "[AGT]"),
                            ("B", "[CGT]"),
                            ("W", "[AT]"),
                            ("S", "[CG]"),
                            ("M", "[AC]"),
                            ("K", "[GT]"),
                            ("R", "[AG]"),
                            ("Y", "[CT]"),
                        ]
                    )
                )
            )

        @staticmethod
        def rem_loc_markers(string):
            return string.replace("^", "").replace("_", "")

        @classmethod
        def is_pal(cls, string, strip=False):
            st = cls.rem_loc_markers(string) if strip else string
            return st == cls.rev_comp(st)

        def __init__(self, logger):
            class Batch(RestrictionBatch):
                def get_wrapped(self, enz):
                    class EnzymeWrapper(object):
                        def __init__(self, enz):
                            self.enz = enz

                        def __getattr__(self, item):
                            return getattr(self.enz, item)

                        def __str__(self):
                            return str(self.enz)

                        @property
                        def restriction_sites(self):
                            sites = set()
                            if self.cut_once():
                                sites.add(self.elucidate())
                            else:

                                class two_cutter_wrapper(object):
                                    def __init__(self, enz):
                                        self.enz = enz
                                        self.cut_twice_fn = None
                                        self.fst5_store = None
                                        self.fst3_store = None

                                    def __enter__(self):
                                        self.cut_twice_fn = self.enz.cut_twice
                                        self.fst5_store = self.enz.fst5
                                        self.fst3_store = self.enz.fst3

                                        def ret_false():
                                            return False

                                        self.enz.cut_twice = ret_false

                                        class ActiveWrapper(object):
                                            def __init__(self, enz):
                                                self.enz = enz
                                                self.locs = (
                                                    (enz.fst5, enz.fst3),
                                                    (enz.scd5, enz.scd3),
                                                )

                                            def elucidate(self, i):
                                                if i == 1:
                                                    self.enz.fst5 = self.locs[0][0]
                                                    self.enz.fst3 = self.locs[0][1]
                                                    return self.enz.elucidate()
                                                else:
                                                    self.enz.fst5 = self.locs[1][0]
                                                    self.enz.fst3 = self.locs[1][1]
                                                    return self.enz.elucidate()

                                        return ActiveWrapper(self.enz)

                                    def __exit__(self, exc_type, exc_val, exc_tb):
                                        self.enz.cut_twice = self.cut_twice_fn
                                        self.enz.fst5 = self.fst5_store
                                        self.enz.fst3 = self.fst3_store

                                with two_cutter_wrapper(self.enz) as enz:
                                    sites.add(enz.elucidate(1))
                                    sites.add(enz.elucidate(2))

                            return sites

                    return EnzymeWrapper(self.get(enz))

            self.logger = logger
            self.all = Batch()
            for enz in AllEnzymes:
                if not enz.is_unknown():
                    self.all.add(enz)

        def create_sites(self, string):
            sites = []

            for st in string.split(","):
                st = st.strip()
                if re.match(r"^[ACGTWSMKRYBDHVN]+$", self.rem_loc_markers(st.upper())):
                    st = st.upper()

                    if "^" not in st:
                        raise self.Exception(
                            "no forward cut-location '^' in cut-site definition '{site}'".format(
                                site=st
                            )
                        )

                    if "_" not in st:
                        if not self.is_pal(st, strip=True):
                            raise self.Exception(
                                "cut-site definition '{site}' has no reverse cut-location '_' and is not palindromic".format(
                                    site=st
                                )
                            )
                        else:
                            idx = len(st) - st.index("^")
                            st = st[:idx] + "_" + st[idx:]

                    if self.is_pal(st):
                        self.logger.info(
                            "Palindromic cut-site '{site}', adding one site-definition".format(
                                site=st
                            )
                        )
                        sites += self.RestrictionSite.create_fwd_site(
                            st, logger=self.logger
                        )
                    else:
                        self.logger.info(
                            "Non-palindromic cut-site '{site}', adding two site-definitions".format(
                                site=st
                            )
                        )
                        sites += self.RestrictionSite.create_fwd_rev_sites(
                            st, logger=self.logger
                        )

                elif st in self.all:
                    enz = self.all.get_wrapped(st)
                    restriction_sites = enz.restriction_sites

                    if enz.cut_twice():
                        self.logger.info(
                            "Dual-cutter '{enz}', adding two cut-sites: {sites}".format(
                                enz=enz, sites=restriction_sites
                            )
                        )
                    else:
                        self.logger.info(
                            "Single-cutter '{enz}', adding one cut-site: {sites}".format(
                                enz=enz, sites=restriction_sites
                            )
                        )

                    sites += list(self.create_sites(string=",".join(restriction_sites)))

                elif (enz := self.enzymes_by_company(st)) is not None:
                    self.logger.info(
                        "'{comp}' enzyme mix: {mix}".format(comp=enz[0], mix=enz[1])
                    )
                    sites += list(self.create_sites(string=",".join(enz[1])))
                else:
                    raise self.Exception(
                        "Unknown enzyme or company '{enz}'".format(enz=st)
                    )

            return set(sites)

        class RestrictionSite(object):
            def __init__(self, site, logger):
                self.loc = site.index("^")
                self.rev_loc = site.index("_")

                if self.loc < self.rev_loc:
                    self.rev_loc -= 1
                else:
                    self.loc -= 1

                self.pattern = Pipeline.RestrictionSites.replace_degen(site)

                logger.info(
                    "Created site-definition '{site}' cutting locally at index {i}".format(
                        site=Pipeline.RestrictionSites.rem_loc_markers(site), i=self.loc
                    )
                )

            def __eq__(self, other):
                return self.loc == other.loc and self.pattern == other.pattern

            def __hash__(self):
                return hash(self.pattern) + hash(self.loc)

            @classmethod
            def create_fwd_rev_sites(cls, string, logger):
                return [
                    cls(string, logger=logger),
                    cls(Pipeline.RestrictionSites.rev_comp(string), logger=logger),
                ]

            @classmethod
            def create_fwd_site(cls, string, logger):
                return [cls(string, logger=logger)]

    class LoggerHandle(object):
        def __init__(self, logger):
            self.threadsAndHandles = []

            self.infoWrite = self.add_logger(log_func=logger.info, id="main")
            self.errorWrite = self.add_logger(log_func=logger.error, id="main")

        def add_logger(self, log_func, id):
            read, write = os.pipe()

            def _thread_func():
                with os.fdopen(read, encoding="utf-8", errors="replace") as file:
                    for line in file:
                        log_func("({id}) {mess}".format(id=id, mess=line[:-1]))

            thread = Thread(target=_thread_func)
            thread.start()
            self.threadsAndHandles.append((thread, write))

            return write

        class HandleWrapper(object):
            def __init__(self, info, error):
                self.info = info
                self.error = error

        def __enter__(self):
            return self.HandleWrapper(info=self.infoWrite, error=self.errorWrite)

        def __exit__(self, exc_type, exc_val, exc_tb):
            for thread, handle in self.threadsAndHandles:
                os.close(handle)
                thread.join()

    class RunParams(object):
        def __init__(
            self,
            name,
            version,
            sam_input,
            sam_output,
            fasta_input,
            restriction_sites,
            n_threads,
            min_map_q,
            info,
            error,
        ):
            self.name = name
            self.version = version
            self.samInput = sam_input
            self.samOutput = sam_output
            self.fastaInput = fasta_input
            self.restrictionSites = restriction_sites
            self.nThreads = n_threads
            self.minMapQ = min_map_q
            self.info = info
            self.error = error

    class IndexRun(RunParams):
        def __init__(self):
            super().__init__(*tuple([None] * 10))

    class LoggerWrapper(object):
        def __init__(self):
            self._logger = None

        @property
        def logger(self):
            return self._logger

        @logger.setter
        def logger(self, logger):
            if logger is not None and not isinstance(logger, logging.Logger):
                raise TypeError("logger must be an instance of 'logging.Logger'")
            self._logger = logger

        def __getattr__(self, attr):
            if attr in self.__dict__:
                return self.__dict__[attr]
            elif "_logger" not in self.__dict__:
                raise AttributeError
            elif self._logger is not None:
                return getattr(self._logger, attr)
            else:
                return self.null

        def null(self, *args, **kwargs):
            pass

    def log_error(self, exception):
        self.logger.error(exception)

        if self.exception_callback is not None:
            self.exception_callback(exception)

    def run(self):
        """
        Runs the pipeline

        :return: (Pipeline.Stats) Mapping statistics object.
        """
        self.logger.info(NAME + " " + VERSION)
        self.logger.info("Starting pipeline...")

        stats = None
        with self as runParams:
            if runParams is not None and not isinstance(runParams, self.IndexRun):
                stats = self.Stats(
                    results=_HiLine_Main(params=runParams), logger=self.logger
                )

        self.logger.info("Pipeline finished")

        if stats is not None and self.save_stats_path is not None:
            self.logger.info("Saving stats...")
            stats.save(
                base_path=self.save_stats_path,
                save_stats_file=True,
                generate_csvs=True,
                generate_figs=True,
            )
            self.logger.info("Stats saved")

        return stats

    class Stats(object):
        """
        Mapping statistics object.
        """

        @staticmethod
        def load_stats_file(path):
            """
            Loads a saved (pickled) statistics object.

            :param path: (str) Path to the saved object.
            :return: (Pipeline.Stats) Statistics object.
            """
            import lzma
            import pickle

            with lzma.open(path, "rb") as file:
                return pickle.load(file)

        class Sequence(object):
            def __init__(
                self,
                name,
                length,
                n_restriction_sites,
                self_cicle,
                dangling_end,
                same_frag_and_strand,
                re_ligation,
                ff,
                fr,
                rf,
                rr,
            ):
                self.name = name
                self.length = length
                self.nRestrictionSites = n_restriction_sites

                self.selfCircle = self_cicle
                self.danglingEnd = dangling_end
                self.sameFragAndStrand = same_frag_and_strand
                self.reLigation = re_ligation

                self._FF = ff
                self._FR = fr
                self._RF = rf
                self._RR = rr

            def __len__(self):
                return self.length

            @property
            def invalid(self):
                return (
                    self.selfCircle
                    + self.danglingEnd
                    + self.sameFragAndStrand
                    + self.reLigation
                )

            @property
            def norm_self_circle(self):
                return self.selfCircle / len(self)

            @property
            def norm_dangling_end(self):
                return self.danglingEnd / len(self)

            @property
            def norm_same_frag_and_strand(self):
                return self.sameFragAndStrand / len(self)

            @property
            def norm_re_ligation(self):
                return self.reLigation / len(self)

            @property
            def norm_n_restriction_sites(self):
                return self.nRestrictionSites / len(self)

            @property
            def norm_invalid(self):
                return self.invalid / len(self)

            def _get_ff(self, other_index):
                return self._FF[other_index]

            def _get_fr(self, other_index):
                return self._FR[other_index]

            def _get_rf(self, other_index):
                return self._RF[other_index]

            def _get_rr(self, other_index):
                return self._RR[other_index]

            @property
            def data_frame(self):
                return pd.DataFrame(
                    {
                        "statistic": ["length"]
                        + (
                            [
                                "restriction sites",
                                "self-circle pairs",
                                "dangling-end pairs",
                                "same fragment and strand pairs",
                                "re-ligated pairs",
                                "total invalid",
                            ]
                            * 2
                        ),
                        "count": [
                            self.length,
                            self.nRestrictionSites,
                            self.selfCircle,
                            self.danglingEnd,
                            self.sameFragAndStrand,
                            self.reLigation,
                            self.invalid,
                            self.norm_n_restriction_sites,
                            self.norm_self_circle,
                            self.norm_dangling_end,
                            self.norm_same_frag_and_strand,
                            self.norm_re_ligation,
                            self.norm_invalid,
                        ],
                        "name": [self.name] * 13,
                        "type": ["basic"] + (["raw count"] * 6) + (["norm count"] * 6),
                        "total": ([False] * 6) + [True] + ([False] * 5) + [True],
                        "pair": ([False] * 2) + ([True] * 5) + ([False]) + ([True] * 5),
                    }
                )

        class ReadCounts(object):
            def __init__(
                self,
                n_invalid_ref,
                n_un_paired,
                n_supplementary,
                n_duplicate,
                n_qc_fail,
                n_secondary,
                n_un_mapped,
                n_below_min_mapq,
                n_too_far_from_restriction_site,
                n_dumped_reads,
                n_good_reads,
            ):
                self.nInvalidRef = n_invalid_ref
                self.nUnPaired = n_un_paired
                self.nSupplementary = n_supplementary
                self.nDuplicate = n_duplicate
                self.nQCFail = n_qc_fail
                self.nSecondary = n_secondary
                self.nUnMapped = n_un_mapped
                self.nBelowMinMapq = n_below_min_mapq
                self.nTooFarFromRestrictionSite = n_too_far_from_restriction_site
                self.nDumpedReads = n_dumped_reads
                self.nGoodReads = n_good_reads

            @property
            def data_frame(self):
                return pd.DataFrame(
                    {
                        "type": [
                            "invalid reference",
                            "un-paired",
                            "supplementary",
                            "duplicate",
                            "QC-failed",
                            "secondary",
                            "un-mapped",
                            "mapping quality too low",
                            "too far from restriction site",
                            "dumped",
                            "good read",
                        ],
                        "counts": [
                            self.nInvalidRef,
                            self.nUnPaired,
                            self.nSupplementary,
                            self.nDuplicate,
                            self.nQCFail,
                            self.nSecondary,
                            self.nUnMapped,
                            self.nBelowMinMapq,
                            self.nTooFarFromRestrictionSite,
                            self.nDumpedReads,
                            self.nGoodReads,
                        ],
                        "class": ["invalid" for _ in range(10)] + ["valid_pairs"],
                    }
                )

        class PairCounts(object):
            def __init__(
                self,
                self_circle,
                dangling_end,
                same_frag_and_strand,
                re_ligation,
                ff,
                fr,
                rf,
                rr,
            ):
                self.selfCircle = self_circle
                self.danglingEnd = dangling_end
                self.sameFragAndStrand = same_frag_and_strand
                self.reLigation = re_ligation
                self.FF = ff
                self.FR = fr
                self.RF = rf
                self.RR = rr

            @property
            def data_frame(self):
                return pd.DataFrame(
                    {
                        "type": [
                            "self-circle",
                            "dangling-end",
                            "same fragment and strand",
                            "re-ligation",
                            "valid_pairs ff",
                            "valid_pairs fr",
                            "valid_pairs rf",
                            "valid_pairs rr",
                        ],
                        "counts": [
                            self.selfCircle,
                            self.danglingEnd,
                            self.sameFragAndStrand,
                            self.reLigation,
                            self.FF,
                            self.FR,
                            self.RF,
                            self.RR,
                        ],
                        "class": ["invalid" for _ in range(4)]
                        + ["valid_pairs" for _ in range(4)],
                    }
                )

            @property
            def valid(self):
                return self.FF + self.FR + self.RF + self.RR

            @property
            def invalid(self):
                return (
                    self.selfCircle
                    + self.danglingEnd
                    + self.sameFragAndStrand
                    + self.reLigation
                )

        class RestrictionSiteDistances(object):
            def __init__(
                self,
                self_circle,
                dangling_end,
                same_frag_and_strand,
                re_ligation,
                ff,
                fr,
                rf,
                rr,
            ):
                self.selfCircle = self_circle
                self.danglingEnd = dangling_end
                self.sameFragAndStrand = same_frag_and_strand
                self.reLigation = re_ligation
                self.FF = ff
                self.FR = fr
                self.RF = rf
                self.RR = rr

            @property
            def valid(self):
                return np.concatenate((self.FF, self.FR, self.RF, self.RR))

            @property
            def invalid(self):
                return np.concatenate(
                    (
                        self.selfCircle,
                        self.danglingEnd,
                        self.sameFragAndStrand,
                        self.reLigation,
                    )
                )

            @property
            def valid_f(self):
                return np.concatenate((self.FF, self.FR[::2], self.RF[1::2]))

            @property
            def valid_r(self):
                return np.concatenate((self.FR[1::2], self.RF[::2], self.RR))

            @property
            def self_circle_and_dangling_end_f(self):
                return np.concatenate((self.selfCircle[1::2], self.danglingEnd[::2]))

            @property
            def self_circle_and_dangling_end_r(self):
                return np.concatenate((self.selfCircle[::2], self.danglingEnd[1::2]))

            @property
            def data_frame(self):
                return pd.DataFrame(
                    {
                        "class": (["ff"] * len(self.FF))
                        + (["fr"] * len(self.FR))
                        + (["rf"] * len(self.RF))
                        + (["rr"] * len(self.RR))
                        + (["sc"] * len(self.selfCircle))
                        + (["de"] * len(self.danglingEnd))
                        + (["ss"] * len(self.sameFragAndStrand))
                        + (["rl"] * len(self.reLigation)),
                        "valid": (
                            [True]
                            * (
                                len(self.FF)
                                + len(self.FR)
                                + len(self.RF)
                                + len(self.RR)
                            )
                        )
                        + (
                            [False]
                            * (
                                len(self.selfCircle)
                                + len(self.danglingEnd)
                                + len(self.sameFragAndStrand)
                                + len(self.reLigation)
                            )
                        ),
                        "dir": (["F"] * len(self.FF))
                        + (["F", "R"] * (len(self.FR) // 2))
                        + (["R", "F"] * (len(self.RF) // 2))
                        + (["R"] * len(self.RR))
                        + (["R", "F"] * (len(self.selfCircle) // 2))
                        + (["F", "R"] * (len(self.danglingEnd) // 2))
                        + (
                            ["UKN"]
                            * (len(self.sameFragAndStrand) + len(self.reLigation))
                        ),
                        "id": np.arange(
                            len(self.FF)
                            + len(self.FR)
                            + len(self.RF)
                            + len(self.RR)
                            + len(self.selfCircle)
                            + len(self.danglingEnd)
                            + len(self.sameFragAndStrand)
                            + len(self.reLigation)
                        )
                        // 2,
                        "dis": np.concatenate(
                            (
                                self.FF,
                                self.FR,
                                self.RF,
                                self.RR,
                                self.selfCircle,
                                self.danglingEnd,
                                self.sameFragAndStrand,
                                self.reLigation,
                            )
                        ),
                    }
                )

        class IntraSeqDistances(object):
            def __init__(self, ff, fr, rf, rr):
                self.FF = ff
                self.FR = fr
                self.RF = rf
                self.RR = rr

            @property
            def all(self):
                return np.concatenate((self.FF, self.FR, self.RF, self.RR))

            @property
            def data_frame(self):
                return pd.DataFrame(
                    {
                        "class": (["ff"] * len(self.FF))
                        + (["fr"] * len(self.FR))
                        + (["rf"] * len(self.RF))
                        + (["rr"] * len(self.RR)),
                        "dis": self.all,
                    }
                )

            def calculate_histograms_and_fit_models(self):
                """
                Calculates pair-separation histograms and fits various models.

                :return:    (pandas.DataFrame) Histograms and model fits.
                            (pandas.DataFrame) Model-vs-model F-tests.
                """

                class Model(object):
                    def __init__(
                        self, poly_coef, prob_data, sep_data, intercept=None, name=None
                    ):
                        if type(poly_coef) == int or type(poly_coef) == float:
                            poly_coef = [poly_coef]

                        self.polyCoef = poly_coef
                        self._intercept = intercept
                        self.probData = prob_data
                        self.sepData = sep_data
                        self.name = name

                    def ftest(self, other):
                        from scipy.stats import f

                        assert np.all(self.sepData == other.sepData) and np.all(
                            self.probData == other.probData
                        )

                        if self.n_param > other.n_param:
                            model1 = other
                            model2 = self
                        else:
                            return np.nan

                        dom = model2.n_param - model1.n_param
                        dof = model2.n_data - model2.n_param
                        if dof < 1:
                            return np.nan

                        return f.sf(
                            ((model1.rss - model2.rss) / dom) / (model2.rss / dof),
                            dom,
                            dof,
                        )

                    @property
                    def poly(self):
                        poly = np.zeros(len(self.sepData))
                        for index, coef in enumerate(self.polyCoef):
                            poly += (self.sepData ** (index + 1)) * coef

                        return poly

                    @property
                    def intercept(self):
                        if self._intercept is not None:
                            return self._intercept
                        else:
                            return np.mean(self.probData - self.poly)

                    @property
                    def n_data(self):
                        return len(self.sepData)

                    @property
                    def n_param(self):
                        return (
                            1 if self._intercept is None else (len(self.polyCoef) + 1)
                        )

                    @property
                    def model_probs(self):
                        return self.intercept + self.poly

                    @property
                    def rss(self):
                        return ((self.probData - self.model_probs) ** 2).sum()

                    # Naumova N, Imakaev M, Fudenberg G, Zhan Y, Lajoie BR, Mirny LA, Dekker J.
                    # Organization of the Mitotic Chromosome. Science. 2013;342:948–53.
                    @classmethod
                    def inter_phase_model(cls, prob_data, sep_data):
                        class InterphaseModel(cls):
                            def __init__(self):
                                super().__init__(
                                    poly_coef=-1,
                                    prob_data=prob_data,
                                    sep_data=sep_data,
                                    name="I-phase model",
                                )

                        return InterphaseModel()

                    # Naumova N, Imakaev M, Fudenberg G, Zhan Y, Lajoie BR, Mirny LA, Dekker J.
                    # Organization of the Mitotic Chromosome. Science. 2013;342:948–53.
                    @classmethod
                    def mitotic_model(cls, prob_data, sep_data):
                        class MitoticModel(cls):
                            def __init__(self):
                                super().__init__(
                                    poly_coef=-0.5,
                                    prob_data=prob_data,
                                    sep_data=sep_data,
                                    name="M-phase model",
                                )

                        return MitoticModel()

                    @classmethod
                    def fit(cls, prob_data, sep_data, deg):
                        from numpy.polynomial.polynomial import Polynomial as Poly

                        poly_coef = np.nan * np.zeros(deg + 1)
                        if len(prob_data) > 1:
                            poly_coef = (
                                Poly.fit(sep_data, prob_data, deg).convert().coef
                            )

                        return cls(
                            poly_coef=poly_coef[1:],
                            intercept=poly_coef[0],
                            prob_data=prob_data,
                            sep_data=sep_data,
                            name="poly model deg {n}".format(n=deg),
                        )

                    @classmethod
                    def fit_all(cls, prob_data, sep_data, deg_max=2):
                        return [
                            cls.fit(prob_data=prob_data, sep_data=sep_data, deg=n)
                            for n in range(deg_max + 1)
                        ] + [
                            cls.inter_phase_model(
                                prob_data=prob_data, sep_data=sep_data
                            ),
                            cls.mitotic_model(prob_data=prob_data, sep_data=sep_data),
                        ]

                def perfrom_fit(df, cls, prob_in=None, sep=None, bins=None):
                    linear_fit_range = (1000, 7000000)
                    n_bins_max = 1024

                    # https://stackoverflow.com/questions/39418380/histogram-with-equal-number-of-points-in-each-bin
                    def histedges_equal_n(x):
                        nbin = np.clip(
                            len(np.histogram_bin_edges(x, bins="auto")) - 1,
                            1,
                            n_bins_max,
                        )
                        npt = len(x)
                        return np.interp(
                            np.linspace(0, npt, nbin + 1), np.arange(npt), np.sort(x)
                        )

                    if prob_in is None:
                        bins = histedges_equal_n(df["dis"])
                        count = np.histogram(df["dis"], bins=bins)[0]
                        prob, sep = np.histogram(df["dis"], bins=bins, density=True)
                        sep = (sep[:-1] + sep[1:]) / 2
                    else:
                        prob = prob_in * np.histogram(df["dis"], bins=bins)[0]
                        count = 1

                    prob_df = pd.DataFrame({"prob": prob, "sep": sep})

                    select = np.logical_and(
                        prob_df["sep"] > linear_fit_range[0],
                        prob_df["sep"] < linear_fit_range[1],
                    )
                    models = Model.fit_all(
                        prob_data=np.log(prob_df["prob"][select]),
                        sep_data=np.log(prob_df["sep"][select]),
                    )

                    prob_df["type"] = "histogram"
                    ftest_df = pd.DataFrame({"model1": [], "model2": [], "p-value": []})
                    for model in models:
                        prob_df = pd.concat(
                            (
                                prob_df,
                                pd.DataFrame(
                                    {
                                        "prob": np.exp(model.model_probs),
                                        "sep": np.exp(model.sepData),
                                        "type": model.name,
                                    }
                                ),
                            )
                        )
                        for model2 in models:
                            ftest_df = pd.concat(
                                (
                                    ftest_df,
                                    pd.DataFrame(
                                        {
                                            "model1": [model.name],
                                            "model2": [model2.name],
                                            "p-value": [model.ftest(other=model2)],
                                        }
                                    ),
                                )
                            )

                    ftest_df["class"] = cls
                    prob_df["class"] = cls

                    return prob_df, ftest_df, prob / count, sep, bins

                df = self.data_frame
                model_df, ftest_df, norm_prob, sep, bins = perfrom_fit(df=df, cls="all")
                for cls in set(df["class"]):
                    dfs = perfrom_fit(
                        df=df[df["class"] == cls],
                        cls=cls,
                        prob_in=norm_prob,
                        sep=sep,
                        bins=bins,
                    )
                    model_df = pd.concat((model_df, dfs[0]))
                    ftest_df = pd.concat((ftest_df, dfs[1]))

                return model_df, ftest_df

        def __init__(self, logger, results):
            self.logger = logger

            (
                type_counts,
                seq_names,
                seq_lengths,
                restriction_sites_per_seq,
                self_circle_pairs_per_seq,
                dangling_end_pairs_per_seq,
                same_fragment_and_strand_pairs_per_seq,
                re_ligated_pairs_per_seq,
                ff_per_seq,
                fr_per_seq,
                rf_per_seq,
                rr_per_seq,
                self_circle_restriction_site_distances,
                dangling_end_restriction_site_distances,
                same_fragment_and_strand_restriction_site_distances,
                re_ligated_restriction_site_distances,
                ff_restriction_site_distances,
                fr_restriction_site_distances,
                rf_restriction_site_distances,
                rr_restriction_site_distances,
                ff_intra_sequence_read_distances,
                fr_intra_sequence_read_distances,
                rf_intra_sequence_read_distances,
                rr_intra_sequence_read_distances,
                ff_intra_sequence_fragment_distances,
                fr_intra_sequence_fragment_distances,
                rf_intra_sequence_fragment_distances,
                rr_intra_sequence_fragment_distances,
            ) = results

            n_good_reads = int(2 * type_counts[-8:].sum())
            self.reads = self.ReadCounts(
                n_invalid_ref=type_counts[0],
                n_un_paired=type_counts[1],
                n_supplementary=type_counts[2],
                n_duplicate=type_counts[3],
                n_qc_fail=type_counts[4],
                n_secondary=type_counts[5],
                n_un_mapped=type_counts[6],
                n_below_min_mapq=type_counts[7],
                n_too_far_from_restriction_site=type_counts[8],
                n_dumped_reads=int(type_counts[9] - n_good_reads),
                n_good_reads=n_good_reads,
            )

            self.pairs = self.PairCounts(
                self_circle=type_counts[10],
                dangling_end=type_counts[11],
                same_frag_and_strand=type_counts[12],
                re_ligation=type_counts[13],
                ff=type_counts[14],
                fr=type_counts[15],
                rf=type_counts[16],
                rr=type_counts[17],
            )

            self.sequences = [
                self.Sequence(
                    name=name,
                    length=length,
                    n_restriction_sites=sites,
                    self_cicle=sc,
                    dangling_end=de,
                    same_frag_and_strand=sfs,
                    re_ligation=rl,
                    ff=ff,
                    fr=fr,
                    rf=rf,
                    rr=rr,
                )
                for name, length, sites, sc, de, sfs, rl, ff, fr, rf, rr in zip(
                    seq_names,
                    seq_lengths,
                    restriction_sites_per_seq,
                    self_circle_pairs_per_seq,
                    dangling_end_pairs_per_seq,
                    same_fragment_and_strand_pairs_per_seq,
                    re_ligated_pairs_per_seq,
                    ff_per_seq,
                    fr_per_seq,
                    rf_per_seq,
                    rr_per_seq,
                )
            ]

            self.restrictionSiteDistances = self.RestrictionSiteDistances(
                self_circle=self_circle_restriction_site_distances,
                dangling_end=dangling_end_restriction_site_distances,
                same_frag_and_strand=same_fragment_and_strand_restriction_site_distances,
                re_ligation=re_ligated_restriction_site_distances,
                ff=ff_restriction_site_distances,
                fr=fr_restriction_site_distances,
                rf=rf_restriction_site_distances,
                rr=rr_restriction_site_distances,
            )

            self.intraSeqReadDistances = self.IntraSeqDistances(
                ff=ff_intra_sequence_read_distances,
                fr=fr_intra_sequence_read_distances,
                rf=rf_intra_sequence_read_distances,
                rr=rr_intra_sequence_read_distances,
            )

            self.intraSeqFragmentDistances = self.IntraSeqDistances(
                ff=ff_intra_sequence_fragment_distances,
                fr=fr_intra_sequence_fragment_distances,
                rf=rf_intra_sequence_fragment_distances,
                rr=rr_intra_sequence_fragment_distances,
            )

        @property
        def n_sequences(self):
            return len(self.sequences)

        @property
        def sequence_lengths(self):
            return np.array([seq.length for seq in self.sequences])

        @property
        def total_genome_length(self):
            return self.sequence_lengths.sum()

        def ff_count(self, seq1, seq2, norm=True):
            index1 = self.sequences.index(seq1)
            index2 = self.sequences.index(seq2)
            min_idx = min(index1, index2)
            max_idx = max(index1, index2)

            count = self.sequences[min_idx]._get_ff(max_idx - min_idx)

            if norm:
                count /= len(self.sequences[min_idx]) * len(self.sequences[max_idx])

            return count

        def fr_count(self, seq1, seq2, norm=True):
            index1 = self.sequences.index(seq1)
            index2 = self.sequences.index(seq2)
            min_idx = min(index1, index2)
            max_idx = max(index1, index2)

            count = self.sequences[min_idx]._get_fr(max_idx - min_idx)

            if norm:
                count /= len(self.sequences[min_idx]) * len(self.sequences[max_idx])

            return count

        def rf_count(self, seq1, seq2, norm=True):
            index1 = self.sequences.index(seq1)
            index2 = self.sequences.index(seq2)
            min_idx = min(index1, index2)
            max_idx = max(index1, index2)

            count = self.sequences[min_idx]._get_rf(max_idx - min_idx)

            if norm:
                count /= len(self.sequences[min_idx]) * len(self.sequences[max_idx])

            return count

        def rr_count(self, seq1, seq2, norm=True):
            index1 = self.sequences.index(seq1)
            index2 = self.sequences.index(seq2)
            min_idx = min(index1, index2)
            max_idx = max(index1, index2)

            count = self.sequences[min_idx]._get_rr(max_idx - min_idx)

            if norm:
                count /= len(self.sequences[min_idx]) * len(self.sequences[max_idx])

            return count

        def valid_count(self, seq1, seq2, norm=True):
            index1 = self.sequences.index(seq1)
            index2 = self.sequences.index(seq2)
            min_idx = min(index1, index2)
            max_idx = max(index1, index2)

            count = (
                self.sequences[min_idx]._get_ff(max_idx - min_idx)
                + self.sequences[min_idx]._get_fr(max_idx - min_idx)
                + self.sequences[min_idx]._get_rf(max_idx - min_idx)
                + self.sequences[min_idx]._get_rr(max_idx - min_idx)
            )

            if norm:
                count /= len(self.sequences[min_idx]) * len(self.sequences[max_idx])

            return count

        @property
        def sequences_data_frame(self):
            return pd.concat([seq.data_frame for seq in self.sequences])

        @property
        def valid_pairs_data_frame(self):
            id1 = []
            id2 = []
            types = []
            norm = []
            count = []

            for index, seq in enumerate(self.sequences):
                ff = seq._FF
                fr = seq._FR
                rf = seq._RF
                rr = seq._RR
                m = len(ff)
                id1 += [seq.name] * 10 * m
                id2 += [s.name for s in self.sequences[index:]] * 10
                types += (
                    (["ff"] * 2 * m)
                    + (["fr"] * 2 * m)
                    + (["rf"] * 2 * m)
                    + (["rr"] * 2 * m)
                    + (["total"] * 2 * m)
                )
                norm += (([False] * m) + ([True] * m)) * 5
                lens = np.array([len(s) for s in self.sequences[index:]]) * len(seq)
                count += (
                    (np.concatenate((ff, ff / lens))).tolist()
                    + (np.concatenate((fr, fr / lens))).tolist()
                    + (np.concatenate((rf, rf / lens))).tolist()
                    + (np.concatenate((rr, rr / lens))).tolist()
                    + (
                        np.concatenate(
                            ((ff + fr + rf + rr), (ff + fr + rf + rr) / lens)
                        )
                    ).tolist()
                )

            return pd.DataFrame(
                {"id1": id1, "id2": id2, "type": types, "norm": norm, "count": count}
            )

        @property
        def inter_sequence_indexes(self):
            df = self.valid_pairs_data_frame
            df = df[(df["type"] == "total") & (df["norm"] == False)]
            df = pd.concat(
                (
                    df,
                    df[df["id1"] != df["id2"]].rename(
                        columns={"id1": "id2", "id2": "id1"}
                    ),
                ),
                sort=False,
            )
            n = df.pivot("id1", "id2", "count").sum()
            ni = df[df["id1"] != df["id2"]].pivot("id1", "id2", "count").sum()

            df = self.sequences_data_frame
            seq_lengths = df[df["statistic"] == "length"].set_index("name")["count"]

            idx = (self.total_genome_length * ni) / (n * seq_lengths)
            return pd.DataFrame({"name": idx.index, "index": np.array(idx)})

        @property
        def inter_intra_valid_pairs_data_frame(self):
            totals = np.array(
                [self.pairs.FF, self.pairs.FR, self.pairs.RF, self.pairs.RR]
            )
            intra = np.array(
                [
                    len(self.intraSeqReadDistances.FF),
                    len(self.intraSeqReadDistances.FR),
                    len(self.intraSeqReadDistances.RF),
                    len(self.intraSeqReadDistances.RR),
                ]
            )

            return pd.DataFrame(
                {
                    "class": ["ff", "fr", "rf", "rr"] * 2,
                    "type": (["intra"] * 4) + (["inter"] * 4),
                    "count": np.concatenate((intra, totals - intra)),
                }
            )

        def save(
            self,
            base_path,
            save_stats_file=False,
            generate_csvs=False,
            generate_figs=True,
        ):
            """
            Saves statistics.

            :param base_path: (str) Base save path.
            :param save_stats_file: (bool) Whether or not to save (pickle) the entire statistics object.
            :param generate_csvs: (bool) Whether or not to generate and save csv files.
            :param generate_figs: (bool) Whether or not to generate and save figures.
            """
            _showwarning = warnings.showwarning
            with warnings.catch_warnings():
                warnings.showwarning = _showwarning
                warnings.filterwarnings("ignore", message="divide by zero")
                warnings.filterwarnings(
                    "ignore", message="invalid value encountered in double_scalars"
                )

                if not (save_stats_file or generate_csvs or generate_figs):
                    return

                if generate_figs:
                    import matplotlib as mpl
                    import matplotlib.pyplot as plt

                    mpl.use("agg")
                    font = {"family": "sans", "weight": "normal", "size": 14}
                    mpl.rc("font", **font)
                    mpl.rcParams["agg.path.chunksize"] = 10000
                    plt.style.use("ggplot")
                    plt.rcParams["figure.figsize"] = (10.0, 7.0)

                    import seaborn as sns

                    sns.set()
                    sns.set(style="darkgrid")
                    sns.set(color_codes=True)

                def log(msg):
                    self.logger.info(f"(stats) {msg}")

                def zero_bandwidth_wrapper(func, *args, **kwargs):
                    try:
                        return func(*args, **kwargs)
                    except RuntimeError as rte:
                        if str(rte).startswith("Selected KDE bandwidth is 0"):
                            kwargs["kde_kws"] = {"bw": 0.1}
                            return func(*args, **kwargs)
                        else:
                            raise rte

                def process_basic_stats(
                    df, dir, name, x="counts", y1="type", y2="class"
                ):
                    if generate_csvs:
                        df.to_csv(
                            path_or_buf=dir(name="{name}.csv.bz2".format(name=name)),
                            index=False,
                        )

                    if generate_figs:
                        figure, axes = plt.subplots(2, 1)
                        sns.barplot(x=x, y=y1, data=df, ax=axes[0])
                        sns.barplot(
                            x=x, y=y2, data=df, ax=axes[1], estimator=sum, ci=None
                        )

                        figure.savefig(
                            dir(name="{name}.png".format(name=name)),
                            dpi=200,
                            bbox_inches="tight",
                        )
                        figure.clf()
                        plt.close()

                    if generate_figs or generate_csvs:
                        log("{name} stats saved".format(name=name))

                def process_sequence_stats(df, dir, name):
                    if generate_figs:
                        figure, axes = plt.subplots(2, 1)

                        df_sub = df[df.type == "raw count"]
                        sns.stripplot(data=df_sub, y="statistic", x="count", ax=axes[0])
                        sns.violinplot(
                            data=df_sub, y="statistic", x="count", ax=axes[0]
                        )
                        axes[0].get_yaxis().label.set_visible(False)

                        df_sub = df[df.type == "norm count"]
                        sns.stripplot(data=df_sub, y="statistic", x="count", ax=axes[1])
                        sns.violinplot(
                            data=df_sub, y="statistic", x="count", ax=axes[1]
                        )
                        pad = 0.1 * df_sub["count"].max()
                        axes[1].set_xlim(
                            df_sub["count"].min() - pad, df_sub["count"].max() + pad
                        )
                        axes[1].set_xlabel("count per base-pair")
                        axes[1].get_yaxis().label.set_visible(False)

                        plt.tight_layout()
                        figure.savefig(
                            dir(name="{name}.png".format(name=name)),
                            dpi=200,
                            bbox_inches="tight",
                        )
                        figure.clf()
                        plt.close()

                        log(
                            "{baseName} stats, {name} figure saved".format(
                                name=name, baseName=dir.path_name
                            )
                        )

                def process_valid_pairs_per_sequence(df, dir, name, totals=False):
                    if generate_figs:
                        if totals:
                            figure, axes = plt.subplots(1, 1)
                            axes = [axes]
                            df = df[df.type == "total"]
                            ids = ["total"]
                        else:
                            figure, axes = plt.subplots(2, 2)
                            axes = axes.flat
                            df = df[df.type != "total"]
                            ids = ("ff", "fr", "rf", "rr")

                        for ax, id in zip(axes, ids):
                            sns.heatmap(
                                df[df.type == id].pivot("id1", "id2", "count"),
                                square=True,
                                ax=ax,
                                vmin=df["count"].min(),
                                vmax=df["count"].max(),
                            )
                            ax.get_yaxis().label.set_visible(False)
                            ax.get_xaxis().label.set_visible(False)
                            ax.set_title(id)

                        plt.tight_layout()
                        figure.savefig(
                            dir(name="{name}.png".format(name=name)),
                            dpi=200,
                            bbox_inches="tight",
                        )
                        figure.clf()
                        plt.close()

                def process_site_distance_plots(
                    df, dir, name, fr_classes, other_classes
                ):
                    if generate_figs:
                        figure, axes = plt.subplots(2, 2)
                        for ax, id in zip(axes.flat, fr_classes + other_classes):
                            zero_bandwidth_wrapper(
                                sns.histplot, df[df["class"] == id].dis, ax=ax
                            )
                            ax.set_xlabel("bp from site")
                            ax.set_title(id)

                        plt.tight_layout()
                        figure.savefig(
                            dir(name="{name}.png".format(name="pair_type")),
                            dpi=200,
                            bbox_inches="tight",
                        )
                        figure.clf()
                        plt.close()

                        figure, axes = plt.subplots(1, 2)
                        for ax, id in zip(axes.flat, ("F", "R")):
                            zero_bandwidth_wrapper(
                                sns.histplot, df[df["dir"] == id].dis, ax=ax
                            )
                            ax.set_xlabel("bp from site")
                            ax.set_title(id)

                        plt.tight_layout()
                        figure.savefig(
                            dir(name="{name}.png".format(name="read_dir")),
                            dpi=200,
                            bbox_inches="tight",
                        )
                        figure.clf()
                        plt.close()

                        df_sub = df[
                            np.array([c in fr_classes for c in df["class"]])
                        ].pivot("id", "dir", "dis")
                        figure = zero_bandwidth_wrapper(
                            sns.jointplot, data=df_sub, x="F", y="R", kind="hex"
                        ).fig
                        figure.axes[0].set_xlabel("F, bp from site")
                        figure.axes[0].set_ylabel("R, bp from site")

                        plt.tight_layout()
                        figure.savefig(
                            dir(
                                name="{name}.png".format(
                                    name="{name}_bivariate".format(
                                        name="_".join(fr_classes)
                                    )
                                )
                            ),
                            dpi=200,
                            bbox_inches="tight",
                        )
                        figure.clf()
                        plt.close()

                        log(
                            "{name1}, {name2} plots saved".format(
                                name1=name, name2=dir.path_name
                            )
                        )

                def process_model_f_tests(df, dir, name, all=False):
                    if generate_figs:
                        if all:
                            figure, axes = plt.subplots(1, 1)
                            axes = [axes]
                            df = df[df["class"] == "all"]
                            ids = ["all"]
                        else:
                            figure, axes = plt.subplots(2, 2)
                            axes = axes.flat
                            df = df[df["class"] != "all"]
                            ids = ("ff", "fr", "rf", "rr")

                        for ax, id in zip(axes, ids):
                            sns.heatmap(
                                df[df["class"] == id].pivot(
                                    "model1", "model2", "p-value"
                                ),
                                square=True,
                                ax=ax,
                                vmin=df["p-value"].min(),
                                vmax=df["p-value"].max(),
                            )
                            ax.get_yaxis().label.set_visible(False)
                            ax.get_xaxis().label.set_visible(False)
                            ax.set_title(id)

                        plt.tight_layout()
                        figure.savefig(
                            dir(name="{name}.png".format(name=name)),
                            dpi=200,
                            bbox_inches="tight",
                        )
                        figure.clf()
                        plt.close()

                class Dir(object):
                    def __init__(self, path=base_path):
                        self.path = path

                    @property
                    def path_name(self):
                        return basename(normpath(self.path))

                    def make_path(self, name):
                        return join(self.path, name)

                    def enter(self, name):
                        return Dir(path=self.make_path(name=name))

                    def __enter__(self):
                        makedirs(name=self.path, exist_ok=True)
                        return self

                    def __exit__(self, exc_type, exc_val, exc_tb):
                        pass

                    def __call__(self, name):
                        return self.make_path(name=name)

                try:
                    with Dir() as dir:
                        import pickle
                        import lzma
                        from multiprocessing.pool import ThreadPool

                        fig_pool = ThreadPool(1)
                        pool = ThreadPool(3)

                        def save_file():
                            if save_stats_file:
                                log("writing stats file...")
                                with lzma.open(
                                    dir(name="stats.pickle.xz"), "wb"
                                ) as file:
                                    pickle.dump(self, file)
                                log("stats file writen")

                        pool.apply_async(save_file)

                        if generate_figs or generate_csvs:
                            with dir.enter(name="reads") as dir2:
                                dir_cache_0 = dir2

                                def thread_fn_0():
                                    process_basic_stats(
                                        df=self.reads.data_frame,
                                        dir=dir_cache_0,
                                        name=dir_cache_0.path_name,
                                    )

                                fig_pool.apply_async(thread_fn_0)

                            with dir.enter(name="HiCPairs") as dir2:
                                dir_cache_1 = dir2

                                def thread_fn_1():
                                    process_basic_stats(
                                        df=self.pairs.data_frame,
                                        dir=dir_cache_1,
                                        name=dir_cache_1.path_name,
                                    )

                                    name = "intra_inter_sequence_valid_pairs"
                                    df = self.inter_intra_valid_pairs_data_frame

                                    if generate_csvs:
                                        df.to_csv(
                                            path_or_buf=dir_cache_1(
                                                name="{name}.csv.bz2".format(name=name)
                                            ),
                                            index=False,
                                        )

                                    if generate_figs:
                                        figure, axes = plt.subplots(2, 1)
                                        sns.barplot(
                                            x="count",
                                            y="class",
                                            data=df,
                                            hue="type",
                                            ax=axes[0],
                                        )
                                        sns.barplot(
                                            x="count",
                                            y="type",
                                            data=df,
                                            ax=axes[1],
                                            estimator=sum,
                                            ci=None,
                                        )

                                        figure.savefig(
                                            dir_cache_1(
                                                name="{name}.png".format(name=name)
                                            ),
                                            dpi=200,
                                            bbox_inches="tight",
                                        )
                                        figure.clf()
                                        plt.close()

                                fig_pool.apply_async(thread_fn_1)

                            with dir.enter(name="sequences") as dir2:
                                with dir2.enter(name="per_sequence") as dir3:
                                    dir_cache_2 = dir3

                                    def thread_fn_2():
                                        df = self.sequences_data_frame

                                        if generate_csvs:
                                            df.to_csv(
                                                path_or_buf=dir_cache_2(
                                                    name="{name}.csv.bz2".format(
                                                        name=dir_cache_2.path_name
                                                    )
                                                ),
                                                index=False,
                                            )

                                        process_sequence_stats(
                                            df=df[df.pair & (df.total == False)],
                                            dir=dir_cache_2,
                                            name="invalid_pairs",
                                        )
                                        process_sequence_stats(
                                            df=df[df.pair & df.total],
                                            dir=dir_cache_2,
                                            name="total_invalid_pairs",
                                        )
                                        process_sequence_stats(
                                            df=df[
                                                (df.statistic == "restriction sites")
                                                & (df.total == False)
                                            ],
                                            dir=dir_cache_2,
                                            name="restriction_sites",
                                        )

                                        name = "inter_sequence_indexes"
                                        df = self.inter_sequence_indexes

                                        if generate_csvs:
                                            df.to_csv(
                                                path_or_buf=dir_cache_2(
                                                    name="{name}.csv.bz2".format(
                                                        name=name
                                                    )
                                                ),
                                                index=False,
                                            )

                                        if generate_figs:
                                            figure, axes = plt.subplots(2, 1)

                                            sns.stripplot(
                                                data=df, x="index", ax=axes[0]
                                            )
                                            sns.violinplot(
                                                data=df, x="index", ax=axes[0]
                                            )
                                            axes[0].get_yaxis().label.set_visible(False)
                                            axes[0].set_xlabel("inter sequence index")

                                            sns.barplot(
                                                data=df, x="index", y="name", ax=axes[1]
                                            )
                                            axes[1].get_yaxis().label.set_visible(False)

                                            axes[1].set_xlabel("inter sequence index")
                                            axes[1].get_yaxis().label.set_visible(False)

                                            plt.tight_layout()
                                            figure.savefig(
                                                dir_cache_2(
                                                    name="{name}.png".format(name=name)
                                                ),
                                                dpi=200,
                                                bbox_inches="tight",
                                            )
                                            figure.clf()
                                            plt.close()

                                            log(
                                                "{baseName} stats, {name} figure saved".format(
                                                    name=name,
                                                    baseName=dir_cache_2.path_name,
                                                )
                                            )

                                    fig_pool.apply_async(thread_fn_2)

                                with dir2.enter(name="valid_pairs") as dir3:
                                    dir_cache_3 = dir3

                                    def thread_fn_3():
                                        df = self.valid_pairs_data_frame

                                        if generate_csvs:
                                            df.to_csv(
                                                path_or_buf=dir_cache_3(
                                                    name="{name}.csv.bz2".format(
                                                        name=dir_cache_3.path_name
                                                    )
                                                ),
                                                index=False,
                                            )

                                        process_valid_pairs_per_sequence(
                                            df=df[df.norm == 0],
                                            dir=dir_cache_3,
                                            name="valid_pairs_raw_counts",
                                        )
                                        process_valid_pairs_per_sequence(
                                            df=df[df.norm == 1],
                                            dir=dir_cache_3,
                                            name="valid_pairs_norm_counts",
                                        )
                                        process_valid_pairs_per_sequence(
                                            df=df[df.norm == 0],
                                            dir=dir_cache_3,
                                            name="total_valid_pairs_raw_count",
                                            totals=True,
                                        )
                                        process_valid_pairs_per_sequence(
                                            df=df[df.norm == 1],
                                            dir=dir_cache_3,
                                            name="total_valid_pairs_norm_counts",
                                            totals=True,
                                        )

                                    fig_pool.apply_async(thread_fn_3)

                            with dir.enter(name="restriction_sites") as dir2:
                                df = self.restrictionSiteDistances.data_frame

                                df_cache_4 = df
                                dir_cache_4 = dir2

                                def thread_fn_4():
                                    if generate_csvs:
                                        df_cache_4.to_csv(
                                            path_or_buf=dir_cache_4(
                                                name="{name}.csv.bz2".format(
                                                    name=dir_cache_4.path_name
                                                )
                                            ),
                                            index=False,
                                        )

                                pool.apply_async(thread_fn_4)

                                if generate_figs:
                                    with dir2.enter(name="valid_pairs") as dir3:
                                        dir_cache_5 = dir3

                                        def thread_fn_5():
                                            process_site_distance_plots(
                                                df=df_cache_4[df_cache_4.valid],
                                                dir=dir_cache_5,
                                                name=dir_cache_4.path_name,
                                                fr_classes=("fr", "rf"),
                                                other_classes=("ff", "rr"),
                                            )

                                        fig_pool.apply_async(thread_fn_5)

                                if generate_figs:
                                    with dir2.enter(name="invalid_pairs") as dir3:
                                        dir_cache_6 = dir3

                                        def thread_fn_6():
                                            process_site_distance_plots(
                                                df=df_cache_4[
                                                    df_cache_4.valid == False
                                                ],
                                                dir=dir_cache_6,
                                                name=dir_cache_4.path_name,
                                                fr_classes=("sc", "de"),
                                                other_classes=("sfs", "rl"),
                                            )

                                        fig_pool.apply_async(thread_fn_6)

                            with dir.enter(
                                name="intra_sequence_pair_separations"
                            ) as dir2:
                                df = self.intraSeqReadDistances.data_frame

                                df_cache_7 = df
                                dir_cache_7 = dir2

                                def thread_fn_7():
                                    if generate_csvs:
                                        df_cache_7.to_csv(
                                            path_or_buf=dir_cache_7(
                                                name="{name}.csv.bz2".format(
                                                    name=dir_cache_7.path_name
                                                )
                                            ),
                                            index=False,
                                        )

                                pool.apply_async(thread_fn_7)

                                (
                                    df_fits,
                                    df_f_tests,
                                ) = (
                                    self.intraSeqReadDistances.calculate_histograms_and_fit_models()
                                )
                                with dir2.enter(name="histograms_plus_models") as dir3:

                                    df_cache_8 = df_fits
                                    dir_cache_8 = dir3

                                    def thread_fn_8():
                                        df = df_cache_8

                                        if generate_csvs:
                                            df.to_csv(
                                                path_or_buf=dir_cache_8(
                                                    name="{name}.csv.bz2".format(
                                                        name=dir_cache_8.path_name
                                                    )
                                                ),
                                                index=False,
                                            )

                                        if generate_figs:
                                            plot = sns.relplot(
                                                x="sep",
                                                y="prob",
                                                kind="line",
                                                data=df[
                                                    (df["class"] == "all")
                                                    & (df["type"] == "histogram")
                                                ],
                                            )
                                            plot.set(xscale="log", yscale="log")
                                            for ax in plot.fig.axes:
                                                ax.set_xlabel("pair separation, bp")
                                                ax.set_ylabel("probability")

                                            plt.tight_layout()
                                            plot.fig.savefig(
                                                dir_cache_8(
                                                    name="{name}.png".format(
                                                        name="pair_separation"
                                                    )
                                                ),
                                                dpi=200,
                                                bbox_inches="tight",
                                            )
                                            plt.close()

                                            plot = sns.relplot(
                                                x="sep",
                                                y="prob",
                                                kind="line",
                                                data=df[
                                                    (df["class"] != "all")
                                                    & (df["type"] == "histogram")
                                                ],
                                                hue="class",
                                            )
                                            plot.set(xscale="log", yscale="log")
                                            for ax in plot.fig.axes:
                                                ax.set_xlabel("pair separation, bp")
                                                ax.set_ylabel("probability")

                                            plt.tight_layout()
                                            plot.fig.savefig(
                                                dir_cache_8(
                                                    name="{name}.png".format(
                                                        name="pair_separation_by_type"
                                                    )
                                                ),
                                                dpi=200,
                                                bbox_inches="tight",
                                            )
                                            plt.close()

                                            plot = sns.relplot(
                                                x="sep",
                                                y="prob",
                                                kind="line",
                                                data=df[df["class"] != "all"],
                                                hue="type",
                                                legend="full",
                                                col="class",
                                                col_wrap=2,
                                                col_order=["ff", "fr", "rf", "rr"],
                                            )
                                            plot.set(xscale="log", yscale="log")
                                            for ax in plot.fig.axes:
                                                ax.set_xlabel("pair separation, bp")
                                                ax.set_ylabel("probability")

                                            plt.tight_layout()
                                            plot.fig.savefig(
                                                dir_cache_8(
                                                    name="{name}.png".format(
                                                        name="model_fits_by_pair_type"
                                                    )
                                                ),
                                                dpi=200,
                                                bbox_inches="tight",
                                            )
                                            plt.close()

                                            plot = sns.relplot(
                                                x="sep",
                                                y="prob",
                                                kind="line",
                                                data=df[df["class"] == "all"],
                                                hue="type",
                                                legend="full",
                                            )
                                            plot.set(xscale="log", yscale="log")
                                            for ax in plot.fig.axes:
                                                ax.set_xlabel("pair separation, bp")
                                                ax.set_ylabel("probability")

                                            plt.tight_layout()
                                            plot.fig.savefig(
                                                dir_cache_8(
                                                    name="{name}.png".format(
                                                        name="model_fits"
                                                    )
                                                ),
                                                dpi=200,
                                                bbox_inches="tight",
                                            )
                                            plt.close()

                                        log(
                                            "{name} stats saved".format(
                                                name=dir_cache_8.path_name
                                            )
                                        )

                                    fig_pool.apply_async(thread_fn_8)

                                with dir2.enter(name="model_F_tests") as dir3:

                                    dir_cache_9 = dir3
                                    df_cache_9 = df_f_tests

                                    def thread_fn_9():
                                        df = df_cache_9

                                        if generate_csvs:
                                            df.to_csv(
                                                path_or_buf=dir_cache_9(
                                                    name="{name}.csv.bz2".format(
                                                        name="model_F_tests"
                                                    )
                                                ),
                                                index=False,
                                            )

                                        process_model_f_tests(
                                            df=df,
                                            dir=dir_cache_9,
                                            name="F_tests_by_pair_type",
                                            all=False,
                                        )
                                        process_model_f_tests(
                                            df=df,
                                            dir=dir_cache_9,
                                            name="F_tests",
                                            all=True,
                                        )

                                        log(
                                            "{name} stats saved".format(
                                                name=dir_cache_9.path_name
                                            )
                                        )

                                    fig_pool.apply_async(thread_fn_9)

                            with dir.enter(
                                name="intra_sequence_fragment_separations"
                            ) as dir2:
                                df = self.intraSeqFragmentDistances.data_frame

                                df_cache_10 = df
                                dir_cache_10 = dir2

                                def thread_fn_10():
                                    if generate_csvs:
                                        df_cache_10.to_csv(
                                            path_or_buf=dir_cache_10(
                                                name="{name}.csv.bz2".format(
                                                    name=dir_cache_10.path_name
                                                )
                                            ),
                                            index=False,
                                        )

                                pool.apply_async(thread_fn_10)

                                (
                                    df_fits2,
                                    df_f_tests2,
                                ) = (
                                    self.intraSeqFragmentDistances.calculate_histograms_and_fit_models()
                                )
                                with dir2.enter(name="histograms_plus_models") as dir3:

                                    df_cache_11 = df_fits2
                                    dir_cache_11 = dir3

                                    def thread_fn_11():
                                        df = df_cache_11

                                        if generate_csvs:
                                            df.to_csv(
                                                path_or_buf=dir_cache_11(
                                                    name="{name}.csv.bz2".format(
                                                        name=dir_cache_11.path_name
                                                    )
                                                ),
                                                index=False,
                                            )

                                        if generate_figs:
                                            plot = sns.relplot(
                                                x="sep",
                                                y="prob",
                                                kind="line",
                                                data=df[
                                                    (df["class"] == "all")
                                                    & (df["type"] == "histogram")
                                                ],
                                            )
                                            plot.set(xscale="log", yscale="log")
                                            for ax in plot.fig.axes:
                                                ax.set_xlabel("fragment separation, bp")
                                                ax.set_ylabel("probability")

                                            plt.tight_layout()
                                            plot.fig.savefig(
                                                dir_cache_11(
                                                    name="{name}.png".format(
                                                        name="pair_separation"
                                                    )
                                                ),
                                                dpi=200,
                                                bbox_inches="tight",
                                            )
                                            plt.close()

                                            plot = sns.relplot(
                                                x="sep",
                                                y="prob",
                                                kind="line",
                                                data=df[
                                                    (df["class"] != "all")
                                                    & (df["type"] == "histogram")
                                                ],
                                                hue="class",
                                            )
                                            plot.set(xscale="log", yscale="log")
                                            for ax in plot.fig.axes:
                                                ax.set_xlabel("fragment separation, bp")
                                                ax.set_ylabel("probability")

                                            plt.tight_layout()
                                            plot.fig.savefig(
                                                dir_cache_11(
                                                    name="{name}.png".format(
                                                        name="pair_separation_by_type"
                                                    )
                                                ),
                                                dpi=200,
                                                bbox_inches="tight",
                                            )
                                            plt.close()

                                            plot = sns.relplot(
                                                x="sep",
                                                y="prob",
                                                kind="line",
                                                data=df[df["class"] != "all"],
                                                hue="type",
                                                legend="full",
                                                col="class",
                                                col_wrap=2,
                                                col_order=["ff", "fr", "rf", "rr"],
                                            )
                                            plot.set(xscale="log", yscale="log")
                                            for ax in plot.fig.axes:
                                                ax.set_xlabel("fragment separation, bp")
                                                ax.set_ylabel("probability")

                                            plt.tight_layout()
                                            plot.fig.savefig(
                                                dir_cache_11(
                                                    name="{name}.png".format(
                                                        name="model_fits_by_pair_type"
                                                    )
                                                ),
                                                dpi=200,
                                                bbox_inches="tight",
                                            )
                                            plt.close()

                                            plot = sns.relplot(
                                                x="sep",
                                                y="prob",
                                                kind="line",
                                                data=df[df["class"] == "all"],
                                                hue="type",
                                                legend="full",
                                            )
                                            plot.set(xscale="log", yscale="log")
                                            for ax in plot.fig.axes:
                                                ax.set_xlabel("fragment separation, bp")
                                                ax.set_ylabel("probability")

                                            plt.tight_layout()
                                            plot.fig.savefig(
                                                dir_cache_11(
                                                    name="{name}.png".format(
                                                        name="model_fits"
                                                    )
                                                ),
                                                dpi=200,
                                                bbox_inches="tight",
                                            )
                                            plt.close()

                                        log(
                                            "{name} stats saved".format(
                                                name=dir_cache_11.path_name
                                            )
                                        )

                                    fig_pool.apply_async(thread_fn_11)

                                with dir2.enter(name="model_F_tests") as dir3:

                                    dir_cache_12 = dir3
                                    df_cache_12 = df_f_tests2

                                    def thread_fn_12():
                                        df = df_cache_12

                                        if generate_csvs:
                                            df.to_csv(
                                                path_or_buf=dir_cache_12(
                                                    name="{name}.csv.bz2".format(
                                                        name="model_F_tests"
                                                    )
                                                ),
                                                index=False,
                                            )

                                        process_model_f_tests(
                                            df=df,
                                            dir=dir_cache_12,
                                            name="F_tests_by_pair_type",
                                            all=False,
                                        )
                                        process_model_f_tests(
                                            df=df,
                                            dir=dir_cache_12,
                                            name="F_tests",
                                            all=True,
                                        )

                                        log(
                                            "{name} stats saved".format(
                                                name=dir_cache_12.path_name
                                            )
                                        )

                                    fig_pool.apply_async(thread_fn_12)

                        fig_pool.close()
                        fig_pool.join()
                        pool.close()
                        pool.join()

                except Exception as ex:
                    self.logger.error("Error saving stats: {msg}".format(msg=ex))
                    raise ex

    class SamReader(object):
        """Class for reading in an existing / external SAM alignment."""

        def __init__(self, name, mark_dups):
            """

            :param name: (str) Path to SAM/BAM/CRAM alignment. May be "<stdin>".
            :param mark_dups: (bool) Whether or not to run the samtools mark_dup pipeline on the input.
            """
            self.name = name
            self.markDups = mark_dups

        def __str__(self):
            return "SAM_input_file: {file}, run mark_dups: {mkdup}".format(
                file=self.name, mkdup="yes" if self.markDups else "no"
            )

    def register_sam_input(self, file_name, mark_dups=True):
        """
        Helper function to assign an existing / external SAM alignment as an input.

        :param file_name: (str) Path to SAM/BAM/CRAM alignment. May be "<stdin>".
        :param mark_dups: (bool) Whether or not to run the samtools mark_dup pipeline on the input.
        """
        self.sam_input = self.SamReader(name=file_name, mark_dups=mark_dups)

    class Align(object):
        """Class for performing an alignment of HiC reads."""

        def __init__(self, reads1, reads2, mark_dups, trim, version):
            """

            :param reads1: (str) Path to first read set.
            :param reads2: (str) Path to second read set. May be None, in which case reads1 will be assumed to be interleaved reads.
            :param mark_dups: (bool) Whether or not to run the samtools mark_dup pipeline on the resulting alignments.
            :param trim: (bool) Run HiC read trimming.
            :param version: (int, = 1, 2 or 3) Use bwa mem (1), bwa-mem2 (2) or minimap2 (3).

            One and only one of reads1 or reads2 may be "<stdin>".
            """
            self.reads1 = reads1
            self.reads2 = reads2
            self.markDups = mark_dups
            self.trim = trim
            self.version = version

        def __str__(self):
            return "Alignment: {version} {trim} read-trimming. {n_reads} fastq read file{s}: {reads}, run mark_dups: {mkdup}".format(
                reads=self.reads1
                if self.reads2 is None
                else "{r1},{r2}".format(r1=self.reads1, r2=self.reads2),
                n_reads=1 if self.reads2 is None else 2,
                s="" if self.reads2 is None else "s",
                mkdup="yes" if self.markDups else "no",
                version="bwa mem"
                if self.version == 1
                else ("bwa-mem2" if self.version == 2 else "minimap2"),
                trim="with" if self.trim else "without",
            )

    class SamReads(Align):
        """Class for performing an alignment of HiC reads in SAM/BAM/CRAM format."""

        def __init__(self, reads, mark_dups, trim, version, tags=[]):
            """

            :param reads: (str) Path to SAM/BAM/CRAM reads. May be "<stdin>".
            :param mark_dups: (bool) Whether or not to run the samtools mark_dup pipeline on the resulting alignments.
            :param tags: (list(str)) List of SAM tags to append to reads.
            :param trim: (bool) Run HiC read trimming.
            :param version: (int, = 1, 2 or 3) Use bwa mem (1), bwa-mem2 (2) or minimap2 (3).
            """
            super().__init__(
                reads1=reads,
                reads2=None,
                mark_dups=mark_dups,
                trim=trim,
                version=version,
            )
            self.tags = tags

        def __str__(self):
            return "Alignment: {version} {trim} read-trimming. SAM input reads: {reads}, run mark_dups: {mkdup}{extra}".format(
                reads=self.reads1,
                mkdup="yes" if self.markDups else "no",
                extra=""
                if len(self.tags) == 0
                else ", extra appended SAM tags: {tags}".format(
                    tags=",".join(self.tags)
                ),
                version="bwa mem"
                if self.version == 1
                else ("bwa-mem2" if self.version == 2 else "minimap2"),
                trim="with" if self.trim else "without",
            )

    class IndexOnly(Align):
        """Class for only performing an alignment reference index"""

        def __init__(self, trim, version):
            """

            :param trim: (bool) Run HiC read trimming.
            :param version: (int, = 1, 2 or 3) Use bwa mem (1), bwa-mem2 (2) or minimap2 (3).
            """
            super().__init__(
                reads1=None, reads2=None, mark_dups=None, trim=trim, version=version
            )

    def register_alignment(self, reads, mark_dups=True, trim=True, version=2):
        """
        Helper function to perform an alignment of HiC reads in (gzipped) FASTQ format.

        :param reads: (list(str)) One- or two-item list of paths to reads. One read path indicates reads are interleaved. One and only one path may be "<stdin>".
        :param mark_dups: (bool) Whether or not to run the samtools mark_dup pipeline on the resulting alignments.
        :param trim: (bool) Run HiC read-trimming, trim sections of reads that align past restriction sites.
        :param version: (int, = 1, 2 or 3) Use bwa mem (1), bwa-mem2 (2) or minimap2 (3).
        """

        if version not in (1, 2, 3):
            self.log_error(
                exception=self.Exception(
                    "Alignment: 'version' needs to equal either 1, 2 or 3"
                )
            )
            return

        if type(reads) is list and 0 < len(reads) < 3:
            reads1 = reads[0]
            reads2 = reads[1] if len(reads) == 2 else None
            self.sam_input = self.Align(
                reads1=reads1,
                reads2=reads2,
                mark_dups=mark_dups,
                trim=trim,
                version=version,
            )
        else:
            self.log_error(
                exception=self.Exception(
                    "Alignment: 'reads' needs to be a list of one or two files"
                )
            )

    def register_alignment_sam_reads(
        self, reads, tags=[], mark_dups=True, trim=True, version=2
    ):
        """
        Helper function to perform an alignment of HiC reads in SAM/BAM/CRAM format.

        :param reads: (str) Path to reads. May be "<stdin>".
        :param tags: (list(str)) List of SAM tags to append to reads
        :param mark_dups: (bool) Whether or not to run the samtools mark_dup pipeline on the resulting alignments.
        :param trim: (bool) Run HiC read-trimming, trim sections of reads that align past restriction sites.
        :param version: (int, = 1, 2 or 3) Use bwa mem (1), bwa-mem2 (2) or minimap2 (3).
        """

        if version not in (1, 2, 3):
            self.log_error(
                exception=self.Exception(
                    "Alignment: 'version' needs to equal either 1, 2 or 3"
                )
            )
            return

        self.sam_input = self.SamReads(
            reads=reads, mark_dups=mark_dups, tags=tags, trim=trim, version=version
        )

    def register_alignment_index_only(self, trim=True, version=2):
        """
        Helper function to perform and save only an alignment index.

        :param trim: (bool) Run HiC read-trimming, trim sections of reads that align past restriction sites.
        :param version: (int, = 1, 2 or 3) Use bwa mem (1), bwa-mem2 (2) or minimap2 (3).
        """

        if version not in (1, 2, 3):
            self.log_error(
                exception=self.Exception(
                    "Alignment: 'version' needs to equal either 1, 2 or 3"
                )
            )
            return

        self.sam_input = self.IndexOnly(trim=trim, version=version)


def documenter(docstring):
    def inner_documenter(f):
        f.__doc__ = docstring
        return f

    return inner_documenter


@click.group(chain=True)
@documenter(
    """{name} version {ver}. {des}
    {licence}
    
    \b
    Usage Example
    -------------
    HiLine params -t 32 Data/reference.fa.bgz Arima align-sam-reads Data/HiC_reads.cram valid-pairs Results/valid_pairs.cram save-stats Results/stats
    
    \b
    Required Commands
    -----------------
        params:                         set pipeline parameters
    
    \b
    Input Commands (one and only-one must be invoked)
    --------------
        read-sam:                       read in external alignment
        align-one-read:                 alignment with one interleaved (gzipped) FASTQ read source
        align-two-reads:                alignment with two (gzipped) FASTQ read sources
        align-sam-reads:                alignment with a SAM/BAM/CRAM read source
        index-only:                     just create and save a reference index
    
    \b    
    Output Commands
    ---------------
        all-reads:                      alias for good-reads + bad-reads
        good-reads:                     alias for valid-pairs + invalid-reads
        valid-pairs:                    alias for valid-ff + valid-fr + valid-rf + valid-rr        
        valid-ff:                       valid ff read pairs
        valid-fr:                       valid fr read pairs
        valid-rf:                       valid rf read pairs
        valid-rr:                       valid rr read pairs
        invalid-reads:                  alias for invalid-pairs + dumped
        invalid-pairs:                  alias for self-circle + dangling-end + same-frag-and-strand + re-ligated
        self-circle:                    self-circular read pairs
        dangling-end:                   dangling-end read pairs
        same-frag-and-strand:           read pairs on the same restriction fragment and strand
        re-ligated:                     re-ligated read pairs
        dumped:                         dumped reads (good reads with a bad mate read)
        bad-reads:                      alias for low-mapq + too-far-from-restriction-site + bad-reference + unmapped + unpaired + supplementary + qc-fail + duplicate + secondary
        low-mapq:                       reads below minimum mapping quality
        too-far-from-restriction-site:  reads too far from a restriction site
        bad-reference:                  reads aligned to an invalid reference
        unmapped:                       unmapped reads
        unpaired:                       unpaired reads
        supplementary:                  supplementary reads
        qc-fail:                        qc-failed reads
        duplicate:                      duplicate reads
        secondary:                      secondary reads
    
    \b    
    Saving Statistics
    -----------------    
        save-stats:                     save alignment statistics to disk
    """.format(
        name=NAME, ver=VERSION, des=DESCRIPTION, licence=LICENCE
    )
)
def cli():
    logger.setLevel(logging.INFO)
    logger.addHandler(
        create_logger_handle(stream=sys.stderr, typeid="status", level=logging.INFO)
    )
    logger.addHandler(
        create_logger_handle(stream=sys.stderr, typeid="error", level=logging.ERROR)
    )
    logger.addHandler(
        create_logger_handle(stream=sys.stderr, typeid="warning", level=logging.WARNING)
    )

    def _showwarning(message, category, filename, lineno, file=None, line=None):
        logger.warning(
            f"[{filename} {lineno}] {message}"
            if line is None
            else f"[{filename} {lineno}] {message} {line}"
        )

    warnings.showwarning = _showwarning


def processor(f):
    def new_func(*args, **kwargs):
        def processor(pipeline):
            return f(pipeline, *args, **kwargs)

        return processor

    return update_wrapper(new_func, f)


@cli.command()
@click.argument("reference", type=click.Path(exists=True))
@click.argument("restriction_sites")
@click.option(
    "-t",
    "--threads",
    type=click.IntRange(3, None, clamp=True),
    default=4,
    help="Number of threads to use, must be at least 3. Default=4",
)
@click.option(
    "-q",
    "--minmapq",
    type=click.IntRange(0, None, clamp=True),
    default=10,
    help="Minimum mapping quality. Default=10",
)
@processor
@documenter(
    """
Sets pipeline parameters.

\b
REFERENCE:
    Path to a reference genome in (gzipped) FASTA format.
    
\b
RESTRICTION_SITES:
    Restriction site(s) specification string, a comma-separated list of restriction site definitions.
    /^<restriction>(, <restriction>)*$/
    
    where <restriction> is either: the name of a HiC kit provider, the name of a restriction enzyme, or the string-definition of the restriction site.
    <restriction> = /(<kit>|<enzyme>|<site>)/
    
    <kit> is the name of a HiC kit provider, currently supported are:
        Arima Genomics version 2: "Arima_v2"
        Arima Genomics: "Arima"
        Dovetail Genomics: "Dovetail"
        Phase Genomics: "Phase"
        Qiagen: "Qiagen"
    <kit> = /((Arima_v2)|(Arima)|(Dovetail)|(Phase)|(Qiagen))/, note: names are case-insensitive.
    
    <enzyme> is the name of a restriction enzyme e.g. "DpnII", currently all known enzymes defined in the 'Biopython.Restriction' package are supported.
    Note: enzyme names are case-sensitive.
    
    <site> is the IUPAC string definition of a restriction site e.g. "^GATC".
    <site> = /([ACGTWSMKRYBDHVN]*\^[ACGTWSMKRYBDHVN]*_?[ACGTWSMKRYBDHVN]*)/
    A caret character "^" is required to appear once and only once, and defines the cut location.
    If the site is not palindromic, an underscore character "_" is required to appear once and only once, and defines the cut location on the reverse strand. Otherwise it is not required.
"""
)
def params(pipeline, reference, restriction_sites, threads, minmapq):
    pipeline.reference = reference
    pipeline.restriction_sites = restriction_sites
    pipeline.threads = threads
    pipeline.min_mapq = minmapq


@cli.command()
@click.argument("sam", type=click.File("rb", lazy=False))
@click.option(
    "--rmdups/--no-rmdups",
    default=True,
    help="Run samtools mark_dup pipeline on alignment. Default=rmdups",
)
@processor
@documenter(
    """
Read in an external alignment.

\b
SAM:
    Path to alignment in SAM/BAM/CRAM format. Use "-" for stdin.
"""
)
def read_sam(pipeline, sam, rmdups):
    pipeline.register_sam_input(file_name=sam.name, mark_dups=rmdups)


@cli.command()
@click.argument("reads", type=click.File("rb", lazy=False))
@click.option(
    "--rmdups/--no-rmdups",
    default=True,
    help="Run samtools mark_dup pipeline on alignment. Default=rmdups",
)
@click.option(
    "--trim/--no-trim",
    default=True,
    help="Run HiC read trimming, trim sections of reads that align past restriction sites. Default=trim",
)
@click.option("--bwa1", "version", flag_value=1, help="Use bwa mem. Default=False")
@click.option(
    "--bwa2", "version", flag_value=2, default=True, help="Use bwa-mem2. Default=True"
)
@click.option(
    "--minimap2",
    "version",
    flag_value=3,
    default=False,
    help="Use minimap2. Default=False",
)
@processor
@documenter(
    """
Alignment with one read source.

\b
READS:
    Path to interleaved reads in (gzipped) FASTQ format. Use "-" for stdin.
"""
)
def align_one_read(pipeline, reads, rmdups, trim, version):
    pipeline.register_alignment(
        reads=[reads.name], mark_dups=rmdups, trim=trim, version=version
    )


@cli.command()
@click.argument("reads", type=click.File("rb", lazy=False), nargs=2)
@click.option(
    "--rmdups/--no-rmdups",
    default=True,
    help="Run samtools mark_dup pipeline on alignment. Default=rmdups",
)
@click.option(
    "--trim/--no-trim",
    default=True,
    help="Run HiC read trimming, trim sections of reads that align past restriction sites. Default=trim",
)
@click.option("--bwa1", "version", flag_value=1, help="Use bwa mem. Default=False")
@click.option(
    "--bwa2", "version", flag_value=2, default=True, help="Use bwa-mem2. Default=True"
)
@click.option(
    "--minimap2",
    "version",
    flag_value=3,
    default=False,
    help="Use minimap2. Default=False",
)
@processor
@documenter(
    """
Alignment with two read sources.

\b
READS:
    Two paths to reads in (gzipped) FASTQ format. Use "-" for stdin (for one and only-one path).
"""
)
def align_two_reads(pipeline, reads, rmdups, trim, version):
    pipeline.register_alignment(
        reads=[r.name for r in reads], mark_dups=rmdups, trim=trim, version=version
    )


@cli.command()
@click.argument("reads", type=click.File("rb", lazy=False))
@click.option(
    "--rmdups/--no-rmdups",
    default=True,
    help="Run samtools mark_dup pipeline on alignment. Default=rmdups",
)
@click.option("--tag", "-t", multiple=True, help="SAM tag(s) to append to reads.")
@click.option(
    "--trim/--no-trim",
    default=True,
    help="Run HiC read trimming, trim sections of reads that align past restriction sites. Default=trim",
)
@click.option("--bwa1", "version", flag_value=1, help="Use bwa mem. Default=False")
@click.option(
    "--bwa2", "version", flag_value=2, default=True, help="Use bwa-mem2. Default=True"
)
@click.option(
    "--minimap2",
    "version",
    flag_value=3,
    default=False,
    help="Use minimap2. Default=False",
)
@processor
@documenter(
    """
Alignment with SAM/BAM/CRAM reads.

\b
READS:
    Path to reads in SAM/BAM/CRAM format. Use "-" for stdin.
"""
)
def align_sam_reads(pipeline, reads, rmdups, tag, trim, version):
    pipeline.register_alignment_sam_reads(
        reads=reads.name,
        mark_dups=rmdups,
        tags=[t for tg in tag for t in tg.split()],
        trim=trim,
        version=version,
    )


@cli.command()
@click.option(
    "--trim/--no-trim",
    default=True,
    help="Run HiC read trimming, trim sections of reads that align past restriction sites. Default=trim",
)
@click.option("--bwa1", "version", flag_value=1, help="Use bwa mem. Default=False")
@click.option(
    "--bwa2", "version", flag_value=2, default=True, help="Use bwa-mem2. Default=True"
)
@click.option(
    "--minimap2",
    "version",
    flag_value=3,
    default=False,
    help="Use minimap2. Default=False",
)
@processor
@documenter(
    """
Alignment index only.
"""
)
def index_only(pipeline, trim, version):
    pipeline.register_alignment_index_only(
        trim=trim,
        version=version,
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write valid-pairs to OUTPUT.

\b
valid-pairs is an alias for (valid-ff, valid-fr, valid-rf, valid-rr).

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def valid_pairs(pipeline, output, sort):
    for handle in pipeline.output.valid_pairs:
        pipeline.register_output_file(file_name=output.name, handle=handle, sort=sort)


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write invalid-pairs to OUTPUT.

\b
invalid-pairs is an alias for (self-circle, dangling-end, same-frag-and-strand, re-ligated).

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def invalid_pairs(pipeline, output, sort):
    for handle in pipeline.output.invalid_pairs:
        pipeline.register_output_file(file_name=output.name, handle=handle, sort=sort)


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write invalid-reads to OUTPUT.

\b
invalid-reads is an alias for (invalid-pairs, dumped).

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def invalid_reads(pipeline, output, sort):
    for handle in pipeline.output.invalid_reads:
        pipeline.register_output_file(file_name=output.name, handle=handle, sort=sort)


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write good-reads to OUTPUT.

\b
good-reads is an alias for (valid-pairs, invalid-reads).

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def good_reads(pipeline, output, sort):
    for handle in pipeline.output.good_reads:
        pipeline.register_output_file(file_name=output.name, handle=handle, sort=sort)


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write bad-reads to OUTPUT.

\b
bad-reads is an alias for (low-mapq, too-far-from-restriction-site, bad-reference, unmapped, unpaired, supplementary, qc-fail, duplicate, secondary).

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def bad_reads(pipeline, output, sort):
    for handle in pipeline.output.bad_reads:
        pipeline.register_output_file(file_name=output.name, handle=handle, sort=sort)


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write all-reads to OUTPUT.

\b
all-reads as an alias for (good-reads, bad-reads).

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def all_reads(pipeline, output, sort):
    for handle in pipeline.output.all_reads:
        pipeline.register_output_file(file_name=output.name, handle=handle, sort=sort)


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write valid-ff reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def valid_ff(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.ff, sort=sort
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write valid-fr reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def valid_fr(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.fr, sort=sort
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write valid-rf reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def valid_rf(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.rf, sort=sort
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write valid-rr reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def valid_rr(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.rr, sort=sort
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write re-ligated reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def re_ligated(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.religated, sort=sort
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write same-frag-and-strand reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def same_frag_and_strand(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.samefragmentandstrand, sort=sort
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write dangling-end reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def dangling_end(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.danglingend, sort=sort
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write self-circle reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def self_circle(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.selfcircle, sort=sort
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write dumped reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def dumped(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.dump, sort=sort
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write too-far-from-restriction-site reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def too_far_from_restriction_site(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name,
        handle=pipeline.output.toofarfromrestrictionsite,
        sort=sort,
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write low-mapq reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def low_mapq(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.belowminmapq, sort=sort
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write bad-reference reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def bad_reference(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.invalidreferencename, sort=sort
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write supplementary reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def supplementary(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.supplementary, sort=sort
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write duplicate reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def duplicate(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.duplicate, sort=sort
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write qc-fail reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def qc_fail(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.qcfail, sort=sort
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write secondary reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def secondary(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.secondary, sort=sort
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write unmapped reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def unmapped(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.unmapped, sort=sort
    )


@cli.command()
@click.argument("output", type=click.File("w"))
@click.option(
    "--sort/--no-sort", default=True, help="Sort and index output. Default=sort"
)
@processor
@documenter(
    """
Write unpaired reads to OUTPUT.

\b
OUTPUT:
    Path to output SAM/BAM/CRAM file.
"""
)
def unpaired(pipeline, output, sort):
    pipeline.register_output_file(
        file_name=output.name, handle=pipeline.output.unpaired, sort=sort
    )


@cli.command()
@click.argument("path", type=click.Path())
@processor
@documenter(
    """
Save alignment statistics to disk.
    
\b
PATH:
    Base path of save location.
"""
)
def save_stats(pipeline, path):
    pipeline.save_stats_path = path


@cli.resultcallback()
def process_args(processors):
    pipeline = Pipeline()
    pipeline.logger = logger

    def exit_on_error(_):
        sys.exit(1)

    pipeline.exception_callback = exit_on_error

    pipeline.command_line = " ".join([basename(normpath(sys.argv[0]))] + sys.argv[1:])

    for processor in processors:
        processor(pipeline)

    pipeline.run()
