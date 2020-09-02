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
import sys
import click
import hashlib
import re
import gzip
import os
from os import makedirs
from os.path import isdir, isfile, join, basename, normpath
from subprocess import Popen, PIPE, STDOUT, check_output
from enum import Enum, auto
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tempfile import NamedTemporaryFile
from binascii import hexlify
from threading import Thread

from HiLine import version
from HiLine import Pipeline as pl

from _Aligner import _Aligner_Main

NAME = "_HiLine_Aligner"
DESCRIPTION = "Restriction digest aware HiC aligner. Part of HiLine."
LICENCE = "Copyright (c) 2020 Ed Harry, Wellcome Sanger Institute."

logger = logging.getLogger(__name__)


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
        logging.Formatter("[Aligner {id}] :: %(message)s".format(id=typeid))
    )
    handle.addFilter(LogFilter(level=level))

    return handle


class LoggerHandle(object):
    def __init__(self):
        self.threadsAndHandles = []

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

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        for thread, handle in self.threadsAndHandles:
            os.close(handle)
            thread.join()


def documenter(docstring):
    def inner_documenter(f):
        f.__doc__ = docstring
        return f

    return inner_documenter


@click.command()
@documenter(
    """
        {name} {ver}. {des}
        {lic}
        """.format(
        name=NAME, ver=version, des=DESCRIPTION, lic=LICENCE
    )
)
@click.argument("reference", type=click.Path(exists=True))
@click.argument("reads", type=click.File("rb", lazy=False), nargs=-1)
@click.option("--site", type=(str, int, int), multiple=True)
@click.option(
    "-t", "--threads", type=click.IntRange(2, None, clamp=True), default=2,
)
@click.option("--tag", type=str, multiple=True)
@click.option("--trim/--no-trim", default=True)
@click.option("--bwa1", "bwa", flag_value=1)
@click.option("--bwa2", "bwa", flag_value=2, default=True)
def cli(reference, reads, site, threads, tag, trim, bwa):
    logger.setLevel(logging.INFO)
    logger.addHandler(
        create_logger_handle(stream=sys.stderr, typeid="status", level=logging.INFO)
    )
    logger.addHandler(
        create_logger_handle(stream=sys.stderr, typeid="error", level=logging.ERROR)
    )

    if len(site) == 0:
        raise click.BadOptionUsage("site", "must specify at least one 'site'")

    if not (3 > len(reads) > 0):
        raise click.BadArgumentUsage("only one or two read sources accepted")

    class Aligner(object):
        class Site(object):
            class Overhang(Enum):
                FIVE = auto()
                THREE = auto()
                BLUNT = auto()

                @classmethod
                def get(cls, fwd, rev):
                    if fwd < rev:
                        return cls.FIVE
                    elif rev < fwd:
                        return cls.THREE
                    else:
                        return cls.BLUNT

            def __init__(self, site, fwd, rev):
                self.site = site
                self.fwd = fwd
                self.rev = rev

                if fwd < 0 or rev < 0:
                    raise click.BadOptionUsage(
                        "site fwd and rev locations must be non-negative integers"
                    )

                self.cut_start = min(fwd, rev)
                self.cut_end = max(fwd, rev)

                self.overhang = self.Overhang.get(fwd, rev)

            def __str__(self):
                return ",".join((self.site, str(self.fwd), str(self.rev)))

        def __init__(self, bwa, trim, reference, reads, sites, threads, sam_tags):
            self.trim = trim
            self.reference = reference
            self.reads_from_stdin = "<stdin>" in [read.name for read in reads]
            self.sites = [self.Site(*site) for site in sites]
            self.sam_tags = sam_tags
            self.threads = threads

            self.name = self.reference + "." + str(self.hash) + ".HiLine_Reference"
            self.bwa = self.get_bwa(
                threads=threads // 2,
                reads=["-" if read.name == "<stdin>" else read.name for read in reads],
                version=bwa,
            )

        @staticmethod
        def get_bwa(threads, reads, version):
            class Bwa(object):
                def index_command(self, fasta, prefix=None):
                    if prefix is None:
                        prefix = basename(normpath(fasta))
                    prefix += "." + self.program_name
                    return self.program_name + " index -p " + prefix + " " + fasta

                def run_command(self, reference, stdin=False):
                    return (
                        self.program_name
                        + " mem -5SPYC{p} -t {threads} -B 8 {ref} {r1}{r2}".format(
                            p="p" if stdin else ("p" if len(reads) == 1 else ""),
                            threads=threads,
                            ref=reference + "." + self.program_name,
                            r1="-" if stdin else reads[0],
                            r2=""
                            if stdin
                            else (" " + reads[1] if len(reads) == 2 else ""),
                        )
                    )

            class Bwa1(Bwa):
                program_name = "bwa"
                index_extensions = ("amb", "ann", "bwt", "pac", "sa")

            class Bwa2(Bwa):
                program_name = "bwa-mem2"
                index_extensions = (
                    "0123",
                    "amb",
                    "ann",
                    "bwt.2bit.64",
                    "bwt.8bit.32",
                    "pac",
                )

            try:
                if version == 1:
                    raise FileNotFoundError("bwa 1 selected")

                with Popen(
                    "bwa-mem2 version".split(), stdout=PIPE, stderr=STDOUT
                ) as process:
                    output = "".join(
                        line.decode("utf-8") for line in process.stdout.readlines()
                    )

                match = re.search(
                    r"(?P<major>\d+)\.(?P<minor1>\d+)\.?(?P<minor2>\d*)", output
                )

                if match is None:
                    raise FileNotFoundError("Could not determine 'bwa-mem2' version")

                major = int(match.group("major"))

                if not (major >= 2):
                    raise FileNotFoundError(
                        "Only 'bwa-mem2' version 2 or higher supported"
                    )

                logger.info("Using bwa-mem2")
                return Bwa2()

            except FileNotFoundError as ex:
                logger.info(str(ex))

                try:
                    with Popen("bwa", stdout=PIPE, stderr=STDOUT) as process:
                        output = "".join(
                            [
                                line.decode("utf-8")
                                for line in process.stdout.readlines()
                            ]
                        )

                    match = re.search(
                        r"Version: (?P<major>\d+)\.(?P<minor1>\d+)\.(?P<minor2>\d+)",
                        output,
                    )

                    if match is None:
                        raise Exception("Could not determine 'bwa' version")

                    major = int(match.group("major"))
                    minor1 = int(match.group("minor1"))
                    minor2 = int(match.group("minor2"))

                    if not (
                        major > 0
                        or (major == 0 and minor1 > 7)
                        or (major == 0 and minor1 == 7 and minor2 >= 17)
                    ):
                        raise Exception(
                            "'bwa' version {major}.{minor1}.{minor2} found, version 0.7.17 or greater required".format(
                                major=major, minor1=minor1, minor2=minor2
                            )
                        )

                    logger.info("Using bwa")
                    return Bwa1()

                except Exception as ex:
                    sys.exit(str(ex))

        def digest(self, record):
            cuts = [
                (0, self.Site.Overhang.BLUNT, ""),
                (len(record), self.Site.Overhang.BLUNT, ""),
            ]
            for site in self.sites:
                for match in re.finditer(site.site, str(record.seq)):
                    cuts.append(
                        (
                            match.start() + site.fwd,
                            site.overhang,
                            match.group()[site.cut_start : site.cut_end],
                        )
                    )

            cuts.sort(key=lambda cut: cut[0])

            for cut_start, cut_end in zip(cuts[:-1], cuts[1:]):
                overhang_start = (
                    cut_start[2] if cut_start[1] == self.Site.Overhang.THREE else ""
                )
                overhang_end = (
                    cut_end[2] if cut_end[1] == self.Site.Overhang.FIVE else ""
                )

                yield SeqRecord(
                    overhang_start
                    + record.seq[cut_start[0] : cut_end[0]]
                    + overhang_end,
                    id=record.id
                    + "_"
                    + str(cut_start[0])
                    + "_"
                    + str(len(cut_start[2]))
                    + ("+" if cut_start[1] == self.Site.Overhang.THREE else "-")
                    + "_"
                    + str(len(cut_end[2]))
                    + ("+" if cut_end[1] == self.Site.Overhang.FIVE else "-"),
                    description="",
                    name="",
                )

            logger.info(
                "{id} digested into {n} fragments".format(id=record.id, n=len(cuts) - 1)
            )

        def index(self, logger_handle):
            threads = []

            if self.trim and not (
                isdir(self.name)
                and all(
                    isfile(join(self.name, "ref." + self.bwa.program_name + "." + ext))
                    for ext in self.bwa.index_extensions
                )
            ):
                logger.info("Creating virtual digestion index...")

                def thread_fn():
                    def get_open(reference):
                        with open(reference, "rb") as file:
                            return (
                                gzip.open(reference, "rt")
                                if hexlify(file.read(2)) == b"1f8b"
                                else open(reference)
                            )

                    makedirs(self.name, exist_ok=True)
                    with NamedTemporaryFile(
                        mode="w", buffering=True, suffix="", prefix="ref", dir=self.name
                    ) as tmp_file:
                        with get_open(self.reference) as ref_file:
                            for rec in SeqIO.parse(ref_file, "fasta"):
                                for digest_rec in self.digest(rec):
                                    SeqIO.write(digest_rec, tmp_file, "fasta")
                        handle = logger_handle.add_logger(
                            logger.info, "index digestion fragments"
                        )
                        with Popen(
                            self.bwa.index_command(
                                prefix="ref", fasta=tmp_file.name
                            ).split(),
                            cwd=self.name,
                            stderr=handle,
                            stdout=handle,
                        ) as process:
                            process.communicate()

                thread = Thread(target=thread_fn)
                thread.start()
                threads.append(thread)

            if not (
                all(
                    isfile(self.reference + "." + self.bwa.program_name + "." + ext)
                    for ext in self.bwa.index_extensions
                )
            ):
                logger.info("Creating reference index...")

                def thread_fn():
                    handle = logger_handle.add_logger(logger.info, "index reference")
                    with Popen(
                        self.bwa.index_command(fasta=self.reference).split(),
                        stderr=handle,
                        stdout=handle,
                    ) as process:
                        process.communicate()

                thread = Thread(target=thread_fn)
                thread.start()
                threads.append(thread)

            for thread in threads:
                thread.join()

        def run(self):
            try:
                try:
                    match = re.match(
                        r"^samtools (?P<major>\d+)\.(?P<minor>\d+)",
                        check_output(
                            "samtools --version".split(), stderr=STDOUT
                        ).decode("utf-8"),
                    )
                    if match is None:
                        raise Exception("Could not determine 'samtools' version")

                    major = int(match.group("major"))
                    minor = int(match.group("minor"))

                    if not (major > 1 or (major == 1 and minor >= 10)):
                        raise Exception(
                            "'samtools' version {major}.{minor} found, version 1.10 or greater required".format(
                                major=major, minor=minor
                            )
                        )

                    samtools_version = "{major}.{minor}".format(
                        major=major, minor=minor
                    )

                except FileNotFoundError:
                    raise Exception("'samtools' not found on $PATH")

            except Exception as ex:
                sys.exit(str(ex))

            with LoggerHandle() as logger_handle:

                self.index(logger_handle)

                if self.trim:
                    handle = logger_handle.add_logger(
                        logger.info, "digestion fragment alignment"
                    )
                    with Popen(
                        "samtools view -@ {threads} -hF 0x900 -".format(
                            threads=self.threads
                        ).split(),
                        stdout=PIPE,
                        stderr=handle,
                        stdin=Popen(
                            self.bwa.run_command(
                                reference=join(self.name, "ref")
                            ).split(),
                            stderr=handle,
                            stdout=PIPE,
                            stdin=sys.stdin if self.reads_from_stdin else None,
                        ).stdout,
                    ) as view_process:

                        read_, write_ = os.pipe()

                        def main_fn():
                            class RunParams(object):
                                samInput = view_process.stdout
                                samOutput = write_
                                nThreads = self.threads
                                info = logger_handle.add_logger(
                                    logger.info, "read trimming"
                                )
                                error = logger_handle.add_logger(
                                    logger.error, "read trimming"
                                )

                            _Aligner_Main(params=RunParams())
                            os.close(write_)

                        thread = Thread(target=main_fn)
                        thread.start()
                        threads = [thread]

                        read, write = os.pipe()
                        global seen_header
                        seen_header = False
                        global pg_lines
                        pg_lines = []

                        def thread_fn():
                            global seen_header
                            global pg_lines

                            with os.fdopen(write, "wb") as f_out, os.fdopen(
                                read_, "rb"
                            ) as f_in:
                                for line in f_in:
                                    if not seen_header:
                                        if (
                                            decoded_line := line.decode(
                                                "utf-8", errors="replace"
                                            )
                                        ).startswith("@"):
                                            if decoded_line.startswith("@PG"):
                                                pg_lines.append(line)
                                        else:
                                            seen_header = True

                                    f_out.write(line)

                        thread = Thread(target=thread_fn)
                        thread.start()
                        threads.append(thread)

                        samtools_fastq_cmd = "samtools fastq -@ {threads} -t{tags} -0 /dev/null -F 0x900 -".format(
                            tags=""
                            if len(self.sam_tags) == 0
                            else "T {tg}".format(tg=",".join(t for t in self.sam_tags)),
                            threads=self.threads,
                        )

                        handle = logger_handle.add_logger(
                            logger.info, "reference alignment"
                        )
                        with Popen(
                            self.bwa.run_command(
                                reference=self.reference, stdin=True
                            ).split(),
                            stdout=PIPE,
                            stderr=handle,
                            stdin=Popen(
                                samtools_fastq_cmd.split(),
                                stdout=PIPE,
                                stderr=handle,
                                stdin=read,
                            ).stdout,
                        ) as bwa_process:

                            def thread_fn_2():
                                global seen_header
                                global pg_lines
                                local_pg_lines = []
                                header_buffer = []
                                header_mode = True
                                for line in bwa_process.stdout:
                                    if header_mode:
                                        if (
                                            decoded_line := line.decode(
                                                "utf-8", errors="replace"
                                            )
                                        ).startswith("@"):
                                            if decoded_line.startswith("@PG"):
                                                local_pg_lines.append(line)
                                            else:
                                                header_buffer.append(line)
                                        else:
                                            header_mode = False
                                            while not seen_header:
                                                pass

                                            chain = pl.SamPGChain(pg_lines)
                                            chain.append(
                                                "@PG\tID:{id}\tPN:{pn}\tCL:{cl}\tDS:{ds}\tVN:{vn}\n".format(
                                                    id=NAME,
                                                    pn=NAME,
                                                    cl=" ".join(
                                                        [
                                                            basename(
                                                                normpath(sys.argv[0])
                                                            )
                                                        ]
                                                        + sys.argv[1:]
                                                    ),
                                                    vn=version,
                                                    ds=DESCRIPTION,
                                                ).encode(
                                                    "utf-8"
                                                )
                                            )
                                            chain.append(
                                                "@PG\tID:samtools\tPN:samtools\tVN:{version}\tCL:{cmd}\n".format(
                                                    version=samtools_version,
                                                    cmd=samtools_fastq_cmd,
                                                ).encode(
                                                    "utf-8"
                                                )
                                            )
                                            for pg_line in local_pg_lines:
                                                chain.append(pg_line)

                                            for header_line in header_buffer:
                                                sys.stdout.buffer.write(header_line)

                                            for pg_line in chain:
                                                sys.stdout.buffer.write(pg_line)

                                            sys.stdout.buffer.write(line)
                                    else:
                                        sys.stdout.buffer.write(line)

                            thread_2 = Thread(target=thread_fn_2)
                            thread_2.start()
                            threads.append(thread_2)
                            for thread in threads:
                                thread.join()
                else:
                    with Popen(
                        self.bwa.run_command(reference=self.reference).split(),
                        stderr=logger_handle.add_logger(
                            logger.info, "reference alignment"
                        ),
                        stdout=sys.stdout,
                        stdin=sys.stdin if self.reads_from_stdin else None,
                    ) as process:
                        process.communicate()

        def __str__(self):
            return ";".join(
                (self.name, self.reference, str([str(site) for site in self.sites]),)
            )

        @property
        def hash(self):
            return sum(
                int(hashlib.sha256(string.encode("utf-8")).hexdigest(), 16)
                for string in (
                    [basename(normpath(self.reference))]
                    + [str(site) for site in self.sites]
                )
            ) % (2 ** 64)

    Aligner(
        bwa=bwa,
        trim=trim,
        reference=reference,
        reads=reads,
        sites=site,
        threads=threads,
        sam_tags=[t for tg in tag for t in tg.split()],
    ).run()
