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

import os
import sys
import shutil
import sysconfig
from subprocess import Popen, PIPE, STDOUT

import numpy as np
from setuptools import setup, Extension, find_packages
from setuptools.command.install import install as setup_install

LIBHS = ["libhs.a", "libhs_runtime.a"]
LIBHS_VERSION = "v5.4.0"
LIBHS_SOURCE = "https://github.com/intel/hyperscan.git"


def get_compiler():
    compilers = {"clang": "fast", "gcc": "3"}
    env = os.environ.copy()

    if cc := env.get("CC"):
        if o_flag := compilers.get(cc):
            return cc, o_flag

    for compiler in compilers.keys():
        with Popen(
            ["which", compiler], stdin=PIPE, stdout=PIPE, stderr=PIPE
        ) as process:
            process.communicate()
            if process.returncode == 0:
                return compiler, compilers[compiler]

    raise SystemError("No suitable compiler found")


class HiLineInstall(setup_install):
    def run(self):
        install_hyperscan()
        setup_install.run(self)


def install_hyperscan():
    if not np.all([os.path.exists(lib) for lib in LIBHS]):
        print("Building hyperscan...")

        env = os.environ.copy()
        cc, o_flag = get_compiler()
        env["CC"] = cc

        c_flags = "-fPIC -O{o} -DNDEBUG -lm".format(o=o_flag)
        cpp_flags = c_flags + " -lstdc++"

        with Popen(["which", "git"], stdin=PIPE, stdout=PIPE, stderr=PIPE) as process:
            process.communicate()
            if process.returncode != 0:
                raise SystemError(
                    "Could not find 'git' on your PATH, cannot build hyperscan"
                )

        with Popen(["which", "cmake"], stdin=PIPE, stdout=PIPE, stderr=PIPE) as process:
            process.communicate()
            if process.returncode != 0:
                raise SystemError(
                    "Could not find 'cmake' on your PATH, cannot build hyperscan"
                )

        with Popen(["which", "ninja"], stdin=PIPE, stdout=PIPE, stderr=PIPE) as process:
            process.communicate()
            if process.returncode != 0:
                raise SystemError(
                    "Could not find 'ninja' on your PATH, cannot build hyperscan"
                )

        with Popen(
            "git clone -b {tag} --single-branch {rep}".format(
                tag=LIBHS_VERSION, rep=LIBHS_SOURCE
            ).split(),
            stdin=PIPE,
            stdout=PIPE,
            stderr=STDOUT,
        ) as process:
            for line in process.stdout:
                print(line.decode("utf-8", errors="ignore")[:-1])
            _, err = process.communicate()
            if process.returncode != 0:
                raise SystemError(
                    "Error cloning hyperscan: {err}".format(
                        err=err.decode("utf-8", errors="ignore")
                    )
                )

        if not os.path.exists("hyperscan"):
            raise SystemError("Error, hyperscan directory does not exist")

        os.mkdir("hyperscan_build")
        os.chdir("hyperscan_build")
        env["CXX"] = cc

        with Popen(
            [
                "cmake",
                "-DCMAKE_BUILD_TYPE=Release",
                "-G",
                "Ninja",
                "-DCMAKE_C_FLAGS_RELEASE='{cflags}'".format(cflags=c_flags),
                "-DCMAKE_CXX_FLAGS_RELEASE='{cppflags}'".format(cppflags=cpp_flags),
                "../hyperscan/",
            ],
            stdin=PIPE,
            stdout=PIPE,
            stderr=STDOUT,
            env=env,
        ) as process:
            for line in process.stdout:
                print(line.decode("utf-8", errors="ignore")[:-1])
            _, err = process.communicate()
            if process.returncode != 0:
                raise SystemError(
                    "Error running cmake: {err}".format(
                        err=err.decode("utf-8", errors="ignore")
                    )
                )

        with Popen(
            ["ninja"], stdin=PIPE, stdout=PIPE, stderr=STDOUT, env=env
        ) as process:
            for line in process.stdout:
                print(line.decode("utf-8", errors="ignore")[:-1])
            _, err = process.communicate()
            if process.returncode != 0:
                raise SystemError(
                    "Error running ninja: {err}".format(
                        err=err.decode("utf-8", errors="ignore")
                    )
                )

        libs = [os.path.join("lib", lib) for lib in LIBHS]
        for lib in libs:
            if not os.path.exists(lib):
                raise SystemError("Error, could not find {lib}".format(lib=lib))
            os.rename(lib, os.path.join("..", lib.split("/")[-1]))

        os.chdir("../")
        try:
            shutil.rmtree("hyperscan")
            shutil.rmtree("hyperscan_build")
        except OSError as e:
            raise SystemError(
                "Error: {name} - {err}.".format(name=e.filename, err=e.strerror)
            )


def main():
    if "--debug" in sys.argv:
        debug = True
        sys.argv.remove("--debug")
    else:
        debug = False

    extra_compile_args = sysconfig.get_config_var("CFLAGS").split()

    cc, o_flag = get_compiler()
    os.environ["CC"] = cc

    if debug:
        o_flag = "0"

    extra_compile_args += [
        "-std=c++17",
        "-O{o}".format(o=o_flag),
        "-lstdc++",
        "-lm",
        "-pthreads",
        "-fPIC",
    ]

    if debug:
        extra_compile_args += ["-g"]

    extra_link_args = ["-lstdc++", "-lm", "-pthread"]

    setup(
        name="HiLine",
        version="0.2.4",
        packages=find_packages(),
        include_package_data=True,
        python_requires=">=3.8.2",
        install_requires=[
            "Click>=7.0",
            "pandas>=1.0.1",
            "numpy>=1.18.1",
            "seaborn>=0.10.0",
            "matplotlib>=3.2.0",
            "biopython>=1.76",
        ],
        description="A HiC alignment and classification pipeline",
        author="Ed Harry",
        author_email="edward.harry@sanger.ac.uk",
        ext_modules=[
            Extension(
                "_HiLine",
                sources=["HiLine.cpp"],
                include_dirs=[os.path.join(os.getcwd(), "include"), np.get_include()],
                library_dirs=[os.getcwd()],
                extra_compile_args=extra_compile_args,
                extra_objects=LIBHS,
                extra_link_args=extra_link_args,
                language="c++17",
            ),
            Extension(
                "_Aligner",
                sources=["ReadTrimmer.cpp"],
                include_dirs=[os.path.join(os.getcwd(), "include")],
                library_dirs=[os.getcwd()],
                extra_compile_args=extra_compile_args,
                extra_link_args=extra_link_args,
                language="c++17",
            ),
        ],
        entry_points={
            "console_scripts": [
                "HiLine=HiLine.main:cli",
                "_HiLine_Aligner=HiLine.aligner:cli",
            ]
        },
        cmdclass={"install": HiLineInstall},
    )


if __name__ == "__main__":
    main()
