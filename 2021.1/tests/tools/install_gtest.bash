#!/bin/bash

# Copyright (c) 2012, Chris Rogers
# All rights reserved.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 3. Neither the name of STFC nor the names of its contributors may be used to
#    endorse or promote products derived from this software without specific
#    prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

here=`pwd`
vers=1.7.0

check_dependencies() {
    dep_list="wget unzip md5sum"
    echo "Checking that dependencies are installed"
    for dep in ${dep_list}; do
        echo "Checking dependency ${dep}..."
        ${dep} --help > /dev/null
        if [ $? -ne 0 ]; then
            echo "FAIL: require ${dep} to run script"
            exit 1
        fi
        echo "                        ok"
    done
}

download_gtest() {
    source_dir=${here}/gtest
    gtest_dir=${source_dir}/googletest-release-${vers}
    file=release-${vers}.zip
    url=https://github.com/google/googletest/archive/${file}
    md5="ef5e700c8a0f3ee123e2e0209b8b4961  ${file}"

    echo "Get a new copy of the gtest library if required"
    mkdir -p ${source_dir}
    cd ${source_dir}

    echo "${md5}" | md5sum -c
    if [ $? -ne 0 ]; then
        wget ${url}
    fi
    echo "${md5}" | md5sum -c
    if [ $? -ne 0 ]; then
        echo "FAIL: failed to download gtest library"
        exit 1
    fi
    echo "Clean existing gtest build"
    rm -rf ${gtest_dir}
    echo "Unzipping gtest source"
    unzip ${file}
}

build_gtest() {
    build_dir=build
    echo "Build the gtest library"
    rm -r ${build_dir}
    mkdir ${build_dir}
    cd ${build_dir}
    GTEST_ROOT=${gtest_dir} cmake -DCMAKE_INSTALL_PREFIX=${source_dir} ${gtest_dir}

    make
    if [ $? -ne 0 ]; then
        echo "FAIL: failed to build gtest"
        exit 1
    fi
}

install_gtest() {

    install_dir=${source_dir}
    inc_dir=${install_dir}/include/
    lib_dir=${install_dir}/lib/

    echo "Installing the gtest library from ${source_dir} to ${install_dir}"

    # I toyed with the idea of doing a cleanup here, decided against it in the
    # end; so we overwrite but don't clean existing
    mkdir -p ${inc_dir}
    cp -rf ${gtest_dir}/include/* ${install_dir}/include/
    mkdir -p ${install_dir}/lib/
    cp -r ${source_dir}/${build_dir}/*.a ${install_dir}/lib/

    echo "add"
    echo "export GTEST_ROOT=${install_dir}"
    echo "to your .bashrc or"
    echo "GTEST_ROOT=${install_dir}"
    echo "when you call cmake"
}

main() {
    check_dependencies
    download_gtest
    build_gtest
    install_gtest
}

main
