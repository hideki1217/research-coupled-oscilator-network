#!/bin/bash

if [ ! -d "external" ]; then
    mkdir external
fi

pushd external > /dev/null

if [ ! -d "boost_1_83_0" ]; then
    wget https://boostorg.jfrog.io/artifactory/main/release/1.83.0/source/boost_1_83_0.tar.gz
    tar -zxvf boost_1_83_0.tar.gz
    rm boost_1_83_0.tar.gz
else 
    echo "SKIP: 'boost' is already installed."
fi

if [ ! -d "thread-pool-4.1.0" ]; then
    wget https://github.com/bshoshany/thread-pool/archive/refs/tags/v4.1.0.tar.gz
    tar -zxvf v4.1.0.tar.gz
    rm v4.1.0.tar.gz
else
    echo "SKIP: 'bshoshany/thread-pool' is already installed."
fi

popd > /dev/null

if [ ! -d "venv" ]; then
    python3 -m venv venv
    source ./venv/bin/activate
    pip install -r requirements.txt
else 
    echo "SKIP: 'venv' is already created."
fi