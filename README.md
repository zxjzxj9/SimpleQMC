# Simple Quantum Monte Carlo Programs
![build status](https://travis-ci.com/zxjzxj9/SimpleQMC.svg?branch=master)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Goals for this repo
To create a series of Quantum Monte Carlo programs. The programs are for demonstration only.

## How to build:
1. Install vcpkg, see the vcpkg [official repo](https://github.com/microsoft/vcpkg).
2. Using cmake and make to build the executables.

```bash
mkdir build
cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=/path/to/your/vckpg/vcpkg.cmake
```