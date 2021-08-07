# Twist-n-Sync CPP Module
## Abstract

![workflow](https://github.com/achains/Twist-n-Sync-CPP-Module/actions/workflows/cmake-ubuntu.yml/badge.svg)

[Twist-n-Sync] is an Android application for synchronizing time on several Android smartphones based on their gyro's data.
During Summer School held by Mathematics & Mechanics Faculty of Saint-Petersburg State University, I wrote this C++ library and [Pavel Mokeev] integrated it into the original application for [Skolteh Mobile Robotics].
Earlier in Twist-n-Sync all calculations were perfomed on external server via Python script. From now our native code can be run directly on smartphones.


## External dependencies
In order to successfully port numpy and scipy methods, which were used in origin, we analyzed several numerical C++ libraries and came to:

- [Eigen] - C++ template library for linear algebra
- [Spline] - library for spline interpolation tasks
- [googletest] - Google's C++ test framework (is needed if you want to build and run tests)

These libraries are fetched by CMake, so you don't need to pre-download them. 

## Prerequisites

| Tool | Version |
| ------ | ------ |
| gcc / clang | 5+ / 3.4+ |
| CMake | 3.14+ |

Project was built and tested on the Ubuntu 21.04 with the above versions, but I believe you may try earlier versions of compilers. You may also try other OS, everthing still should work.
CMake version is mandatory in order to FetchContent command may work properly.

## Installation

First, clone this library on your device. And then build it with CMake:

```sh
cmake -B build [-options]
cmake --build build
```
Custom CMake options are supported:
| Flag | Value | Meaning |
| ------ | ------ | ------ |
| -DBUILD_TESTS | ON / OFF | Build test executable, OFF by default |
| -DUSE_ASAN | ON / OFF | Build with [ASan] verification, OFF by default |

By default project is built in Debug mode. Specify -DCMAKE_BUILD_TYPE option if you want to build in Release mode with -O2 optimization on.
## License

MIT License

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [Twist-n-Sync]: https://github.com/MobileRoboticsSkoltech/twist-n-sync 
   [Pavel Mokeev]: https://github.com/pmokeev
   [Eigen]: https://gitlab.com/libeigen/eigen.git
   [Spline]: https://github.com/ttk592/spline.git
   [googletest]: https://github.com/google/googletest
   [ASan]: https://github.com/google/sanitizers/wiki/AddressSanitizer
   [Skolteh Mobile Robotics]: https://sites.skoltech.ru/mobilerobotics/
