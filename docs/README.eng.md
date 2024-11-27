<p align="center">
  <img src="./logo.svg" width="400" alt="non-parametric-dsp Logo">
</p>

[ [Русский](../README.md) ]

non-parametric-dsp - A C++ library for processing non-stationary signals.

# Table of Contents
- [Features](#features)
- [Getting Started](#getting-started)
- [Credits](#credits)
- [License](#license)
- [Contacts](#contacts)

# Features
non-parametric-dsp is a C++ library for digital signal processing of non-stationary signals, providing the following features:
- Signal approximation using various methods.
- Numerical differentiation and integration (including fractional orders).
- A wide range of filtering algorithms.
- Instantaneous frequency and amplitude computation using Hilbert and Tikhonov methods.
- Signal phase rotation.
- Flexible empirical mode decomposition.
- Time-frequency analysis.
- Signal tokenization for neural network training.

non-parametric-dsp is a research project that offers many experimental and unique approaches.

For more details on its functionality, please refer to the examples in the `examples` directory.

> [!IMPORTANT]
> To integrate into your projects, please use the "for-linkage" branch; it provides a header-only library without built-in visualization capabilities.

## Getting Started
To begin, install the following dependencies:

#### <img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/git/git-original.svg" width="24" height="24" alt="Git Logo"/> Git

#### <img src="https://xmake.io/assets/img/logo.svg" width="24" height="24" alt="XMake Logo"/> XMake

#### <img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/cplusplus/cplusplus-original.svg" width="24" height="24" alt="C++ Logo"/> C++ Compilers:
- Clang
- GCC
- MSVC

#### <img src="https://images.sftcdn.net/images/t_app-icon-m/p/4f6f9692-96da-11e6-9846-00163ed833e7/2948317542/gnuplot-gnuplot-logo.png" width="24" height="24" alt="Gnuplot Logo"/> Gnuplot

Then, run the following commands:
```bash
git clone https://github.com/Dimandragon/non-parametric-dsp.git
cd non-parametric-dsp
git submodule update --init --recursive
./build_matplot.sh
xmake
```
Finally, try running the examples:
```bash
xmake r [example_name]
```
You can find the example names in the `xmake.lua` file; all examples are "binary" targets in the xmake build system.

# Credits
non-parametric-dsp uses the following projects:
- [pocketfft](https://github.com/mreineck/pocketfft) - A fast and lightweight C++ FFT implementation for different data sizes.
- [icecream](https://github.com/renatoGarcia/icecream-cpp) - A C++ library for simple and convenient debug output formatting.
- [matplotplusplus](https://github.com/alandefreitas/matplotplusplus) - A C++ grapchic library for data visualization.
- [gnuplot](http://www.gnuplot.info/) - A command-line tool for data visualization.
- [boost](https://www.boost.org/) - A large collection of useful C++ libraries.
- [alglib](https://www.alglib.net/) - A cross-platform numerical analysis and data processing library.
- [xmake](https://github.com/xmake-io/xmake) - A user-friendly declarative cross-platform build system written in Lua.

# License
[MIT license](LICENSE)

---

# Contacts
[![Email](https://img.shields.io/badge/Email-D14836?style=flat&logo=gmail&logoColor=white)](mailto:dkuznetsov071105@gmail.com)
```
dkuznetsov071105@gmail.com
```

[![Telegram](https://img.shields.io/badge/Telegram-2CA5E0?style=flat&logo=telegram&logoColor=white)](https://t.me/diman_botan)
```
@diman_botan
```