<p align="center">
  <img src="./logo.svg" width="400" alt="non-parametric-dsp Logo">
</p>

[ [Русский](../README.md) ]

non-parametric-dsp - A C++ library for digital processing of non-stationary signals

# Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Getting Started](#getting-started)
- [License](#license)
- [Contact](#contact)
- [Credits](#credits)

# Features
non-parametric-dsp is a C++ library for digital processing of non-stationary signals, offering the following features:
- signal approximation using various methods
- numerical differentiation and integration (including fractional degrees)
- a wide range of filtering algorithms
- computation of instantaneous frequencies and amplitudes using Hilbert and Tikhonov methods
- signal phase shifting
- highly customizable empirical mode decomposition
- time-frequency analysis
- signal tokenization for training neural networks

non-parametric-dsp is a research project offering numerous experimental and unique approaches.

For a detailed overview of non-parametric-dsp's capabilities, please refer to the examples directory.

> [!IMPORTANT]
> When using this library in your projects, please rely on the `for-linkage` vector. It provides a header-only library without built-in visualization features.

## Getting Started
To get started, you need to install the following dependencies:
- git
- xmake
- C++ compiler (clang, gcc, or msvc)
- gnuplot

Then, execute the following commands:
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

You can find example names in the `xmake.lua` file; all examples are `binary` targets in the xmake build system.

# License
[MIT license](LICENSE)

---

# Contact
[![Email](https://img.shields.io/badge/Email-D14836?style=flat&logo=gmail&logoColor=white)](mailto:dkuznetsov071105@gmail.com)
```
dkuznetsov071105@gmail.com
```

[![Telegram](https://img.shields.io/badge/Telegram-2CA5E0?style=flat&logo=telegram&logoColor=white)](https://t.me/diman_botan)
```
@diman_botan
```

# Credits
non-parametric-dsp leverages the following projects:
- [pocketfft](https://github.com/mreineck/pocketfft) - a fast and lightweight C++ FFT implementation for arbitrary input sizes
- [icecream](https://github.com/renatoGarcia/icecream-cpp) - a C++ library for simple and convenient output formatting
- [matplotplusplus](https://github.com/alandefreitas/matplotplusplus) - a C++ data visualization library
- [gnuplot](http://www.gnuplot.info/) - a console-based data visualization
