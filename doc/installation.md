# Installation

[Documentation index](README.md)

## Dependencies

- A C11 compiler and GNU Make; GCC is commonly used on Linux, Clang on macOS.
- FLINT, preferably version 3.6.0, with GMP and MPFR.
- The required PML subset is bundled; no separate PML installation is needed.
- OpenMP enables supported parallel stages. Builds without OpenMP use serial paths.
- SageMath is needed only for the Sage interface.

See the [FLINT repository](https://github.com/flintlib/flint) for installation instructions. When using distribution packages, install the development headers as well.

## Build with Makefile

Run from the repository root:

```bash
./configure
make -j4
./drsolve --version
./drsolve -h
```

For a custom FLINT installation:

```bash
./configure --flint-include=/path/to/flint/include --flint-lib=/path/to/flint/lib
make -j4
```

Run `./configure --help` and `make help` for further options. Reconfigure and clean old build artifacts when changing compilers or dependency versions.

```bash
make check
# Optional installation into a personal prefix
make install PREFIX="$HOME/.local"
```

`make check` runs regression checks, not a performance benchmark. Runtime depends on the machine and build configuration. Local development uses Makefile; CMake remains maintained for CI and other build environments, including Windows.

## Linking and platforms

The default link mode depends on configuration. `make static` statically links project libraries; `make static-all` attempts to link dependencies statically and requires their static libraries. These names do not guarantee independence from every system runtime on every platform. Check actual dynamic dependencies before deployment.

Linux builds follow the commands above. On macOS, OpenMP requires a runtime compatible with the compiler; serial builds remain usable. Keep compiler and library architectures consistent.

Windows users can use the [release packages](https://github.com/drsolve/drsolve/releases), preserving the executable/DLL layout. Source-build settings are available in the [MinGW toolchain](../cmake/toolchain-mingw64.cmake) and [CI workflow](../.github/workflows/CI.yml). CLI syntax is shared across platforms; in PowerShell, use `./drsolve.exe`.

## Installation check

```bash
./drsolve "x^2-1, y-x" 257 -o /tmp/drsolve-check.dr
```

The system has two solutions over F_257: `(x,y)=(1,1)` and `(256,256)`. On Windows, replace the output path with a local filename.
