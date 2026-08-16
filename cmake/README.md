## Library & Platform Support

DRSolve is distributed as both a **shared library** (`libdrsolve.so` / `libdrsolve.dylib` / `libdrsolve-1.dll`) and a **static library** (`libdrsolve.a` / `libdrsolve-1.a`), alongside the `drsolve` command-line executable.

| Platform | Shared lib | Static lib | CLI |
|---|---|---|---|
| Linux (x86-64, ARM64) | `libdrsolve.so` | `libdrsolve.a` | `drsolve` |
| macOS (x86-64, Apple Silicon) | `libdrsolve.dylib` | `libdrsolve.a` | `drsolve` |
| Windows (x86-64) | `libdrsolve-1.dll` | `libdrsolve-1.a` | `drsolve.exe` + `drsolve_win_gui.exe` |


## Build (quick start)

DRSolve uses **CMake** (≥ 3.16) as its primary build system.

### Linux / macOS

```bash
git clone https://github.com/drsolve/drsolve.git && cd drsolve
cmake -B build
cmake --build build -j$(nproc)
ctest --test-dir build          # optional: run tests
sudo cmake --install build      # optional: install to /usr/local
```

FLINT must be installed (e.g. via your package manager) before configuring.  
If FLINT is in a non-standard location, pass `-DFLINT_ROOT=/path/to/flint`.

```bash
# Example: FLINT installed under $HOME/.local
cmake -B build -DFLINT_ROOT=$HOME/.local
cmake --build build -j$(nproc)
```

### Windows

**Option A — MSYS2/UCRT64**

```bash
pacman -S mingw-w64-ucrt-x86_64-cmake \
          mingw-w64-ucrt-x86_64-gcc \
          mingw-w64-ucrt-x86_64-flint
cmake -B build -G "MinGW Makefiles"
cmake --build build -j$(nproc)
```

**Option B — cross-compile from Linux/macOS (MinGW-w64)**  
Bundled `mingw/` and `pml_det/` dependencies are used automatically:

```bash
cmake -B build-win \
      -DCMAKE_TOOLCHAIN_FILE="$(pwd)/cmake/toolchain-mingw64.cmake"
cmake --build build-win -j$(nproc)
```

---

## Common CMake options

| Option | Default | Description |
|--------|---------|-------------|
| `CMAKE_BUILD_TYPE` | `Release` | `Release` / `Debug` / `RelWithDebInfo` |
| `CMAKE_INSTALL_PREFIX` | `/usr/local` | Install prefix |
| `FLINT_ROOT` | *(auto)* | Root of FLINT installation (contains `include/` and `lib/`) |
| `PML_ROOT` | *(auto)* | Root of PML installation |
| `DRSOLVE_ENABLE_PML` | `ON` | Use PML if found |
| `DRSOLVE_ENABLE_OPENMP` | `ON` | Use OpenMP if available |
| `DRSOLVE_ENABLE_LTO` | `ON` | Link-time optimisation |
| `DRSOLVE_ENABLE_ASAN` | `OFF` | AddressSanitizer |
| `DRSOLVE_STATIC` | `OFF` | Link FLINT/PML statically into `drsolve` |
| `DRSOLVE_BUILD_SHARED_LIB` | `ON` | Build `libdrsolve.so` / `.dylib` |
| `DRSOLVE_BUILD_STATIC_LIB` | `ON` | Build `libdrsolve.a` |
| `DRSOLVE_BUILD_GUI` | `ON` | Build Windows GUI targets |
| `DRSOLVE_BUILD_ATTACK` | `ON` | Build `../Attack/*.c` programs |
| `DRSOLVE_USE_BUNDLED_DEPS` | `ON` | Use bundled third-party dependencies (cross-compile only) |

Examples:

```bash
# Non-standard FLINT location
cmake -B build -DFLINT_ROOT=$HOME/.local

# Fully static binary (no .so deps at runtime)
cmake -B build -DDRSOLVE_STATIC=ON

# Debug build with AddressSanitizer
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DDRSOLVE_ENABLE_ASAN=ON

# Shared library only (skip static library)
cmake -B build -DDRSOLVE_BUILD_STATIC_LIB=OFF

# Install to custom prefix
cmake -B build -DCMAKE_INSTALL_PREFIX=/opt/drsolve
cmake --build build && cmake --install build
```

---

## Build outputs

### Linux / macOS

```
build/
  drsolve                    ← CLI executable
  lib/
    libdrsolve.so.1.0.0        ← shared library (Linux)
    libdrsolve.so → .so.1    ← symlinks
    libdrsolve.dylib           ← shared library (macOS)
    libdrsolve.a               ← static library
```

### Windows

```
build-win/
  drsolve.exe                  ← launcher (double-click or from Explorer)
  drsolve_win_gui.exe          ← GUI frontend
  bin/drsolve_cli_real.exe     ← actual CLI
  dll/libdrsolve-1.dll         ← Dixon shared library
  dll/*.dll                  ← runtime DLL dependencies
  lib/libdrsolve-1.a           ← static archive
  lib/libdrsolve-1.dll.a       ← import library
```

---

## Running tests

```bash
ctest --test-dir build --output-on-failure
```

---

## Uninstall

CMake does not provide a built-in uninstall target.
Use the generated install manifest:

```bash
xargs rm -f < build/install_manifest.txt
```
