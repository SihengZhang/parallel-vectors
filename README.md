# parallel-vectors
An implementation of 3D parallel vector curves extractor.

## VTK Submodule
This project integrates the **Visualization Toolkit (VTK)** as a Git submodule. 
To ensure a clean development environment and fast daily compile times, we use a "Vendored" workflow: we build VTK once into a local directory and link our application against it.

## 0. Prerequisites
### macOS
Xcode Command Line Tools:
```bash
xcode-select --install
```
Build Tools (via Homebrew):
```bash
brew install cmake ninja
```
### Linux (Ubuntu/Debian)
```bash
sudo apt-get update
sudo apt-get install build-essential cmake ninja-build libgl1-mesa-dev libxt-dev
```
### Windows
1. Visual Studio 2022: Install with the "Desktop development with C++" workload.
2. CMake: Install from [cmake.org](https://cmake.org/download/)
## 1. Clone & Initialize

Clone the repository and initialize the VTK submodule:

```bash
git clone https://github.com/SihengZhang/parallel-vectors.git
cd parallel-vectors

# Initialize the VTK submodule
git submodule update --init --recursive
```
## 2. Build VTK
### macOS
```bash
mkdir -p external/VTK-build && cd external/VTK-build

# Configure (Enable Cocoa for Mac GUI)
cmake ../VTK \
  -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DVTK_BUILD_TESTING=OFF \
  -DVTK_BUILD_EXAMPLES=OFF \
  -DVTK_USE_COCOA=ON

# Build
ninja
```
### Linux (Ubuntu/Debian)
```bash
mkdir -p external/VTK-build && cd external/VTK-build

cmake ../VTK \
  -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DVTK_BUILD_TESTING=OFF \
  -DVTK_BUILD_EXAMPLES=OFF

ninja
```
### Windows
```PowerShell
mkdir external/VTK-build
cd external/VTK-build

# Note: We do NOT use Ninja here to ensure Visual Studio compatibility
cmake ../VTK -DCMAKE_BUILD_TYPE=Release -DVTK_BUILD_TESTING=OFF -DVTK_BUILD_EXAMPLES=OFF

# Build (Release mode is critical on Windows)
cmake --build . --config Release --parallel 8
```
## 3. Build
```bash
mkdir build && cd build

# Tell CMake where to find our local VTK
cmake .. -GNinja -DVTK_DIR=$PWD/../external/vtk-build

ninja
./MyApp
```

