cmake_minimum_required(VERSION 3.12)
project(parallel-vectors)

# 1. Tell CMake where to find the local VTK build
# (We will also set this in VS Code settings, but this is a fallback)
set(VTK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/VTK-build" CACHE PATH "Path to local VTK build")

# 2. Find the VTK package
find_package(VTK REQUIRED)

if (APPLE)
    # Tell CMake to use RPATH (path to libs) in the build tree
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    
    # Add the path to your local VTK dylibs so the app finds them
    set(CMAKE_BUILD_RPATH "${CMAKE_CURRENT_SOURCE_DIR}/external/VTK/build/lib")
    
    # Enable standard macOS RPATH behavior
    set(CMAKE_MACOSX_RPATH ON)
endif()

# 3. Add your executable
add_executable(MyApp src/main.cpp
                     src/interval.h)

# 4. Link VTK to your executable
target_link_libraries(MyApp PRIVATE ${VTK_LIBRARIES})

# 5. Initialize VTK modules (Required for VTK 9+)
vtk_module_autoinit(
  TARGETS MyApp
  MODULES ${VTK_LIBRARIES}
)