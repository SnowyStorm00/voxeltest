cmake_minimum_required(VERSION 3.20)
project(VoxelRPG CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable debug symbols by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

add_executable(voxel WIN32
    src/main.cpp
    src/d3d_app.cpp
    src/camera.cpp
    src/voxel_world.cpp
    src/voxel_mesher.cpp
)

target_include_directories(voxel PRIVATE include)

# Link Windows and D3D11
if(WIN32)
    target_link_libraries(voxel PRIVATE d3d11 dxgi d3dcompiler)
endif()

# Copy shaders to build dir (root) and next to the executable (per-config)
add_custom_target(copy_shaders ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_SOURCE_DIR}/shaders
            ${CMAKE_CURRENT_BINARY_DIR}/shaders)
add_dependencies(voxel copy_shaders)

add_custom_command(TARGET voxel POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_SOURCE_DIR}/shaders
            $<TARGET_FILE_DIR:voxel>/shaders)
