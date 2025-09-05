cmake_minimum_required(VERSION 3.20)
project(VoxelRPG CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Link MSVC runtime statically to avoid requiring VC++ Redistributable
if(MSVC)
    # Requires CMake 3.15+
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

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
    src/obj_loader.cpp
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

# Copy assets (models/textures) next to the executable
add_custom_command(TARGET voxel POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory $<TARGET_FILE_DIR:voxel>/assets
    COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_SOURCE_DIR}/assets
            $<TARGET_FILE_DIR:voxel>/assets)

# Install rules for packaging/distribution
install(TARGETS voxel RUNTIME DESTINATION .)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/shaders DESTINATION .)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/assets DESTINATION .)
install(FILES ${CMAKE_SOURCE_DIR}/README.md ${CMAKE_SOURCE_DIR}/LICENSE DESTINATION .)

# CPack ZIP package configuration
set(CPACK_PACKAGE_NAME "VoxelRPG")
set(CPACK_PACKAGE_VENDOR "SnowyStorm00")
set(CPACK_PACKAGE_VERSION "0.2.0")
set(CPACK_GENERATOR "ZIP")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-win64")
include(CPack)
