cmake_minimum_required(VERSION 3.22)
project(research-network-exploration)

# Compiler settings
set(CMAKE_CXX_STANDARD 17)


# External library
add_library(boost INTERFACE)
target_include_directories(boost INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/boost_1_83_0)
target_compile_options(boost INTERFACE "-Wno-tautological-constant-compare" "-Wno-deprecated-declarations")

add_library(thread_pool INTERFACE)
target_include_directories(thread_pool INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/thread-pool-4.1.0/include)

# Owned library
add_library(mine INTERFACE)
target_include_directories(mine INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_definitions(mine INTERFACE MINE_UNIT_N=${MINE_UNIT_N})

# Sub projects
add_subdirectory(network-exploration)