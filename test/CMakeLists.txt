cmake_minimum_required(VERSION 3.18)
project(spscore_test)

set(CMAKE_CXX_STANDARD 17)

# Find zlib
find_package(ZLIB REQUIRED)

# Find OpenMP
find_package(OpenMP)

# Include parent source directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Create the test executable
add_executable(spscore_test
    test_spscore.cpp
    ../src/spscore.cpp
)

# Fix: Link libstdc++ and libgcc statically to avoid GLIBCXX version mismatch
# often caused by Conda or older system libraries.

target_link_libraries(spscore_test PRIVATE -static-libstdc++ -static-libgcc)


# Link libraries
target_link_libraries(spscore_test PRIVATE ZLIB::ZLIB)
if(OpenMP_CXX_FOUND)
    target_link_libraries(spscore_test PRIVATE OpenMP::OpenMP_CXX)
endif()

# Exclude main() when building test
target_compile_definitions(spscore_test PRIVATE NO_MAIN)

# Enable testing
enable_testing()
add_test(NAME spscore_test COMMAND spscore_test)
