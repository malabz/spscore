cmake_minimum_required(VERSION 3.18)
project(spscore)

set(CMAKE_CXX_STANDARD 17)

# Find zlib
find_package(ZLIB REQUIRED)

# Find OpenMP
find_package(OpenMP)

# Create a library with the core functionality (no main)
add_library(spscore_lib STATIC src/spscore.cpp)
target_include_directories(spscore_lib PUBLIC src)
target_link_libraries(spscore_lib PUBLIC ZLIB::ZLIB)
if(OpenMP_CXX_FOUND)
    target_link_libraries(spscore_lib PUBLIC OpenMP::OpenMP_CXX)
endif()
# Define NO_MAIN for library to exclude main function
target_compile_definitions(spscore_lib PRIVATE NO_MAIN)

# Create the main executable
add_executable(spscore src/spscore.cpp)
target_include_directories(spscore PRIVATE src)

# Fix: Link libstdc++ and libgcc statically to avoid GLIBCXX version mismatch
# often caused by Conda or older system libraries.
if(UNIX)
    target_link_libraries(spscore PRIVATE -static-libstdc++ -static-libgcc)
endif()

target_link_libraries(spscore PRIVATE ZLIB::ZLIB)
if(OpenMP_CXX_FOUND)
    target_link_libraries(spscore PRIVATE OpenMP::OpenMP_CXX)
endif()

