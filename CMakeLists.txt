cmake_minimum_required(VERSION 3.10)
project(nbody_cuda C)
list(APPEND CMAKE_C_SOURCE_FILE_EXTENSIONS cu)

set(CMAKE_CXX_STANDARD 11)

#include_directories(SYSTEM "/usr/include")

add_executable(nbody_cuda
        nbody.cu
        rng.cu
        run)
SET_SOURCE_FILES_PROPERTIES(nbody.cu PROPERTIES LANGUAGE C)
SET_SOURCE_FILES_PROPERTIES(rng.cu PROPERTIES LANGUAGE C)

SET_TARGET_PROPERTIES(nbody_cuda PROPERTIES LINKER_LANGUAGE C)

