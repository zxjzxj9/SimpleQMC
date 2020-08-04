cmake_minimum_required(VERSION 3.10)
project("SimpleQMC")

# See https://stackoverflow.com/questions/6594796/how-do-i-make-cmake-output-into-a-bin-dir
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
#add_library(stub "")
add_subdirectory(src)
# find_package(fmt CONFIG REQUIRED)
# add_definitions(-std=c++11)
# target_link_libraries(main PRIVATE fmt::fmt fmt::fmt-header-only)