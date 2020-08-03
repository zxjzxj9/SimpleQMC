add_definitions(-std=c++11)
find_package(fmt CONFIG REQUIRED)

set(SIMPLE_QMC_SRC simple_qmc.cc)
add_executable(simple_qmc.x ${SIMPLE_QMC_SRC})
target_link_libraries(simple_qmc.x PRIVATE fmt::fmt fmt::fmt-header-only)