find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
add_executable(EC_Sec ecuacionCalor2D.cpp)
target_link_libraries(EC_Sec ${ARMADILLO_LIBRARIES})
