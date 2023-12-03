cmake_minimum_required(VERSION 3.10)

find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(use_kdl euler_demo.cpp)
target_link_libraries(use_kdl orocos-kdl)

add_executable(fanuc_kine fanuc.cpp)
target_link_libraries(fanuc_kine orocos-kdl)

# add_executable(demo demo.cpp)
# target_link_libraries(demo orocos-kdl)
add_executable(ik ik_demo.cpp KineSolve.cpp)