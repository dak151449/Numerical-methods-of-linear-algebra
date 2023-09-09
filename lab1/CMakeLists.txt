cmake_minimum_required(VERSION 3.20)
project(lab1)

set(CMAKE_CXX_STANDARD 11)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")



find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
include_directories(${Python3_INCLUDE_DIRS})
add_executable(result_plot main.cpp "Matrix.h" "Matrix.cpp" "Func_image.h" "Func_image.cpp")
target_link_libraries(result_plot ${Python3_LIBRARIES} -larmadillo )


#target_link_libraries(lab1 -I /usr/include/python3.10 -lpython3.10)
# add_executable(lab1
#                     main.cpp "Matrix.h" "Matrix.cpp")