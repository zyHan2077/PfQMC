cmake_minimum_required(VERSION 3.13)

PROJECT(MQMC)

set(ENV{MKLROOT} "/DATA.c9/opt/intel/oneapi/mkl")

execute_process(COMMAND which icpc)

include_directories(
${PROJECT_SOURCE_DIR}/inc
)
INCLUDE_DIRECTORIES ( ${EIGEN3_INCLUDE_DIR} )

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(SRC_LIST src/skewMatUtils.cpp src/operator.cpp src/pfqmc.cpp)

set(MKL_INTERFACE lp64)

# add_library( libyaml STATIC IMPORTED )

# set_target_properties( libyaml PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/lib/libyaml-cpp.a )

set(CMAKE_CXX_COMPILER "mpiicpc -mkl")
set(CMAKE_C_LINK)
set(FP_MODEL_FLAGS "-fp-model=precise -fimf-arch-consistency=true -no-fma")

# set(WARNING_FLAGS "-Wall")

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${FP_MODEL_FLAGS} ${WARNING_FLAGS}")

add_executable(main ./main.cpp ${SRC_LIST})
#target_compile_options(main PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
#target_include_directories(main PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
#target_link_libraries(main PUBLIC $<LINK_ONLY:MKL::MKL>)

# add_test(NAME mytest COMMAND myapp)
# if(MKL_ENV)
#   set_tests_properties(mytest PROPERTIES ENVIRONMENT "${MKL_ENV}")
# endif()
# add_executable(rapidtest src/rapidjsontest.cpp)
# target_link_libraries(main libyaml)