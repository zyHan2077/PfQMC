cmake_minimum_required(VERSION 3.13)

PROJECT(MQMC)

include_directories(
${PROJECT_SOURCE_DIR}/inc
)
INCLUDE_DIRECTORIES ( ${EIGEN3_INCLUDE_DIR} )

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(SRC_LIST src/skewMatUtils.cpp src/pfqmc.cpp)

set(MKL_INTERFACE lp64)

set(MKL_ROOT $ENV{MKLROOT})
enable_testing()
find_package(MKL CONFIG REQUIRED)
find_package(GTest CONFIG REQUIRED)
#message(STATUS "${MKL_IMPORTED_TARGETS}") #Provides available list of targets based on input

# add_library( libyaml STATIC IMPORTED )

# set_target_properties( libyaml PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/lib/libyaml-cpp.a )

set(CMAKE_CXX_COMPILER "mpiicpx")
set(CMAKE_C_LINK)
set(FP_MODEL_FLAGS "-fp-model=precise -fimf-arch-consistency=true -no-fma")

set(WARNING_FLAGS "-O3")

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${FP_MODEL_FLAGS} ${WARNING_FLAGS}")

add_executable(main ./main.cpp ${SRC_LIST})
target_compile_options(main PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(main PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(main PUBLIC $<LINK_ONLY:MKL::MKL>)
target_link_libraries(main PUBLIC ${PROJECT_SOURCE_DIR}/inc/pfapack/c_interface/libcpfapack.a ${PROJECT_SOURCE_DIR}/inc/pfapack/fortran/libpfapack.a)

set(TEST_LIST ./test/main_test.cpp ./test/squareLatticeTest.cpp)
add_executable(main_test ${TEST_LIST}  ${SRC_LIST})
target_compile_options(main_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(main_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(main_test PUBLIC $<LINK_ONLY:MKL::MKL> GTest::gtest_main)
target_link_libraries(main_test PUBLIC ${PROJECT_SOURCE_DIR}/inc/pfapack/c_interface/libcpfapack.a ${PROJECT_SOURCE_DIR}/inc/pfapack/fortran/libpfapack.a)

include(GoogleTest)
gtest_discover_tests(main_test)

# add_test(NAME mytest COMMAND myapp)
# if(MKL_ENV)
#   set_tests_properties(mytest PROPERTIES ENVIRONMENT "${MKL_ENV}")
# endif()
# add_executable(rapidtest src/rapidjsontest.cpp)
# target_link_libraries(main libyaml)