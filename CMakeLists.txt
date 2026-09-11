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

set(TEST_LIST ./test/main_test.cpp ./test/squareLatticeTest.cpp ./test/performanceTest.cpp)
add_executable(main_test ${TEST_LIST}  ${SRC_LIST})
target_compile_options(main_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(main_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(main_test PUBLIC $<LINK_ONLY:MKL::MKL> GTest::gtest_main)
target_link_libraries(main_test PUBLIC ${PROJECT_SOURCE_DIR}/inc/pfapack/c_interface/libcpfapack.a ${PROJECT_SOURCE_DIR}/inc/pfapack/fortran/libpfapack.a)

include(GoogleTest)
gtest_discover_tests(main_test)

# Independent four-mode Fock-space / Pfaffian validation.
add_executable(sign_free_check "sign_free_verify/four_mode.cpp" src/skewMatUtils.cpp)
target_compile_features(sign_free_check PRIVATE cxx_std_17)
target_link_libraries(sign_free_check PRIVATE MKL::MKL
    ${PROJECT_SOURCE_DIR}/inc/pfapack/c_interface/libcpfapack.a
    ${PROJECT_SOURCE_DIR}/inc/pfapack/fortran/libpfapack.a)
add_test(NAME SignFree.FourMode COMMAND sign_free_check --samples 24)

# add_test(NAME mytest COMMAND myapp)
# if(MKL_ENV)
#   set_tests_properties(mytest PROPERTIES ENVIRONMENT "${MKL_ENV}")
# endif()
# add_executable(rapidtest src/rapidjsontest.cpp)
# target_link_libraries(main libyaml)

# Reproducible simulation benchmark; output is CSV, validation is outside timing.
add_executable(pfqmc_bench bench/pfqmc_bench.cpp ${SRC_LIST})
target_link_libraries(pfqmc_bench PRIVATE MKL::MKL
    ${PROJECT_SOURCE_DIR}/inc/pfapack/c_interface/libcpfapack.a
    ${PROJECT_SOURCE_DIR}/inc/pfapack/fortran/libpfapack.a)

add_library(pfqmc_catalog src/catalog_core.cpp src/catalog_models.cpp src/catalog_worldline.cpp src/catalog_search.cpp src/catalog_projected.cpp ${SRC_LIST})
target_compile_features(pfqmc_catalog PUBLIC cxx_std_17)
target_link_libraries(pfqmc_catalog PUBLIC MKL::MKL
    ${PROJECT_SOURCE_DIR}/inc/pfapack/c_interface/libcpfapack.a
    ${PROJECT_SOURCE_DIR}/inc/pfapack/fortran/libpfapack.a)
add_executable(sign_free_models sign_free_catalog/main.cpp sign_free_catalog/checks.cpp sign_free_catalog/specials.cpp)
target_link_libraries(sign_free_models PRIVATE pfqmc_catalog)
add_test(NAME SignFree.CatalogControls COMMAND sign_free_models --model checks --samples 2 --sweeps 2)
set(CATALOG_MODELS
    open_hopping odd_pair_ring odd_pair_ring_three four_mode
    adjoint_d3 adjoint_d5 adjoint_d7 spin_cell_j1 spin_cell_j2
    spin_tensor_3x3 spin_tensor_5x5 spin_lattice_j1 spin_lattice_j2
    spin_lattice_mixed five_orbital pairing_lattice pairing_lattice_2d
    sl_exterior_4_2 sl_exterior_6_3 sl_exterior_8_4
    sp_exterior_2_2 sp_exterior_4_3 sp_exterior_5_5
    sl2_u2 sl2_u3 sl2_mixed local_unitary_hopping)
foreach(model IN LISTS CATALOG_MODELS)
    add_test(NAME SignFree.Catalog.${model}
        COMMAND sign_free_models --model ${model} --samples 2 --sweeps 2)
endforeach()

add_executable(sign_free_search sign_free_catalog/search.cpp sign_free_catalog/checks.cpp)
target_link_libraries(sign_free_search PRIVATE pfqmc_catalog)
add_test(NAME SignFree.SearchControls COMMAND sign_free_search --model controls --samples 2 --sweeps 2)
foreach(model lorentz4 lorentz5 lorentz8 lorentz4_lattice conjugate_tensor2
    conjugate_tensor3 conjugate_tensor3_lattice symmetric_power3_2
    symmetric_power3_3 symmetric_power4_2 symmetric_power3_2_lattice
    quaternion_herm2 quaternion_herm3 quaternion_herm3_lattice
    quaternion_anti2 quaternion_anti3 quaternion_anti2_lattice
    lorentz_sym5 lorentz_sym6 lorentz_sym5_lattice)
    add_test(NAME SignFree.Search.${model} COMMAND sign_free_search --model ${model} --samples 2 --sweeps 2)
endforeach()

add_executable(sign_free_projected sign_free_catalog/projected.cpp)
target_link_libraries(sign_free_projected PRIVATE pfqmc_catalog)
add_test(NAME SignFree.Projected COMMAND sign_free_projected --samples 2 --seed 731)

target_sources(pfqmc_catalog PRIVATE src/catalog_dynamic_ising.cpp)
add_executable(sign_free_dynamic_ising sign_free_catalog/dynamic_ising.cpp)
target_link_libraries(sign_free_dynamic_ising PRIVATE pfqmc_catalog)
add_test(NAME SignFree.DynamicIsing COMMAND sign_free_dynamic_ising --samples 4 --seed 160507420)

target_sources(pfqmc_catalog PRIVATE src/catalog_poisson.cpp src/catalog_albert.cpp)
add_executable(sign_free_poisson sign_free_catalog/poisson.cpp)
target_link_libraries(sign_free_poisson PRIVATE pfqmc_catalog)
add_test(NAME SignFree.Poisson COMMAND sign_free_poisson --samples 8 --seed 731)
add_executable(sign_free_albert sign_free_catalog/albert.cpp)
target_link_libraries(sign_free_albert PRIVATE pfqmc_catalog)
add_test(NAME SignFree.Albert COMMAND sign_free_albert --model all --samples 2 --sweeps 2)

target_sources(pfqmc_catalog PRIVATE src/catalog_quantum_link.cpp)
add_executable(sign_free_quantum_link sign_free_catalog/quantum_link.cpp)
target_link_libraries(sign_free_quantum_link PRIVATE pfqmc_catalog)
add_test(NAME SignFree.QuantumLink COMMAND sign_free_quantum_link --samples 16 --seed 230508917)

target_sources(pfqmc_catalog PRIVATE src/catalog_mobile_spin.cpp)
add_executable(sign_free_mobile_spin sign_free_catalog/mobile_spin.cpp)
target_link_libraries(sign_free_mobile_spin PRIVATE pfqmc_catalog)
add_test(NAME SignFree.MobileSpin COMMAND sign_free_mobile_spin --samples 8 --seed 160507421)

target_sources(pfqmc_catalog PRIVATE src/catalog_meron.cpp)
add_executable(sign_free_meron sign_free_catalog/meron.cpp)
target_link_libraries(sign_free_meron PRIVATE pfqmc_catalog)
add_test(NAME SignFree.Meron COMMAND sign_free_meron --samples 8 --seed 731)

target_sources(pfqmc_catalog PRIVATE src/catalog_ring_exchange.cpp)
add_executable(sign_free_ring_exchange sign_free_catalog/ring_exchange.cpp)
target_link_libraries(sign_free_ring_exchange PRIVATE pfqmc_catalog)
add_test(NAME SignFree.RingExchange COMMAND sign_free_ring_exchange --samples 16 --seed 409112)

target_sources(pfqmc_catalog PRIVATE src/catalog_color_projector.cpp src/catalog_frustrated_impurity.cpp)
add_executable(sign_free_color_projector sign_free_catalog/color_projector.cpp)
target_link_libraries(sign_free_color_projector PRIVATE pfqmc_catalog)
add_test(NAME SignFree.ColorProjector COMMAND sign_free_color_projector --samples 8 --seed 240602926)
add_executable(sign_free_frustrated_impurity sign_free_catalog/frustrated_impurity.cpp)
target_link_libraries(sign_free_frustrated_impurity PRIVATE pfqmc_catalog)
add_test(NAME SignFree.FrustratedImpurity COMMAND sign_free_frustrated_impurity --samples 8 --seed 160805144)

target_sources(pfqmc_catalog PRIVATE src/catalog_stabilizer.cpp)
add_executable(sign_free_stabilizer sign_free_catalog/stabilizer.cpp)
target_link_libraries(sign_free_stabilizer PRIVATE pfqmc_catalog)
add_test(NAME SignFree.Stabilizer COMMAND sign_free_stabilizer --samples 8 --seed 731)

target_sources(pfqmc_catalog PRIVATE src/catalog_tj_chain.cpp src/catalog_jq.cpp)
add_executable(sign_free_tj_chain sign_free_catalog/tj_chain.cpp)
target_link_libraries(sign_free_tj_chain PRIVATE pfqmc_catalog)
add_test(NAME SignFree.TJChain COMMAND sign_free_tj_chain --samples 8 --seed 13106496)
add_executable(sign_free_jq sign_free_catalog/jq.cpp)
target_link_libraries(sign_free_jq PRIVATE pfqmc_catalog)
add_test(NAME SignFree.JQ COMMAND sign_free_jq --samples 8 --seed 611343)

target_sources(pfqmc_catalog PRIVATE src/catalog_bilayer.cpp src/catalog_zigzag.cpp)
add_executable(sign_free_bilayer sign_free_catalog/bilayer.cpp)
target_link_libraries(sign_free_bilayer PRIVATE pfqmc_catalog)
add_test(NAME SignFree.Bilayer COMMAND sign_free_bilayer --samples 8 --seed 731)
add_executable(sign_free_zigzag sign_free_catalog/zigzag.cpp)
target_link_libraries(sign_free_zigzag PRIVATE pfqmc_catalog)
add_test(NAME SignFree.Zigzag COMMAND sign_free_zigzag --samples 8 --seed 9707019)

target_sources(pfqmc_catalog PRIVATE src/catalog_fredkin.cpp)
add_executable(sign_free_fredkin sign_free_catalog/fredkin.cpp)
target_link_libraries(sign_free_fredkin PRIVATE pfqmc_catalog)
add_test(NAME SignFree.Fredkin COMMAND sign_free_fredkin --samples 8 --seed 160503842)

target_sources(pfqmc_catalog PRIVATE src/catalog_pxp.cpp)
add_executable(sign_free_pxp sign_free_catalog/pxp.cpp)
target_link_libraries(sign_free_pxp PRIVATE pfqmc_catalog)
add_test(NAME SignFree.PXP COMMAND sign_free_pxp --samples 8 --seed 731)
