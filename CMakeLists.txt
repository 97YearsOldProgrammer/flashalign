cmake_minimum_required(VERSION 3.18)

# The release version is read from pyproject.toml, which must hold exactly one.
file(STRINGS "${CMAKE_CURRENT_LIST_DIR}/pyproject.toml" FA_VERSION_LINES
    REGEX "^version = \"[0-9.]+\"$")
list(LENGTH FA_VERSION_LINES FA_VERSION_LINE_COUNT)
if(NOT FA_VERSION_LINE_COUNT EQUAL 1)
    message(FATAL_ERROR
        "pyproject.toml must contain exactly one numeric project version")
endif()
string(REGEX REPLACE "^version = \"([0-9.]+)\"$" "\\1"
    FA_PROJECT_VERSION "${FA_VERSION_LINES}")

project(FlashAlign VERSION "${FA_PROJECT_VERSION}" LANGUAGES C CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Version headers: cmake/GenerateVersion.cmake adds a "+g<sha>[.dirty]" suffix unless the tree
# is the release tag. It runs at configure time and again on every build (fa_version).
option(FLASHALIGN_VERSION_SUFFIX
    "Append a git commit suffix (+g<sha>[.dirty]) to the reported version" ON)
set(FA_GENERATED_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
set(FA_GENERATE_VERSION_COMMAND
    ${CMAKE_COMMAND}
    -DFA_VERSION_BASE=${FA_PROJECT_VERSION}
    -DFA_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DFA_VERSION_H_IN=${CMAKE_CURRENT_SOURCE_DIR}/csrc/version.h.in
    -DFA_VERSION_H_OUT=${FA_GENERATED_INCLUDE_DIR}/fa_version.h
    -DFA_VERSION_HPP_IN=${CMAKE_CURRENT_SOURCE_DIR}/include/flashalign/version.hpp.in
    -DFA_VERSION_HPP_OUT=${FA_GENERATED_INCLUDE_DIR}/flashalign/version.hpp
    -DFA_VERSION_SUFFIX=${FLASHALIGN_VERSION_SUFFIX}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/GenerateVersion.cmake)
file(MAKE_DIRECTORY ${FA_GENERATED_INCLUDE_DIR}/flashalign)
execute_process(COMMAND ${FA_GENERATE_VERSION_COMMAND}
    RESULT_VARIABLE FA_GENERATE_VERSION_RC)
if(NOT FA_GENERATE_VERSION_RC EQUAL 0)
    message(FATAL_ERROR "failed to generate the version headers")
endif()
add_custom_target(fa_version ALL
    COMMAND ${FA_GENERATE_VERSION_COMMAND}
    COMMENT "Refreshing the generated version headers"
    VERBATIM)


include(cmake/Dependencies.cmake)

# ksw2 DP kernels, adapted from minimap2 (ksw2_ll_sse.c unmodified). SSE4.1 on x86; on ARM
# through the bundled sse2neon, with -fsigned-char because ksw2 assumes a signed char.
set(KSW2_DIR ${CMAKE_CURRENT_SOURCE_DIR}/csrc/dp)
set(KSW2_SOURCES
    ${KSW2_DIR}/ksw2_extz2_sse.c
    ${KSW2_DIR}/ksw2_extd2_sse.c
    ${KSW2_DIR}/ksw2_exts2_sse.c
    ${KSW2_DIR}/ksw2_ll_sse.c)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64|arm64|ARM64)$")
    set_source_files_properties(${KSW2_SOURCES} PROPERTIES
        COMPILE_OPTIONS "-DKSW_SSE2_ONLY;-D__SSE2__;-fsigned-char;-I${KSW2_DIR}/sse2neon")
    message(STATUS "FlashAlign: ksw2 built for ARM/NEON via bundled sse2neon shim.")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64|i[3-6]86)$")
    include(CheckCCompilerFlag)
    check_c_compiler_flag("-msse4.1" FLASHALIGN_C_HAS_MSSE41)
    if(NOT FLASHALIGN_C_HAS_MSSE41)
        message(FATAL_ERROR
            "The selected x86 C compiler does not support ksw2's -msse4.1 flag")
    endif()
    set_source_files_properties(${KSW2_SOURCES} PROPERTIES COMPILE_OPTIONS "-msse4.1")
    message(STATUS "FlashAlign: ksw2 built for x86 via SSE4.1.")
else()
    message(FATAL_ERROR
        "Unsupported ksw2 architecture: ${CMAKE_SYSTEM_PROCESSOR}")
endif()

# ksw2 and its arena allocator as one PIC object library with hidden symbols, shared by the
# CLI and the Python module. -O3 whatever the build type.
add_library(fa_ksw2 OBJECT ${KSW2_SOURCES} ${KSW2_DIR}/ksw_arena.c)
target_include_directories(fa_ksw2 PRIVATE csrc)
target_compile_options(fa_ksw2 PRIVATE -O3)
set_target_properties(fa_ksw2 PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    C_VISIBILITY_PRESET hidden)

include(cmake/CompilerOptions.cmake)
include(cmake/Modules.cmake)

# The public C++ API, compiled once and linked into both the shared libflashalign and the
# Python extension, so the wheel holds a single shared object.
add_library(flashalign_public OBJECT csrc/api/public.cpp csrc/api/fastx.cpp)
set_target_properties(flashalign_public PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(flashalign_public PUBLIC cxx_std_17)
target_include_directories(flashalign_public
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/csrc)
target_link_libraries(flashalign_public PUBLIC
    flashalign_api flashalign_io)

add_library(flashalign_library SHARED $<TARGET_OBJECTS:flashalign_public>)
set_target_properties(flashalign_library PROPERTIES
    OUTPUT_NAME flashalign
    EXPORT_NAME flashalign)
target_compile_features(flashalign_library PUBLIC cxx_std_17)
target_include_directories(flashalign_library
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/generated>
      $<INSTALL_INTERFACE:include>)
target_link_libraries(flashalign_library PRIVATE
    flashalign_api flashalign_io)
add_library(FlashAlign::flashalign ALIAS flashalign_library)

# The command-line tool.
add_executable(flashalign
    csrc/cli/main.cpp
    csrc/cli/parse.cpp
    csrc/cli/help.cpp
    csrc/cli/option_registry.cpp
    csrc/cli/align_command.cpp
    csrc/cli/output_writer.cpp
    csrc/cli/index_command.cpp
    csrc/cli/report.cpp
)
target_include_directories(flashalign PRIVATE csrc ${FA_GENERATED_INCLUDE_DIR})
# Both targets compile fa_version.h, which fa_version rewrites at build time.
add_dependencies(flashalign fa_version)
add_dependencies(flashalign_io fa_version)
target_link_libraries(flashalign PRIVATE
    flashalign_api flashalign_io)

# The nanobind Python extension. pyproject.toml turns it on for `pip install .`.
option(FLASHALIGN_BUILD_PYTHON "Build the nanobind Python extension (_flashalign)" OFF)
if(FLASHALIGN_BUILD_PYTHON)
    find_package(
        Python 3.12 REQUIRED
        COMPONENTS Interpreter Development.Module Development.SABIModule)
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -c "import nanobind; print(nanobind.cmake_dir())"
        OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE _nb_cmake_dir)
    list(APPEND CMAKE_PREFIX_PATH "${_nb_cmake_dir}")
    find_package(nanobind CONFIG REQUIRED)
    # Stable ABI: one wheel per platform for Python 3.12+.
    nanobind_add_module(
        _flashalign STABLE_ABI NB_STATIC
        csrc/python/module.cpp
        csrc/python/bind_index.cpp
        csrc/python/bind_config.cpp
        csrc/python/bind_aligner.cpp
        csrc/python/bind_alignment.cpp
        csrc/python/bind_fastx.cpp)
    target_include_directories(
        _flashalign PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_BINARY_DIR}/generated)
    # Linked in statically rather than against libflashalign, which auditwheel cannot bundle.
    target_link_libraries(_flashalign PRIVATE flashalign_public)
    # module.cpp includes the generated flashalign/version.hpp.
    add_dependencies(_flashalign fa_version)
    # Export no zlib or engine symbols, so they cannot clash with another zlib loaded in the
    # same interpreter.
    set_target_properties(_flashalign PROPERTIES
        C_VISIBILITY_PRESET hidden CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
    if(NOT APPLE AND NOT WIN32)
        target_link_options(_flashalign PRIVATE -Wl,--exclude-libs,ALL)
    endif()
    install(
        TARGETS _flashalign
        LIBRARY DESTINATION flashalign COMPONENT Python)
endif()

install(TARGETS flashalign_library
    EXPORT FlashAlignTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT SharedLibrary
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT SharedLibrary
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(
    TARGETS flashalign
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT CLI)
install(FILES flashalign.1
    DESTINATION ${CMAKE_INSTALL_MANDIR}/man1
    COMPONENT CLI)
install(DIRECTORY include/flashalign/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/flashalign
    COMPONENT Development
    FILES_MATCHING PATTERN "*.hpp")
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/generated/flashalign/version.hpp"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/flashalign
    COMPONENT Development)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FlashAlignConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/FlashAlignConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/FlashAlign")
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/FlashAlignConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion)
install(EXPORT FlashAlignTargets
    FILE FlashAlignTargets.cmake
    NAMESPACE FlashAlign::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/FlashAlign"
    COMPONENT Development)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/FlashAlignConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/FlashAlignConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/FlashAlign"
    COMPONENT Development)
