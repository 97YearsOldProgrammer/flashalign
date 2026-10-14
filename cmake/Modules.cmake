include_guard(GLOBAL)

add_library(flashalign_core STATIC
    csrc/core/core.cpp)
flashalign_configure_internal_target(flashalign_core)

add_library(flashalign_threading INTERFACE)

add_library(flashalign_index STATIC
    csrc/index/build.cpp
    csrc/index/faix.cpp)
flashalign_configure_internal_target(flashalign_index)
target_link_libraries(flashalign_index PRIVATE
    flashalign_core flashalign_threading)

add_library(flashalign_chain STATIC
    csrc/chain/colinear_chain.cpp
    csrc/chain/dense_chain.cpp)
flashalign_configure_internal_target(flashalign_chain)
target_link_libraries(flashalign_chain PRIVATE
    flashalign_core)
# dense_chain.cpp's vector code must round exactly as its scalar code, so no FMA contraction.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set_source_files_properties(csrc/chain/dense_chain.cpp PROPERTIES
        COMPILE_OPTIONS "-ffp-contract=off")
endif()
add_library(flashalign_seed STATIC
    csrc/seed/seed.cpp)
flashalign_configure_internal_target(flashalign_seed)
target_link_libraries(flashalign_seed PRIVATE
    flashalign_core flashalign_index flashalign_chain)

add_library(flashalign_vote STATIC
    csrc/vote/vote.cpp
    csrc/vote/partition.cpp)
flashalign_configure_internal_target(flashalign_vote)
target_link_libraries(flashalign_vote PUBLIC
    flashalign_core flashalign_index flashalign_seed)

add_library(flashalign_dp STATIC
    csrc/dp/dp.cpp
    $<TARGET_OBJECTS:fa_ksw2>)
flashalign_configure_internal_target(flashalign_dp)
target_link_libraries(flashalign_dp PRIVATE
    flashalign_core)

add_library(flashalign_dna STATIC
    csrc/dna/map.cpp
    csrc/dna/seeds.cpp
    csrc/dna/catalogue.cpp
    csrc/dna/placement.cpp
    csrc/dna/ownership.cpp
    csrc/dna/ava.cpp
    csrc/dna/projection.cpp
    csrc/dna/anchor_path.cpp
    csrc/dna/realization.cpp
    csrc/dna/align.cpp
    csrc/dna/mapq.cpp)
flashalign_configure_internal_target(flashalign_dna)
target_link_libraries(flashalign_dna PRIVATE
    flashalign_core flashalign_index flashalign_dp flashalign_chain
    flashalign_seed flashalign_vote)
target_compile_definitions(flashalign_dna PRIVATE FLASHALIGN_BUILDING_DNA)
add_library(flashalign_rna STATIC
    csrc/rna/map.cpp
    csrc/rna/vote.cpp
    csrc/rna/harvest.cpp
    csrc/rna/anchor_path.cpp
    csrc/rna/realization.cpp
    csrc/rna/align.cpp
    csrc/rna/hit.cpp
    csrc/rna/leftover.cpp)
flashalign_configure_internal_target(flashalign_rna)
target_link_libraries(flashalign_rna PRIVATE
    flashalign_core flashalign_index flashalign_dp flashalign_chain
    flashalign_seed flashalign_vote fa_zlib)
target_compile_definitions(flashalign_rna PRIVATE FLASHALIGN_BUILDING_RNA)

add_library(flashalign_options STATIC csrc/options/options.cpp)
flashalign_configure_internal_target(flashalign_options)
target_link_libraries(flashalign_options PRIVATE
    flashalign_core flashalign_seed)

add_library(flashalign_engine STATIC csrc/engine/engine.cpp)
flashalign_configure_internal_target(flashalign_engine)
target_link_libraries(flashalign_engine PUBLIC
    flashalign_core flashalign_index flashalign_dna flashalign_rna
    flashalign_threading)

add_library(flashalign_io STATIC
    csrc/io/reads.cpp
    csrc/io/gzip.cpp
    csrc/io/deflate_split.cpp
    csrc/io/bgzf.cpp
    csrc/io/format.cpp)
flashalign_configure_internal_target(flashalign_io)
# Threads for the decompression teams (gzip.cpp, bgzf.cpp).
target_link_libraries(flashalign_io PRIVATE
    flashalign_core fa_zlib Threads::Threads)

add_library(flashalign_api STATIC
    csrc/api/aligner.cpp)
flashalign_configure_internal_target(flashalign_api)
target_link_libraries(flashalign_api PRIVATE
    flashalign_engine flashalign_index flashalign_options flashalign_io
    flashalign_threading Threads::Threads fa_zlib)
