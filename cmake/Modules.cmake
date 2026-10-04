include_guard(GLOBAL)

add_library(flashalign_core STATIC
    csrc/core/cigar.cpp
    csrc/core/sequence.cpp)
flashalign_configure_internal_target(flashalign_core)

add_library(flashalign_threading INTERFACE)
add_library(flashalign_split INTERFACE)

add_library(flashalign_index STATIC
    csrc/index/builder.cpp
    csrc/index/faix_format.cpp
    csrc/index/faix_storage_posix.cpp)
flashalign_configure_internal_target(flashalign_index)
target_link_libraries(flashalign_index PRIVATE
    flashalign_core flashalign_threading)

add_library(flashalign_chaining STATIC
    csrc/chaining/colinear_chain.cpp
    csrc/chaining/dense_chain.cpp
    csrc/chaining/partition.cpp)
flashalign_configure_internal_target(flashalign_chaining)
target_link_libraries(flashalign_chaining PRIVATE
    flashalign_core flashalign_split)
# dense_chain.cpp's vector code must round exactly as its scalar code, so no FMA contraction.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set_source_files_properties(csrc/chaining/dense_chain.cpp PROPERTIES
        COMPILE_OPTIONS "-ffp-contract=off")
endif()
add_library(flashalign_seeding STATIC
    csrc/seeding/posting_density.cpp)
flashalign_configure_internal_target(flashalign_seeding)
target_link_libraries(flashalign_seeding PRIVATE
    flashalign_core flashalign_index flashalign_chaining
    flashalign_split)

add_library(flashalign_voting STATIC
    csrc/voting/query_tiles.cpp
    csrc/voting/candidate_catalogue.cpp
    csrc/voting/query_partition.cpp)
flashalign_configure_internal_target(flashalign_voting)
target_link_libraries(flashalign_voting PUBLIC
    flashalign_core flashalign_index flashalign_seeding)

add_library(flashalign_dp STATIC
    csrc/dp/splice_kernel.cpp
    $<TARGET_OBJECTS:fa_ksw2>)
flashalign_configure_internal_target(flashalign_dp)
target_link_libraries(flashalign_dp PRIVATE
    flashalign_core)

add_library(flashalign_dna STATIC
    csrc/dna/alternative_hypothesis.cpp
    csrc/dna/chain_mapq.cpp
    csrc/dna/postdp_scoring.cpp
    csrc/dna/mapq_routing.cpp
    csrc/dna/retained_seed_density.cpp
    csrc/dna/inv_local_chain.cpp
    csrc/dna/placement_chaining.cpp
    csrc/dna/placement_family_adapter.cpp
    csrc/dna/record_family.cpp
    csrc/dna/family_projection.cpp
    csrc/dna/family_realization.cpp
    csrc/dna/residue_emission.cpp
    csrc/dna/residue_trigger.cpp
    csrc/dna/dp_runner.cpp
    csrc/dna/realization_controller.cpp
    csrc/dna/ordered_anchor_path.cpp
    csrc/dna/backend.cpp)
flashalign_configure_internal_target(flashalign_dna)
target_link_libraries(flashalign_dna PRIVATE
    flashalign_core flashalign_index flashalign_dp flashalign_chaining
    flashalign_seeding flashalign_voting flashalign_split)
target_compile_definitions(flashalign_dna PRIVATE FLASHALIGN_BUILDING_DNA)
add_library(flashalign_rna STATIC
    csrc/rna/annotation/known_junctions.cpp
    csrc/rna/placement/aggregate.cpp
    csrc/rna/placement/fused_capture.cpp
    csrc/rna/placement/coarse_chain.cpp
    csrc/rna/placement/output.cpp
    csrc/rna/anchoring/skeleton_regions.cpp
    csrc/rna/anchoring/skeleton_harvest.cpp
    csrc/rna/anchoring/selected_locus_anchoring.cpp
    csrc/rna/anchoring/exact_anchor_path.cpp
    csrc/rna/anchoring/select_exact_anchor_path.cpp
    csrc/rna/anchoring/frame_diagonal_support.cpp
    csrc/rna/config.cpp
    csrc/rna/chain_mapq.cpp
    csrc/rna/query_partition.cpp
    csrc/rna/rival_pricing.cpp
    csrc/rna/output.cpp
    csrc/rna/realization/splice_anchor_view.cpp
    csrc/rna/realization/splice_controller.cpp
    csrc/rna/realization/splice_realizer.cpp
    csrc/rna/realization/rank2_arbitration.cpp
    csrc/rna/realization/rival_lifecycle.cpp
    csrc/rna/backend.cpp)
flashalign_configure_internal_target(flashalign_rna)
target_link_libraries(flashalign_rna PRIVATE
    flashalign_core flashalign_index flashalign_dp flashalign_chaining
    flashalign_seeding flashalign_voting fa_zlib)
target_compile_definitions(flashalign_rna PRIVATE FLASHALIGN_BUILDING_RNA)

add_library(flashalign_options STATIC csrc/options/resolve.cpp)
flashalign_configure_internal_target(flashalign_options)
target_link_libraries(flashalign_options PRIVATE
    flashalign_core flashalign_seeding)

add_library(flashalign_engine INTERFACE)
target_link_libraries(flashalign_engine INTERFACE
    flashalign_core flashalign_index flashalign_dna flashalign_rna
    flashalign_threading)

add_library(flashalign_io STATIC
    csrc/io/byte_source.cpp
    csrc/io/deflate_split.cpp
    csrc/io/gz_member_source.cpp
    csrc/io/paf.cpp
    csrc/io/sam.cpp)
flashalign_configure_internal_target(flashalign_io)
# Threads for gz_member_source.cpp's decompression threads.
target_link_libraries(flashalign_io PRIVATE
    flashalign_core fa_zlib Threads::Threads)

add_library(flashalign_api STATIC
    csrc/api/aligner.cpp
    csrc/api/reference_kind.cpp)
flashalign_configure_internal_target(flashalign_api)
target_link_libraries(flashalign_api PRIVATE
    flashalign_engine flashalign_index flashalign_options flashalign_io
    flashalign_threading Threads::Threads fa_zlib)
