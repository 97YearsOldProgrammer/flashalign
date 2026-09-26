include_guard(GLOBAL)

find_package(Threads REQUIRED)

# zlib for gzip, BGZF and BAM; 1.2.9 is the first release with crc32_z.
find_package(ZLIB 1.2.9 REQUIRED)
add_library(fa_zlib INTERFACE)
target_link_libraries(fa_zlib INTERFACE ZLIB::ZLIB)
