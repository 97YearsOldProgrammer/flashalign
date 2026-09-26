# Writes the two version headers (cmake -P). Adds a PEP 440 local suffix to the base version:
#
#   0.1.0                     the v<base> tag is checked out, git is missing or the suffix is OFF
#   0.1.0+g27b015c7           any other commit
#   0.1.0+g27b015c7.dirty     ... with uncommitted changes to tracked files
#
# Outputs are only replaced when they change, so an unchanged version rebuilds nothing.
#
# Required -D arguments:
#   FA_VERSION_BASE     the version parsed out of pyproject.toml
#   FA_SOURCE_DIR       the repository/source root (may be a git worktree)
#   FA_VERSION_H_IN     csrc/version.h.in
#   FA_VERSION_H_OUT    <build>/generated/fa_version.h
#   FA_VERSION_HPP_IN   include/flashalign/version.hpp.in
#   FA_VERSION_HPP_OUT  <build>/generated/flashalign/version.hpp
#   FA_VERSION_SUFFIX   ON/OFF

cmake_minimum_required(VERSION 3.18)

set(FA_VERSION_FULL "${FA_VERSION_BASE}")

if(FA_VERSION_SUFFIX)
    find_program(FA_GIT_EXECUTABLE NAMES git)
    # In a git worktree .git is a file.
    if(FA_GIT_EXECUTABLE AND EXISTS "${FA_SOURCE_DIR}/.git")
        execute_process(
            COMMAND "${FA_GIT_EXECUTABLE}" -C "${FA_SOURCE_DIR}"
                    rev-parse --short HEAD
            RESULT_VARIABLE FA_GIT_RC
            OUTPUT_VARIABLE FA_GIT_SHA
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(FA_GIT_RC EQUAL 0 AND FA_GIT_SHA)
            # The release tag gets the bare version, as a tarball build does.
            execute_process(
                COMMAND "${FA_GIT_EXECUTABLE}" -C "${FA_SOURCE_DIR}"
                        describe --tags --exact-match --match "v*" HEAD
                RESULT_VARIABLE FA_TAG_RC
                OUTPUT_VARIABLE FA_TAG
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE)
            if(NOT (FA_TAG_RC EQUAL 0 AND FA_TAG STREQUAL "v${FA_VERSION_BASE}"))
                set(FA_VERSION_FULL "${FA_VERSION_BASE}+g${FA_GIT_SHA}")
                # Untracked files do not count.
                execute_process(
                    COMMAND "${FA_GIT_EXECUTABLE}" -C "${FA_SOURCE_DIR}"
                            diff --quiet HEAD
                    RESULT_VARIABLE FA_DIRTY_RC
                    OUTPUT_QUIET ERROR_QUIET)
                if(NOT FA_DIRTY_RC EQUAL 0)
                    set(FA_VERSION_FULL "${FA_VERSION_FULL}.dirty")
                endif()
            endif()
        endif()
    endif()
endif()

function(fa_write_if_different template destination)
    set(staged "${destination}.tmp")
    configure_file("${template}" "${staged}" @ONLY)
    execute_process(COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                    "${staged}" "${destination}")
    file(REMOVE "${staged}")
endfunction()

fa_write_if_different("${FA_VERSION_H_IN}" "${FA_VERSION_H_OUT}")
fa_write_if_different("${FA_VERSION_HPP_IN}" "${FA_VERSION_HPP_OUT}")
