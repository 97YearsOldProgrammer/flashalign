#pragma once

#include <string>
#include <string_view>

namespace flashalign {

/**
 * Reverse complement, as mappy's revcomp(): IUPAC codes are complemented (S, W and N map to
 * themselves), case is kept, U becomes A, and any other character is kept as is.
 */
std::string revcomp(std::string_view seq);

}  // namespace flashalign
