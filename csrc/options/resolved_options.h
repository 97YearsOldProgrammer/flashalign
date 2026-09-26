// The resolved mapping configuration. The variant holds the mode's options;
// RNA options embed the DNA ones in RnaLongOptions::base.
#pragma once

#include "common_profile.h"
#include "dna_profile.h"
#include "rna_profile.h"

#include <variant>

namespace fa::cpu::options {

struct ResolvedOptions {
  CommonOptions common;
  IndexIdentity index;
  std::variant<::fa::cpu::lr::DnaLongOptions,
               ::fa::cpu::lr::rna::RnaLongOptions>
      mapping;

  bool is_rna() const noexcept {
    return std::holds_alternative<
        ::fa::cpu::lr::rna::RnaLongOptions>(mapping);
  }

  ::fa::cpu::lr::DnaLongOptions& long_read() noexcept {
    if (auto* rna =
            std::get_if<::fa::cpu::lr::rna::RnaLongOptions>(&mapping)) {
      return rna->base;
    }
    return std::get<::fa::cpu::lr::DnaLongOptions>(mapping);
  }

  const ::fa::cpu::lr::DnaLongOptions& long_read() const noexcept {
    if (const auto* rna =
            std::get_if<::fa::cpu::lr::rna::RnaLongOptions>(&mapping)) {
      return rna->base;
    }
    return std::get<::fa::cpu::lr::DnaLongOptions>(mapping);
  }

  ::fa::cpu::lr::rna::RnaLongOptions* rna() noexcept {
    return std::get_if<::fa::cpu::lr::rna::RnaLongOptions>(&mapping);
  }

  const ::fa::cpu::lr::rna::RnaLongOptions* rna() const noexcept {
    return std::get_if<::fa::cpu::lr::rna::RnaLongOptions>(&mapping);
  }
};

}  // namespace fa::cpu::options
