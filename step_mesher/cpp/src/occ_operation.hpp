#pragma once

#include "step_mesher/occ_context.hpp"
#include <optional>
#include <utility>

namespace stepmesher::detail {

// The former executable returned a nonzero exit status when an OCC operation
// threw; its caller then used the original numerical fallback. Keep that
// boundary for the six extracted operations, but never hide native identity
// or runtime-contract failures introduced by the in-process bridge.
template <class Operation>
auto call_legacy_occ_operation(Operation &&operation)
    -> std::optional<decltype(std::forward<Operation>(operation)())> {
  try {
    return std::forward<Operation>(operation)();
  } catch (const occ::NativeContractError &) {
    throw;
  } catch (...) {
    return std::nullopt;
  }
}

} // namespace stepmesher::detail
