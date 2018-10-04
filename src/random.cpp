#include <random>

#include "maxis/random.hpp"

namespace maxis {
namespace random {

// Returns a random integer in [start, end)
size_t index(size_t start, size_t end) {
  std::uniform_int_distribution<size_t> range{start, end - 1};
  return range(engine);
}

double probability() {
  static std::uniform_real_distribution<double> prob_dist{0.0, 1.0};
  return prob_dist(engine);
}

} // namespace random
} // namespace maxis
