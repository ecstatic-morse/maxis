#ifndef MAXIS_RNG_H
#define MAXIS_RNG_H

#include <random>

namespace maxis {
namespace random {

thread_local static std::mt19937_64 engine{};

// Seeds the thread's random number generator
template <typename Seed> void seed(Seed seed) { engine.seed(seed); }

// Returns an index into
size_t index(size_t start, size_t end);
double probability();

} // namespace random
} // namespace maxis

#endif
