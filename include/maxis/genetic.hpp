#ifndef MAXIS_GENETIC_H
#define MAXIS_GENETIC_H

#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

#include "boost/dynamic_bitset.hpp"

#include "./bitvec.hpp"
#include "./random.hpp"
#include "./solver.hpp"

namespace maxis {
namespace genetic {

using Genome = maxis::BitVec;
using Member = maxis::Member;

class Mutator {
public:
  virtual ~Mutator() {}
  virtual void mutate(Genome &) = 0;
  virtual void update(uint64_t iterations) = 0;
};

class Recombinator {
public:
  virtual ~Recombinator() {}
  virtual void breed(const Member &, const Member &, Genome &) = 0;
};

class Selector {
public:
  virtual ~Selector() {}
  virtual size_t select(const std::vector<Member> &population,
                        double total_fitness) = 0;
};

// Mutators
// ===========================================================================

class FixedMutationRate {
public:
  FixedMutationRate(double rate) : rate_{rate} {}

  double average_mutations_per_genome() const { return rate_; }
  void update(uint64_t) {}

private:
  const double rate_;
};

// The adaptive mutator uses a sigmoid function to set the mutation rate.
// {ss} is the average number of flipped bits as the number of iterations
// approaches infinity. {halfway} is the number of iterations before the
// mutation rate reaches half of its steady state value. {slope} is the
// slope at the halfway point.
class AdaptiveMutationRate {
public:
  AdaptiveMutationRate(double rate, size_t halfway, double slope)
      : steady_state_{rate}, halfway_point_{static_cast<double>(halfway)},
        slope_halfway_{slope} {}

public:
  double average_mutations_per_genome() const { return rate_; }

  void update(uint64_t iterations) {
    double steepness = 4 * slope_halfway_;

    // S(x) = (1 + e^-x)^-1
    rate_ = exp(-steepness * (iterations - halfway_point_) / halfway_point_);
    rate_ = steady_state_ / (1 + rate_);
  }

private:
  double rate_;

  const double steady_state_;
  const double halfway_point_;
  const double slope_halfway_;
};

template <typename Rate> class SimpleMutator : public virtual Mutator {
public:
  SimpleMutator(Rate r) : rate_{r} {}

  virtual void mutate(Genome &gene) override {
    double threshold = rate_.average_mutations_per_genome() / gene.size();
    FOR_EACH_BV(v, gene) {
      if (random::probability() < threshold) {
        gene[v].flip();
      }
    }
  }

  virtual void update(uint64_t iterations) override {
    rate_.update(iterations);
  }

private:
  Rate rate_;
};

// A mutator where the chance of setting a bit is a fixed multiple of the
// chance of unsetting a bit, regardless of the ratio of set bits.
template <typename Rate> class SetResetMutator : public virtual Mutator {
public:
  SetResetMutator(Rate r, double unset_mul = 0.5)
      : rate_{r}, unset_mul_{unset_mul} {}

  // avg_mutations = p_flip_if_set * set_bits + p_flip_if_unset * unset_bits
  virtual void mutate(Genome &gene) override {
    auto set_bits = gene.count();
    auto unset_bits = gene.size() - set_bits;

    double set_threshold = rate_.average_mutations_per_genome() /
                           (set_bits + unset_mul_ * unset_bits);
    double unset_threshold = unset_mul_ * set_threshold;

    FOR_EACH_BV(v, gene) {
      auto threshold = gene[v] ? set_threshold : unset_threshold;
      if (random::probability() < threshold) {
        gene[v].flip();
      }
    }
  }

  virtual void update(uint64_t iterations) override {
    rate_.update(iterations);
  }

private:
  Rate rate_;

  // The likelihood of an unset bit flipping relative to an set bit.
  // 1.0 means a flip is equally likely.
  double unset_mul_;
};

// Recombinators
// ===========================================================================

class BlendingRecombinator : public virtual Recombinator {
public:
  virtual void breed(const Member &a, const Member &b, Genome &out) override {
    assert(a.genome().size() == a.genome().size());

    // If the two genomes differ, use b's bit with a probability proportional to
    // its fitness.
    auto f1 = a.fitness();
    auto f2 = b.fitness();
    double threshold = f2 / (f1 + f2);

    for (size_t i = 0; i < a.genome().size(); ++i) {
      if (a.genome()[i] == b.genome()[i] || random::probability() < threshold) {
        out[i] = b.genome()[i];
      } else {
        out[i] = a.genome()[i];
      }
    }
  }
};

// Selectors
// ===========================================================================

// The tournament selector picks `size` cantidates out of the population
// at random and returns the best one.
class TournamentSelector : public virtual Selector {
public:
  TournamentSelector(size_t size) : size_{size} {}

  virtual size_t select(const std::vector<Member> &population,
                        double total_fitness) override {
    auto max_fitness = std::numeric_limits<double>::lowest();
    auto best = -1;
    for (size_t i = 0; i < size_; ++i) {
      auto idx = static_cast<int>(random::index(0, population.size()));
      auto fitness = population[idx].fitness();
      if (fitness > max_fitness) {
        max_fitness = fitness;
        best = idx;
      }
    }

    assert(best != -1);
    return best;
  }

private:
  size_t size_;
};

} // namespace genetic
} // namespace maxis

#endif // MAXIS_GENETIC_H
