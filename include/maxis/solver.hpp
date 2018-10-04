#ifndef MAXIS_SOLVER_H
#define MAXIS_SOLVER_H

#include <memory>
#include <unordered_set>
#include <vector>

#include "boost/graph/adjacency_list.hpp"

#include "./bitvec.hpp"

namespace maxis {

using Weight = double;

using Graph =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Weight>;
// using Graph = boost::adjacency_matrix<boost::undirectedS, Weight>;

// Parses a test graph in ASCII (not binary) DIMACS format
Graph parse_ascii_dimacs(std::istream &input);

// A genome along with its fitness store
struct Member {
public: // Constructors
  Member() {}

  // Computes the fitness value of the given bit vector.
  Member(std::unique_ptr<BitVec> bv, const Graph &g);

  // Creates a random, valid independent set.
  static Member random(const Graph &g);

public: // Getters
  Weight fitness() const { return fitness_; }
  const BitVec &genome() const { return *genome_; }

public:
  // Applies the functor to the genome, updating its fitness value
  template <typename F> Weight modify(F functor, const Graph &g) {
    modify_skip_fitness_update(functor);
    recompute_fitness(g);
  }

  template <typename F> void modify_skip_fitness_update(F functor) {
    functor(*genome_);
  }

  // Updates the fitness of the genome
  void recompute_fitness(const Graph &g);

  // Returns true if the member contains a valid independent set on the given
  // graph. Should only be necessary for testing.
  bool is_valid(const Graph &g) const;

private:
  std::unique_ptr<BitVec> genome_;
  Weight fitness_;
};

std::ostream &operator<<(std::ostream &os, const Member &m);

// Forward declare the various genetic operators.
namespace genetic {
class Mutator;
class Recombinator;
class Selector;
} // namespace genetic

class GeneticSolver {
private:
  // The hash table used to test for duplicate genomes.
  using HashTable =
      std::unordered_set<std::reference_wrapper<const BitVec>,
                         std::hash<BitVec>, std::equal_to<BitVec>>;

public: // Constructors
  GeneticSolver(const Graph &g, size_t pop_size,
                std::unique_ptr<genetic::Mutator> mutator,
                std::unique_ptr<genetic::Recombinator> recombinator,
                std::unique_ptr<genetic::Selector> selector);

public: // Getters
  uint64_t iterations() const { return iterations_; }
  size_t size() const { return population_.size(); }

  double average_fitness() const { return total_fitness_ / size(); }

public:
  Member &iterate();

  // Moves {n} members into {dst}, removing them from the population.
  void migrate_to(size_t n, std::vector<Member> &dst);

  // Moves all the members from {src} into the population, skipping any
  // duplicates. Returns the number of members which were successfully added to
  // the population. After this is called {src} will be empty.
  size_t migrate_from(std::vector<Member> &src);

private:
  bool insert_member(Member m);
  void local_search(BitVec &bv);

private:
  const Graph &g_;

  // Genetic operators
  std::unique_ptr<genetic::Mutator> mutator_;
  std::unique_ptr<genetic::Recombinator> recombinator_;
  std::unique_ptr<genetic::Selector> selector_;

  std::vector<Member> population_;
  HashTable duplicates_;

  Weight total_fitness_;
  size_t iterations_;

  // The set of vertices in the graph, sorted by weight / out_degree.
  std::vector<size_t> promising_vertices_;
};

} // namespace maxis

#endif // MAXIS_SOLVER_H
