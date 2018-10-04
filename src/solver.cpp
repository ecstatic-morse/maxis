#include <algorithm>
#include <cstdint>
#include <numeric>
#include <vector>

#include "maxis/bitvec.hpp"
#include "maxis/genetic.hpp"
#include "maxis/solver.hpp"

namespace maxis {

// Finds the index of the nth bit in the bitset, or npos if there aren't that
// many bits.
ssize_t find_nth_bit(const BitVec &bv, size_t n) {
  FOR_EACH_BV(v, bv) {
    if (n == 0) {
      return v;
    }

    n -= 1;
  }

  return boost::dynamic_bitset<>::npos;
}

bool is_independent_set(const BitVec &bv, const Graph &g) {
  FOR_EACH_BV(v, bv) {
    for (auto [e, end] = boost::out_edges(v, g); e != end; ++e) {
      auto t = boost::target(*e, g);
      if (bv[t]) {
        return false;
      }
    }
  }

  return true;
}

Weight fitness(const BitVec &bv, const Graph &g) {
  Weight total = 0.0;
  FOR_EACH_BV(v, bv) { total += g[v]; }
  return total;
}

void Member::recompute_fitness(const Graph &g) {
  fitness_ = maxis::fitness(*genome_, g);
}

// Creates a random independent set on the graph.
Member Member::random(const Graph &g) {
  auto order = boost::num_vertices(g);
  ssize_t uncovered = order;

  // At the beginning, all vertices can be chosen.
  BitVec candidates(order);
  candidates.flip();

  auto bv = std::make_unique<BitVec>(order);

  while (uncovered > 0) {
    // Add a random vertex to the set
    auto v = find_nth_bit(candidates, random::index(0, uncovered));
    (*bv)[v] = true;

    // Remove it and its neighbors from the candidates set
    candidates[v] = false;
    --uncovered;
    for (auto [e, end] = boost::out_edges(v, g); e != end; ++e) {
      auto u = boost::target(*e, g);
      if (candidates[u]) {
        candidates[u] = false;
        --uncovered;
      }
    }
  }

  assert(is_independent_set(*bv, g));
  return Member(std::move(bv), g);
}

bool Member::is_valid(const Graph &g) const {
  return is_independent_set(genome(), g);
}

Member::Member(std::unique_ptr<BitVec> bv, const Graph &g)
    : genome_{std::move(bv)} {
  recompute_fitness(g);
}

std::ostream &operator<<(std::ostream &os, const Member &m) {
  os << "{";
  FOR_EACH_BV(v, m.genome()) { os << v << ", "; }
  os << "}";
  return os;
}

// Ripple carry increment for a collection of bits.
// Returns false on overflow.
bool increment_bit_vector(BitVec &bv) {
  bool carry = true;
  FOR_EACH_BV(v, bv) {
    if (bv[v] == 0) {
      bv[v] = 1;
      carry = false;
      break;
    } else {
      bv[v] = 0;
    }
  }

  return !carry;
}

BitVec solve_brute_force(const Graph &graph) {
  int max_weight = 0;
  BitVec max_set;

  BitVec bv(boost::num_vertices(graph), 0);
  do {
    if (is_independent_set(bv, graph)) {
      auto weight = fitness(bv, graph);
      if (weight > max_weight) {
        max_weight = weight;
        max_set = bv;
      }
    }
  } while (increment_bit_vector(bv));

  return max_set;
}

GeneticSolver::GeneticSolver(
    const Graph &old_g, size_t pop_size, std::unique_ptr<genetic::Mutator> mutator,
    std::unique_ptr<genetic::Recombinator> recombinator,
    std::unique_ptr<genetic::Selector> selector)
    : g_{0}, mutator_{std::move(mutator)},
      recombinator_{std::move(recombinator)}, selector_{std::move(selector)},
      total_fitness_{0.0}, iterations_{0} {

  auto order = boost::num_vertices(old_g);

  // Sort vertices by weight / out_degree
  std::vector promising_vertices(order, 0);
  std::iota(promising_vertices.begin(), promising_vertices.end(), 0);
  std::sort(promising_vertices.begin(), promising_vertices.end(),
            [&old_g](const auto &a, const auto &b) {
              return old_g[a] / boost::out_degree(a, old_g) >
                     old_g[b] / boost::out_degree(b, old_g);
            });

  // Map the vertex indices in the input graph to the vertices in the sorted graph
  std::vector old_vertices(order, 0);
  for (int i = 0; i < order; ++i) {
    old_vertices[promising_vertices[i]] = i;
  }

  g_ = Graph(order);
  for (int i = 0; i < order; ++i) {
    auto v = promising_vertices[i];
    g_[i] = old_g[v]; // Copy the vertex weight

    // Copy the edges
    for (auto [e, end] = boost::out_edges(v, old_g); e != end; ++e) {
      auto u = boost::target(*e, old_g);
      boost::add_edge(i, old_vertices[u], g_);
    }
  }

  // Populate with random genomes
  for (int i = 0; i < pop_size; ++i) {
    for (;;) {
      if (insert_member(Member::random(g_))) {
        break;
      }
    }
  }
}

Member &GeneticSolver::iterate() {
  // Periodically update the mutator so it can adaptively select the mutation
  // rate.
  const uint64_t UPDATE_PERIOD = 128;
  if (iterations_ % UPDATE_PERIOD == 0) {
    mutator_->update(iterations_);
  }

  // Find the worst genome out of 10 randomly chosen to replace
  auto min_fitness = std::numeric_limits<double>::max();
  ssize_t worst = -1;
  for (int i = 0; i < 20; ++i) {
    auto idx = static_cast<ssize_t>(random::index(0, population_.size()));
    auto fitness = population_[idx].fitness();
    if (fitness < min_fitness) {
      min_fitness = fitness;
      worst = idx;
    }
  }

  // The member to be replaced
  auto &member = population_[worst];
  duplicates_.erase(std::ref(member.genome()));
  total_fitness_ -= member.fitness();

  for (;;) {
    auto a = selector_->select(population_, total_fitness_);
    auto b = selector_->select(population_, total_fitness_);

    member.modify_skip_fitness_update([this, a, b](BitVec &genome) {
      recombinator_->breed(population_[a], population_[b], genome);
      mutator_->mutate(genome);
      greedy_local_search(genome);
    });

    auto [_, ok] = duplicates_.emplace(member.genome());
    if (!ok) {
      continue;
    }

    member.recompute_fitness(g_);
    total_fitness_ += member.fitness();
    break;
  }

  iterations_ += 1;
  return member;
}

// Adds a new member to the set, failing if it already exists in the duplicate
// table. Used only during initialization.
bool GeneticSolver::insert_member(Member m) {
  auto [_, ok] = duplicates_.emplace(m.genome());
  if (!ok) {
    return false;
  }

  population_.emplace_back(std::move(m));
  total_fitness_ += population_.back().fitness();
  return true;
}

// Implements the heuristic feasability described by Beasley and Chu.  The
// operator ensures that genotypes created from mutation and breeding are
// valid independent sets. It is essentially a greedy solver. It makes a
// higher mutation rate necessary because some of the mutated bits will be
// immediately cancelled out to make the set independent.
//
// Currently it is the performance bottleneck for the genetic algorithm as it
// runs in O(|V|^2) time where |V| is the order of the graph.
void GeneticSolver::greedy_local_search(BitVec &bv) {
  // Exclude the neighbors of all bits set in the bitvec, starting at the most
  // promising and ending with the least.
  FOR_EACH_BV(v, bv) {
    for (auto [e, end] = boost::out_edges(v, g_); e != end; ++e) {
      auto u = boost::target(*e, g_);
      bv.reset(u);
    }
  }

  // Invert the bit vector so that it now contains all UNSET vertices.
  // (This would not be necessary if we could iterate over unset bits)
  bv.flip();

  // Add back vertices whch have no neighbors in the bitset.
  FOR_EACH_BV(v, bv) {
    auto [e, end] = boost::out_edges(v, g_);
    auto has_neighbors = std::any_of(e, end, [this, &bv](auto e) {
      auto u = boost::target(e, g_);
      return !bv[u];
    });

    if (!has_neighbors) {
      bv.reset(v);
    }
  }

  // Restore the original semantics of the bitset.
  bv.flip();
}

size_t GeneticSolver::migrate_from(std::vector<Member> &src) {
  population_.reserve(src.size() + population_.size());

  auto it = std::make_move_iterator(src.rbegin());
  auto end = std::make_move_iterator(src.rend());
  auto transferred = 0;

  for (; it != end; ++it) {
    auto [_, ok] = duplicates_.insert(std::ref(it->genome()));
    if (!ok) {
      continue;
    }

    transferred += 1;
    total_fitness_ += it->fitness();
    population_.emplace_back(*it);
  }

  src.clear();
  return transferred;
}

void GeneticSolver::migrate_to(size_t n, std::vector<Member> &dst) {
  dst.reserve(n);

  // Move the last `n` memebers of the population into dst, erasing them from
  // the duplicates table and updating the total fitness.
  auto it = std::make_move_iterator(population_.rbegin());
  auto end = it + n;
  for (; it != end; ++it) {
    duplicates_.erase(std::ref(it->genome()));
    total_fitness_ -= it->fitness();
    dst.emplace_back(*it);
  }

  // Truncate the population, calling the destructors of the newly moved from
  // members
  population_.resize(population_.size() - n);
}

// DIMACS parsing

using ParseError = std::runtime_error;

Graph parse_ascii_dimacs(std::istream &input) {
  std::string line, token;
  char c;
  size_t n1, n2;
  bool is_initialized = false;

  Graph g(0);

  while (std::getline(input, line)) {
    if (line.empty()) {
      continue;
    }

    std::istringstream ss{line};

    ss >> c;
    switch (c) {
    case 'p': // p FORMAT NODES EDGES (FORMAT=edge)
      ss >> token;
      if (token != "edge")
        throw ParseError("'p' must have FORMAT of 'edge'");

      ss >> n1; // Number of nodes
      ss >> n2; // Number of edges

      g = Graph(n1);
      is_initialized = true;

      // Set all node weights to 1 by default
      for (int i = 0; i < n1; ++i) {
        g[i] = 1.0;
      }

      break;
    case 'e': // e W V (edge connecting w and v)
      if (!is_initialized)
        throw ParseError("'p' must be defined before edge definitions");
      ss >> n1;
      ss >> n2;
      boost::add_edge(n1 - 1, n2 - 1, g);
      break;
    case 'n': // n ID VALUE (weight of given node)
      if (!is_initialized)
        throw ParseError("'p' must be defined before node definitions");

      ss >> n1; // nodeid
      ss >> n2; // weight

      g[n1 - 1] = n2;
      break;
    case 'c':
    default:
      continue;
    }
  }

  if (!is_initialized)
    throw ParseError("'p' definition not found");
  return g;
}

} // namespace maxis
