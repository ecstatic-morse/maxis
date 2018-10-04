#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <thread>

#include "boost/graph/adjacency_matrix.hpp"
#include "boost/program_options.hpp"

#include "maxis/genetic.hpp"
#include "maxis/parallel.hpp"
#include "maxis/solver.hpp"

using namespace maxis;

void print_heading() {
  std::cout << "Thread" << '\t' << "Iter" << '\t' << "Max" << '\t' << "Avg"
            << std::endl;
}

auto parse_arguments(int argc, char *argv[]) {
  namespace po = boost::program_options;
  po::options_description opts{"Options"};

  // clang-format off
  opts.add_options()
      ("help,h", "Display usage information")
      ( "file", po::value<std::string>(), "Input file")
      ("population,p", po::value<size_t>(), "Population size")
      ( "target,w", po::value<double>(), "Target fitness value")
      ( "mutation,m", po::value<double>()->default_value(9), "Mutation rate")
      ( "mutation-start", po::value<size_t>()->default_value(10000),
        "(Variable Rate Mutator) Point at which mutation reaches half of its final rate")
      ("mutation-gradient", po::value<double>()->default_value(0.05, "0.05"),
       "(Variable Rate Mutator) Mutation growth rate")
      ( "migration-ratio", po::value<double>()->default_value(0.25, "0.25"),
        "The percentage of the population migrated between threads per migration")
      ( "migration-period", po::value<size_t>()->default_value(4096),
        "The number of genetic operations done before a migration is attempted")
      ( "selector,s", po::value<size_t>()->default_value(4),
        "(Tournament Selector) Tournament Size");
  // clang-format on

  po::positional_options_description p;
  p.add("file", -1);
  po::variables_map vm;
  po::store(
      po::command_line_parser(argc, argv).options(opts).positional(p).run(),
      vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << opts << std::endl;
  }

  return vm;
}

auto build_solver(const Graph &g,
                  const boost::program_options::variables_map &opts) {
  // Population size
  auto pop_size = boost::num_vertices(g) / 2;
  if (opts.count("population")) {
    pop_size = opts["population"].as<size_t>();
  }

  // Mutator
  // auto rate = genetic::FixedMutationRate(opts["mutation"].as<double>());
  auto rate = genetic::AdaptiveMutationRate(
      opts["mutation"].as<double>(), opts["mutation-start"].as<size_t>(),
      opts["mutation-gradient"].as<double>());
  auto mut =
      std::make_unique<genetic::SetResetMutator<decltype(rate)>>(rate, 3.0);

  // Selector
  auto sel = std::make_unique<genetic::TournamentSelector>(
      opts["selector"].as<size_t>());

  // Recombinator
  auto rec = std::make_unique<genetic::BlendingRecombinator>();

  return maxis::GeneticSolver(g, pop_size, std::move(mut), std::move(rec),
                              std::move(sel));
}

int main(int argc, char *argv[]) {
  auto opts = parse_arguments(argc, argv);

  if (opts.count("help"))
    return 0;

  if (!opts.count("file")) {
    std::cerr << "No input file specified" << std::endl;
    return 1;
  }

  // Parse some options for the parallel solver
  auto migration_period = opts["migration-period"].as<size_t>();
  auto migration_ratio = opts["migration-ratio"].as<double>();

  // Parse and initialize graph
  std::ifstream is(opts["file"].as<std::string>());
  if (!is) {
    std::cerr << "Input file not found" << std::endl;
    return 1;
  }
  auto graph = parse_ascii_dimacs(is);

  // A lock around cout and cerr for worker threads.
  std::mutex io_lock;

  // We must reserve the proper number of workers ahead of time so we don't
  // invalidate references as we emplace them.
  auto num_threads = std::thread::hardware_concurrency();
  std::vector<maxis::ParallelGeneticSolver> workers;
  workers.reserve(num_threads);

  // Create as many channels as there are workers.
  // We arrange the workers in a circle, and each worker sends to the clockwise
  // worker and receives from the counter-clockwise worker.
  std::vector<maxis::OwnedQueue> channels{};
  for (int i = 0; i < num_threads; ++i) {
    auto q = std::make_shared<maxis::Queue>(8);
    channels.push_back(q);
  }

  // Initalize a parallel solver for each thread.
  std::vector<std::thread> threads;
  for (int i = 0; i < num_threads; ++i) {
    workers.emplace_back(build_solver(graph, opts), channels[i],
                         channels[(i + 1) % num_threads], migration_ratio,
                         migration_period);
    auto &worker = workers.back();

    // Move each worker onto a new thread and begin iterating.
    threads.emplace_back([i, &worker, &io_lock]() {
      auto max_fitness = std::numeric_limits<double>::lowest();
      for (;;) {
        auto &member = worker.iterate();
        if (member.fitness() > max_fitness) {
          max_fitness = member.fitness();
          std::lock_guard g{io_lock};
          std::cout << "New best: " << max_fitness << '\n';
          std::cout << member << std::endl;
          print_heading();
        }

        if (worker.iterations() % 4096 == 0) {
          auto avg = worker.average_fitness();
          std::lock_guard g{io_lock};
          std::cout << i << '\t' << worker.iterations() << '\t' << max_fitness
                    << '\t' << avg << std::endl;
        }
      }
    });
  }

  // This join will never complete
  for (auto &t : threads) {
    t.join();
  }

  return 0;
}
