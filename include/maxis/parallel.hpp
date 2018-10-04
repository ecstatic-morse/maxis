#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>
#include <vector>

#include "maxis/solver.hpp"
#include "moodycamel/readerwriterqueue.h"

namespace maxis {

using Queue = ::moodycamel::BlockingReaderWriterQueue<std::vector<Member>>;
using OwnedQueue = std::shared_ptr<Queue>;

class ParallelGeneticSolver {
public:
  ParallelGeneticSolver(GeneticSolver &&solver, const OwnedQueue &send,
                        const OwnedQueue &recv, double migration_ratio,
                        unsigned migration_period)
      : solver_{std::move(solver)}, send_{send}, recv_{recv},
        migration_size_{static_cast<size_t>(migration_ratio * solver.size())},
        migration_period_{static_cast<int64_t>(migration_period)},
        migration_timer_{migration_period_} {
    assert(migration_ratio > 0.0 && migration_ratio < 1.0);
    buf_.reserve(migration_period);
  }

public:
  Member &iterate() {
    migration_timer_ -= 1;

    // It's time to send some of our population to another thread.
    if (migration_timer_ == 0) {
      solver_.migrate_to(migration_size_, buf_);
      assert(send_->try_enqueue(std::move(buf_)));
    }

    // A negative migration timer indicates that we've sent data but not yet
    // received any.
    if (migration_timer_ < 0) {
      // Check to see if our partner thread has migrated some data.
      auto recv_ok = recv_->try_dequeue(buf_);

      // If the scheduler is unfair, or we have scheduled more solvers than CPU
      // cores, it is likely that we will starve the thread which is supposed
      // to migrate its population to us. This is won't cause incorrectness,
      // as we will simply continue running on a smaller subset of the
      // population, but results in a less diverse population.
      /*
      if (!recv_ok && migration_timer_ <= 10*migration_period_) {
        recv_->wait_dequeue(buf_);
        recv_ok = true;
      }
      */

      // Reset the timer once we've received data
      if (recv_ok) {
        migration_timer_ = migration_period_;
        solver_.migrate_from(buf_);
      }
    }

    auto &member = solver_.iterate();
    return member;
  }

  int64_t iterations() const { return solver_.iterations(); }

  double average_fitness() const { return solver_.average_fitness(); }

private:
  GeneticSolver solver_;
  OwnedQueue send_;
  OwnedQueue recv_;

  const size_t migration_size_;
  const int64_t migration_period_;
  int64_t migration_timer_;
  std::vector<Member> buf_;
};

} // namespace maxis
