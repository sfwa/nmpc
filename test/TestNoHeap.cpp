#include <atomic>
#include <cstdlib>
#include <new>

#include "NMPC/Solver.h"
#include "Problems.h"
#include "Typed.h"

/*
A solve must not touch the heap. Every allocation in the test program is
counted, and the count must not move across a reset, an advance and a solve.
EIGEN_NO_MALLOC is defined for this target as well, which makes Eigen assert on
any allocation of its own in a debug build.
*/

namespace {
std::atomic<long> allocations(0);
}

void* operator new(std::size_t size) {
  allocations++;
  void* p = std::malloc(size == 0 ? 1 : size);
  if (!p) {
    std::abort();
  }
  return p;
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

template <typename T>
class NoHeapTest : public ::testing::Test {};

typedef ::testing::Types<float, double> Precisions;
TYPED_TEST_SUITE(NoHeapTest, Precisions);

TYPED_TEST(NoHeapTest, ASolveAllocatesNothing) {
  typedef TypeParam T;
  typedef Attitude_Problem<T> Problem;
  typedef NMPC::Solver<Problem> Solver;
  typedef Eigen::Matrix<T, 3, 1> V3;
  static typename Solver::Workspace workspace;
  Problem problem;
  typename Problem::State x;
  x.q = NMPC::so3_exp<T>(V3(T(0.6), T(-0.4), T(1.2)));
  x.w = V3(T(0.2), T(0), T(-0.3));
  NMPC::Options<T> options;
  options.max_iterations = 3;

  const long before = allocations.load();
  Solver solver(workspace);
  solver.reset(x, Problem::Control::Zero());
  for (int i = 0; i < 5; i++) {
    solver.advance(problem, 1);
    solver.solve(problem, x, options);
  }
  const long after = allocations.load();
  EXPECT_EQ(after, before);
}

TYPED_TEST(NoHeapTest, ASolveWithPerStageJacobiansAllocatesNothing) {
  typedef TypeParam T;
  typedef Attitude_Problem<T> Problem;
  typedef NMPC::Solver<Problem, Problem::Stages, NMPC::Jacobians::PerStage>
      Solver;
  typedef Eigen::Matrix<T, 3, 1> V3;
  static typename Solver::Workspace workspace;
  Problem problem;
  typename Problem::State x;
  x.q = NMPC::so3_exp<T>(V3(T(0.6), T(-0.4), T(1.2)));
  x.w = V3(T(0.2), T(0), T(-0.3));
  NMPC::Options<T> options;
  options.max_iterations = 3;

  const long before = allocations.load();
  Solver solver(workspace);
  solver.reset(x, Problem::Control::Zero());
  for (int i = 0; i < 5; i++) {
    solver.advance(problem, 1);
    solver.solve(problem, x, options);
  }
  const long after = allocations.load();
  EXPECT_EQ(after, before);
}
