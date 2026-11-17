#include "NMPC/Check.h"
#include "Problems.h"
#include "Typed.h"

/*
Every test problem's analytic derivatives against central differences, and the
checks themselves against a deliberately wrong Jacobian.
*/

template <typename T>
class CheckTest : public ::testing::Test {};

typedef ::testing::Types<float, double> Precisions;
TYPED_TEST_SUITE(CheckTest, Precisions);

namespace {
template <class T>
T tolerance() {
  return std::is_same<T, float>::value ? T(5e-3) : T(1e-7);
}
}  // namespace

TYPED_TEST(CheckTest, StopProblem) {
  typedef TypeParam T;
  Stop_Problem<T> problem;
  const typename Stop_Problem<T>::State x(T(0.7), T(-1.3));
  typename Stop_Problem<T>::Control u;
  u(0) = T(2);
  EXPECT_TRUE(
      NMPC::check_step_jacobian(problem, 0, x, u, tolerance<T>()).passed);
  EXPECT_TRUE(
      NMPC::check_residual_jacobians(problem, 0, x, u, tolerance<T>()).passed);
  EXPECT_TRUE(NMPC::check_terminal_jacobian(problem, x, tolerance<T>()).passed);
}

TYPED_TEST(CheckTest, LinearProblem) {
  typedef TypeParam T;
  Linear_Problem<T> problem;
  const typename Linear_Problem<T>::State x(T(0.2), T(-0.1), T(0.4));
  const typename Linear_Problem<T>::Control u(T(0.3), T(-0.6));
  EXPECT_TRUE(
      NMPC::check_step_jacobian(problem, 0, x, u, tolerance<T>()).passed);
  EXPECT_TRUE(
      NMPC::check_residual_jacobians(problem, 0, x, u, tolerance<T>()).passed);
  EXPECT_TRUE(NMPC::check_terminal_jacobian(problem, x, tolerance<T>()).passed);
}

TYPED_TEST(CheckTest, AttitudeProblemOnTheManifold) {
  typedef TypeParam T;
  typedef Eigen::Matrix<T, 3, 1> V3;
  Attitude_Problem<T> problem;
  typename Attitude_Problem<T>::State x;
  x.q = NMPC::so3_exp<T>(V3(T(0.3), T(-0.2), T(0.25)));
  x.w = V3(T(1.5), T(-0.7), T(2.0));
  const typename Attitude_Problem<T>::Control u(T(0.01), T(-0.03), T(0.02));
  const NMPC::CheckResult<T> step =
      NMPC::check_step_jacobian(problem, 0, x, u, tolerance<T>());
  EXPECT_TRUE(step.passed) << "worst " << step.error << " at " << step.row
                           << ", " << step.column;
  /*
  The attitude residual's Jacobian uses the inverse left Jacobian of SO(3) to
  second order, whose error grows as the fourth power of the attitude error:
  about 3e-5 here.
  */
  const T looser = std::is_same<T, float>::value ? T(5e-3) : T(2e-4);
  const NMPC::CheckResult<T> residual =
      NMPC::check_residual_jacobians(problem, 0, x, u, looser);
  EXPECT_TRUE(residual.passed) << "worst " << residual.error << " at "
                               << residual.row << ", " << residual.column;
  EXPECT_TRUE(NMPC::check_terminal_jacobian(problem, x, looser).passed);
}

namespace {
template <class T>
struct Wrong_Problem : Stop_Problem<T> {
  template <class A, class B>
  void step_jacobian(int k, const typename Stop_Problem<T>::State& x,
                     const typename Stop_Problem<T>::Control& u, A& a,
                     B& b) const {
    Stop_Problem<T>::step_jacobian(k, x, u, a, b);
    a(0, 1) *= T(1.1);
  }
};
}  // namespace

TYPED_TEST(CheckTest, FindsAWrongJacobian) {
  typedef TypeParam T;
  Wrong_Problem<T> problem;
  const NMPC::CheckResult<T> result = NMPC::check_step_jacobian(
      problem, 0, typename Stop_Problem<T>::State(T(0), T(1)),
      Stop_Problem<T>::Control::Zero(), tolerance<T>());
  EXPECT_FALSE(result.passed);
  EXPECT_EQ(result.row, 0);
  EXPECT_EQ(result.column, 1);
}
