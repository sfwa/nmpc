/*
Part of NMPC. Distributed under the MIT licence: see LICENSE in the root of
the repository.
*/

#ifndef NMPC_PROBLEM_H
#define NMPC_PROBLEM_H

#include "NMPC/Types.h"
#include <type_traits>
#include <utility>

/*
The problem concept.

The solver is a template on a problem type. The problem owns its model, its
reference and its weights, and the library never sees a vehicle. The cost is
half the sum of squared residuals, plus an optional term that is linear in
the command, so the Gauss-Newton Hessian needs no second derivatives.

A problem provides:

    using Scalar = T;
    using State   = ...;    any copyable type; it may live on a manifold
    using Control = Eigen::Matrix<T, NU, 1>;
    enum { NX, NU, NR, NT, Stages };
        NX      size of the state's tangent (for a state with a quaternion, one less than its numbers)
        NU      commands
        NR, NT  stage and terminal residuals
        Stages  the plan's stages; a solver may be given another count

    static State retract(const State& x, const Tangent& d);
    static Tangent difference(const State& x, const State& y);
        The only two operations the solver needs on a state, with
        retract(y, difference(x, y)) == x. Types.h has the attitude versions.

    State step(int k, const State& x, const Control& u) const;
        One stage of the model. The problem decides how many model steps a
        stage takes and how long the stage is.

    void step_jacobian(int k, const State& x, const Control& u, A& a, B& b) const;
        The derivatives of step in the tangent: a is NX x NX and b is NX x NU,
        with a(i, j) = d difference(step(retract(x, e_j)), step(x)) / d e_j.
        Analytic. Check.h compares it with central differences in the tests.

    void stage_residual(int k, const State& x, const Control& u, R& r) const;
    void terminal_residual(const State& x, R& r) const;
        Stage residuals carry the square root of the stage length, so that
        weights are per second and do not change with the stage length.

    void residual_jacobians(int k, const State& x, const Control& u, JX& jx, JU& ju) const;
    void terminal_jacobian(const State& x, JX& jx) const;
        Analytic, in the tangent.

    void control_bounds(int k, Control& lower, Control& upper) const;
        Data that may change at every solve: an inner loop's share of the
        actuators comes off it, and so do the limits of a failed actuator.
        The solver reads it once a solve. Use infinities for no bound.

The arguments named A, B, R, JX and JU above are the fixed-size Eigen types of
ProblemTraits below. A problem may take them as those types or as templates.

Optionally, a problem provides:

    Tangent state_scale;  Control control_scale;
        Typical sizes of each state and command. The solver works in units of
        these, which conditions its regularisation and the stage QP in single
        precision. Without them every scale is 1.

    T linear_cost(int k, const State& x, const Control& u) const;
    void linear_cost_gradient(int k, const State& x, const Control& u, Tangent& gx, Control& gu) const;
        A cost term with no Gauss-Newton Hessian, such as contouring's reward
        for progress along a path. Both or neither.

    void cost_expansion(int k, const State& x, const Control& u,
                        Tangent& lx, Control& lu, MatA& Lxx, Gain& Lux, MatU& Luu) const;
    void terminal_expansion(const State& x, Tangent& lx, MatA& Lxx) const;
        The Gauss-Newton expansion of the cost, for a problem that can form it
        faster than the generic product of residual Jacobians by knowing their
        structure. It includes the linear term. Only the upper triangle of
        Lxx is read. In the problem's own units: the solver scales it.

    void tail_control(int k, Control& u) const;
        The command for a stage that enters the plan when it is shifted, such
        as the reference's command for that stage. Without it the entering
        stage repeats the last command, which made a 10 s horizon fail in tight
        turns on the fixed wing and is fine over the multirotor's 0.6 s.

    static bool finite(const State& x);
        Required when State is not an Eigen type.

Context such as a wind, a thrust scale or a force error is a member of the
problem that the application sets between solves. Weights are members too, so
a cost can change between solves: a degraded mode drops what can no longer be
held. The library has no context or weight type of its own.
*/

namespace NMPC {

template <class Problem>
struct ProblemTraits {
    using T = typename Problem::Scalar;
    using State = typename Problem::State;
    using Control = typename Problem::Control;
    enum {
        NX = Problem::NX,
        NU = Problem::NU,
        NR = Problem::NR,
        NT = Problem::NT
    };
    using Tangent = Vector<T, NX>;
    using MatA = Matrix<T, NX, NX>;
    using MatB = Matrix<T, NX, NU>;
    using Gain = Matrix<T, NU, NX>;
    using MatU = Matrix<T, NU, NU>;
    using StageResidual = Vector<T, NR>;
    using TerminalResidual = Vector<T, NT>;
    using StageJx = Matrix<T, NR, NX>;
    using StageJu = Matrix<T, NR, NU>;
    using TerminalJx = Matrix<T, NT, NX>;

    static_assert(NX > 0 && NU > 0, "A problem needs a state and a command");
    static_assert(NR >= 0 && NT >= 0, "Residual counts cannot be negative");
    static_assert(std::is_same<Control, Vector<T, NU>>::value,
        "Problem::Control must be Eigen::Matrix<Scalar, NU, 1>");
};

namespace detail {

template <class P> using Tr = ProblemTraits<P>;

template <class P, class = void>
struct HasStateScale : std::false_type {};
template <class P>
struct HasStateScale<P, std::void_t<decltype(typename Tr<P>::Tangent(std::declval<const P&>().state_scale))>>
    : std::true_type {};

template <class P, class = void>
struct HasControlScale : std::false_type {};
template <class P>
struct HasControlScale<P, std::void_t<decltype(typename Tr<P>::Control(std::declval<const P&>().control_scale))>>
    : std::true_type {};

template <class P, class = void>
struct HasLinearCost : std::false_type {};
template <class P>
struct HasLinearCost<P, std::void_t<
        decltype(std::declval<const P&>().linear_cost(0, std::declval<const typename Tr<P>::State&>(),
            std::declval<const typename Tr<P>::Control&>())),
        decltype(std::declval<const P&>().linear_cost_gradient(0, std::declval<const typename Tr<P>::State&>(),
            std::declval<const typename Tr<P>::Control&>(), std::declval<typename Tr<P>::Tangent&>(),
            std::declval<typename Tr<P>::Control&>()))>>
    : std::true_type {};

template <class P, class = void>
struct HasCostExpansion : std::false_type {};
template <class P>
struct HasCostExpansion<P, std::void_t<decltype(std::declval<const P&>().cost_expansion(0,
        std::declval<const typename Tr<P>::State&>(), std::declval<const typename Tr<P>::Control&>(),
        std::declval<typename Tr<P>::Tangent&>(), std::declval<typename Tr<P>::Control&>(),
        std::declval<typename Tr<P>::MatA&>(), std::declval<typename Tr<P>::Gain&>(),
        std::declval<typename Tr<P>::MatU&>()))>>
    : std::true_type {};

template <class P, class = void>
struct HasTerminalExpansion : std::false_type {};
template <class P>
struct HasTerminalExpansion<P, std::void_t<decltype(std::declval<const P&>().terminal_expansion(
        std::declval<const typename Tr<P>::State&>(), std::declval<typename Tr<P>::Tangent&>(),
        std::declval<typename Tr<P>::MatA&>()))>>
    : std::true_type {};

template <class P, class = void>
struct HasTailControl : std::false_type {};
template <class P>
struct HasTailControl<P, std::void_t<decltype(std::declval<const P&>().tail_control(0,
        std::declval<typename Tr<P>::Control&>()))>>
    : std::true_type {};

template <class P, class = void>
struct HasFinite : std::false_type {};
template <class P>
struct HasFinite<P, std::void_t<decltype(P::finite(std::declval<const typename Tr<P>::State&>()))>>
    : std::true_type {};

template <class S, class = void>
struct IsEigen : std::false_type {};
template <class S>
struct IsEigen<S, std::void_t<decltype(std::declval<const S&>().allFinite())>> : std::true_type {};

/* The dispatch of the optional parts of the concept. */

template <class P>
typename Tr<P>::Tangent state_scale(const P& problem) {
    if constexpr (HasStateScale<P>::value) {
        return typename Tr<P>::Tangent(problem.state_scale);
    } else {
        (void)problem;
        return Tr<P>::Tangent::Ones();
    }
}

template <class P>
typename Tr<P>::Control control_scale(const P& problem) {
    if constexpr (HasControlScale<P>::value) {
        return typename Tr<P>::Control(problem.control_scale);
    } else {
        (void)problem;
        return Tr<P>::Control::Ones();
    }
}

template <class P>
bool finite(const typename Tr<P>::State& x) {
    if constexpr (HasFinite<P>::value) {
        return P::finite(x);
    } else {
        static_assert(IsEigen<typename Tr<P>::State>::value,
            "A problem whose State is not an Eigen type must provide static bool finite(const State&)");
        return x.allFinite();
    }
}

/* Half the squared residuals, and the linear term where there is one. */
template <class P>
typename Tr<P>::T stage_cost(const P& problem, int k, const typename Tr<P>::State& x,
        const typename Tr<P>::Control& u) {
    using T = typename Tr<P>::T;
    T cost = T(0);
    if constexpr (Tr<P>::NR > 0) {
        typename Tr<P>::StageResidual r;
        problem.stage_residual(k, x, u, r);
        cost = T(0.5) * r.squaredNorm();
    }
    if constexpr (HasLinearCost<P>::value) {
        cost += problem.linear_cost(k, x, u);
    }
    return cost;
}

template <class P>
typename Tr<P>::T terminal_cost(const P& problem, const typename Tr<P>::State& x) {
    using T = typename Tr<P>::T;
    if constexpr (Tr<P>::NT > 0) {
        typename Tr<P>::TerminalResidual r;
        problem.terminal_residual(x, r);
        return T(0.5) * r.squaredNorm();
    } else {
        (void)problem;
        (void)x;
        return T(0);
    }
}

}

}

#endif
