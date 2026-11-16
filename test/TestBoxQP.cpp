#include "Typed.h"
#include <cmath>
#include <limits>
#include "NMPC/BoxQP.h"

/*
The stage QP against enumeration of its active sets. For each of the 3^M ways
of holding each variable free, at its lower bound or at its upper bound, the
minimiser over the free variables is found directly, and the best one that is
inside the bounds is the solution. A positive definite H makes it unique.
*/

namespace {

template <class T, int M>
struct Case {
    NMPC::Matrix<T, M, M> H;
    NMPC::Vector<T, M> g, lower, upper;
};

/* A deterministic pseudo-random sequence, so that the cases are the same in every precision and on every compiler. */
struct Sequence {
    unsigned state = 12345u;
    double next() {
        state = state * 1103515245u + 12345u;
        return double((state >> 8) & 0xFFFFu) / 65535.0 * 2.0 - 1.0;
    }
};

template <class T, int M>
Case<T, M> make_case(Sequence& random) {
    Case<T, M> c;
    NMPC::Matrix<double, M, M> L;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            L(i, j) = random.next();
        }
    }
    const NMPC::Matrix<double, M, M> H = L * L.transpose() + 0.2 * NMPC::Matrix<double, M, M>::Identity();
    c.H = H.template cast<T>();
    for (int i = 0; i < M; i++) {
        c.g(i) = T(3.0 * random.next());
        const double a = random.next(), b = random.next();
        c.lower(i) = T(std::min(a, b) - 0.05);
        c.upper(i) = T(std::max(a, b) + 0.05);
    }
    return c;
}

template <int M>
NMPC::Vector<double, M> enumerate(const Case<double, M>& c) {
    double best = std::numeric_limits<double>::infinity();
    NMPC::Vector<double, M> best_d = NMPC::Vector<double, M>::Zero();
    int combinations = 1;
    for (int i = 0; i < M; i++) {
        combinations *= 3;
    }
    for (int code = 0; code < combinations; code++) {
        /* 0 free, 1 at the lower bound, 2 at the upper bound. */
        int which[M];
        int rest = code;
        for (int i = 0; i < M; i++) {
            which[i] = rest % 3;
            rest /= 3;
        }
        NMPC::Matrix<double, M, M> Hf = c.H;
        NMPC::Vector<double, M> d = NMPC::Vector<double, M>::Zero();
        for (int i = 0; i < M; i++) {
            if (which[i] == 1) d(i) = c.lower(i);
            if (which[i] == 2) d(i) = c.upper(i);
        }
        NMPC::Vector<double, M> rhs = -(c.g + c.H * d);
        for (int i = 0; i < M; i++) {
            if (which[i] != 0) {
                Hf.row(i).setZero();
                Hf.col(i).setZero();
                Hf(i, i) = 1.0;
                rhs(i) = 0.0;
            }
        }
        const NMPC::Vector<double, M> free_part = Hf.llt().solve(rhs);
        for (int i = 0; i < M; i++) {
            if (which[i] == 0) d(i) = free_part(i);
        }
        bool inside = true;
        for (int i = 0; i < M; i++) {
            inside = inside && d(i) >= c.lower(i) - 1e-12 && d(i) <= c.upper(i) + 1e-12;
        }
        const double value = c.g.dot(d) + 0.5 * d.dot(c.H * d);
        if (inside && value < best) {
            best = value;
            best_d = d;
        }
    }
    return best_d;
}

template <class T, int M>
void check_cases(T tolerance) {
    Sequence random;
    for (int n = 0; n < 200; n++) {
        const Case<T, M> c = make_case<T, M>(random);
        Case<double, M> cd;
        cd.H = c.H.template cast<double>();
        cd.g = c.g.template cast<double>();
        cd.lower = c.lower.template cast<double>();
        cd.upper = c.upper.template cast<double>();
        const NMPC::Vector<double, M> expected = enumerate<M>(cd);

        NMPC::Vector<T, M> d = NMPC::Vector<T, M>::Zero();
        Eigen::Array<bool, M, 1> held;
        Eigen::LLT<NMPC::Matrix<T, M, M>> factor;
        const int passes = NMPC::solve_box_qp<T, M>(c.H, c.g, c.lower, c.upper, d, held, factor, 64);
        ASSERT_GT(passes, 0);
        for (int i = 0; i < M; i++) {
            EXPECT_NEAR(double(d(i)), expected(i), double(tolerance)) << "case " << n << ", variable " << i;
            EXPECT_GE(d(i), c.lower(i));
            EXPECT_LE(d(i), c.upper(i));
            /* A held variable is on a bound, and a variable inside the bounds is free. */
            if (held(i)) {
                EXPECT_TRUE(d(i) == c.lower(i) || d(i) == c.upper(i));
            }
        }
    }
}

}

template <typename T>
class BoxQPTest : public ::testing::Test {};

typedef ::testing::Types<float, double> Precisions;
TYPED_TEST_SUITE(BoxQPTest, Precisions);

TYPED_TEST(BoxQPTest, MatchesEnumerationOfActiveSets) {
    const TypeParam tolerance = std::is_same<TypeParam, float>::value ? TypeParam(2e-4) : TypeParam(1e-10);
    check_cases<TypeParam, 1>(tolerance);
    check_cases<TypeParam, 2>(tolerance);
    check_cases<TypeParam, 4>(tolerance);
    check_cases<TypeParam, 5>(tolerance);
}

TYPED_TEST(BoxQPTest, StartsAHairInsideABound) {
    /*
    The case that stopped a projected Newton method short: a variable starts just inside a bound that
    the minimiser lies beyond. The method must still put it on the bound and hold it there.
    */
    typedef TypeParam T;
    NMPC::Matrix<T, 2, 2> H;
    H << T(2), T(0.3), T(0.3), T(1);
    const NMPC::Vector<T, 2> g(T(-5), T(0.2));
    const NMPC::Vector<T, 2> lower(T(-1), T(-1)), upper(T(1), T(1));
    NMPC::Vector<T, 2> d(T(1) - T(1e-5), T(0));
    Eigen::Array<bool, 2, 1> held;
    Eigen::LLT<NMPC::Matrix<T, 2, 2>> factor;
    ASSERT_GT((NMPC::solve_box_qp<T, 2>(H, g, lower, upper, d, held, factor, 32)), 0);
    EXPECT_EQ(d(0), T(1));
    EXPECT_TRUE(held(0));
    EXPECT_FALSE(held(1));
    EXPECT_NEAR(double(d(1)), double((-g(1) - H(1, 0) * d(0)) / H(1, 1)), 1e-5);
}

TYPED_TEST(BoxQPTest, ReportsAFailedFactorisation) {
    typedef TypeParam T;
    NMPC::Matrix<T, 2, 2> H;
    H << T(1), T(2), T(2), T(1);    /* indefinite */
    const NMPC::Vector<T, 2> g(T(1), T(1));
    const NMPC::Vector<T, 2> lower(T(-1), T(-1)), upper(T(1), T(1));
    NMPC::Vector<T, 2> d = NMPC::Vector<T, 2>::Zero();
    Eigen::Array<bool, 2, 1> held;
    Eigen::LLT<NMPC::Matrix<T, 2, 2>> factor;
    EXPECT_EQ((NMPC::solve_box_qp<T, 2>(H, g, lower, upper, d, held, factor, 32)), -1);
}
