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

TYPED_TEST(BoxQPTest, EqualBoundsAreAlwaysHeld) {
    /*
    A variable whose bounds meet cannot move, so it is held and gets no feedback even where its gradient
    is zero. It was left free when its gradient was zero, and the gain solved with it was not zero.
    */
    typedef TypeParam T;
    NMPC::Matrix<T, 2, 2> H;
    H << T(2), T(0.5), T(0.5), T(1);
    const NMPC::Vector<T, 2> g = NMPC::Vector<T, 2>::Zero();
    const NMPC::Vector<T, 2> lower(T(0), T(-1)), upper(T(0), T(1));
    NMPC::Vector<T, 2> d = NMPC::Vector<T, 2>::Zero();
    Eigen::Array<bool, 2, 1> held;
    Eigen::LLT<NMPC::Matrix<T, 2, 2>> factor;
    ASSERT_GT((NMPC::solve_box_qp<T, 2>(H, g, lower, upper, d, held, factor, 32)), 0);
    EXPECT_TRUE(held(0));
    EXPECT_FALSE(held(1));
    EXPECT_EQ(d(0), T(0));
    /* The factor has the held variable's row and column replaced by the identity. */
    const NMPC::Vector<T, 2> rhs(T(0), T(1));
    const NMPC::Vector<T, 2> solved = factor.solve(rhs);
    EXPECT_EQ(solved(0), T(0));
    EXPECT_NEAR(double(solved(1)), 1.0, 1e-6);

    /* Also from a start with a gradient that would let it go. */
    const NMPC::Vector<T, 2> pulling(T(-3), T(0.5));
    d = NMPC::Vector<T, 2>::Zero();
    ASSERT_GT((NMPC::solve_box_qp<T, 2>(H, pulling, lower, upper, d, held, factor, 32)), 0);
    EXPECT_TRUE(held(0));
    EXPECT_EQ(d(0), T(0));
}

TYPED_TEST(BoxQPTest, AnswerDoesNotDependOnTheScale) {
    /*
    The stop is relative: the same problem multiplied through by a small or a large number has the same
    answer. With an absolute test on the gradient, 0.5 (1e-4 (u - 1))^2 stopped at its start, u = 0.
    */
    typedef TypeParam T;
    for (const double scale : {1e-12, 1e-8, 1e-4, 1.0, 1e4, 1e8}) {
        NMPC::Matrix<T, 1, 1> H;
        H << T(scale);
        NMPC::Vector<T, 1> g, lower, upper, d;
        g << T(-scale);
        lower << T(-10);
        upper << T(10);
        d << T(0);
        Eigen::Array<bool, 1, 1> held;
        Eigen::LLT<NMPC::Matrix<T, 1, 1>> factor;
        ASSERT_GT((NMPC::solve_box_qp<T, 1>(H, g, lower, upper, d, held, factor, 32)), 0) << scale;
        EXPECT_NEAR(double(d(0)), 1.0, 1e-5) << scale;
        EXPECT_FALSE(held(0)) << scale;
    }
    /* The random cases multiplied through: the same answers as unscaled. */
    Sequence random;
    for (int n = 0; n < 50; n++) {
        const Case<T, 4> c = make_case<T, 4>(random);
        NMPC::Vector<T, 4> reference = NMPC::Vector<T, 4>::Zero();
        Eigen::Array<bool, 4, 1> held_reference;
        Eigen::LLT<NMPC::Matrix<T, 4, 4>> factor;
        ASSERT_GT((NMPC::solve_box_qp<T, 4>(c.H, c.g, c.lower, c.upper, reference, held_reference, factor, 64)), 0);
        for (const double scale : {1e-6, 1e6}) {
            const NMPC::Matrix<T, 4, 4> H(c.H * T(scale));
            const NMPC::Vector<T, 4> g(c.g * T(scale));
            NMPC::Vector<T, 4> d = NMPC::Vector<T, 4>::Zero();
            Eigen::Array<bool, 4, 1> held;
            ASSERT_GT((NMPC::solve_box_qp<T, 4>(H, g, c.lower, c.upper, d, held, factor, 64)), 0);
            const double tolerance = std::is_same<T, float>::value ? 1e-4 : 1e-10;
            for (int i = 0; i < 4; i++) {
                EXPECT_NEAR(double(d(i)), double(reference(i)), tolerance) << "case " << n << ", scale " << scale;
                EXPECT_EQ(held(i), held_reference(i)) << "case " << n << ", scale " << scale;
            }
        }
    }
}

TYPED_TEST(BoxQPTest, WarmStartDoesNotChangeTheAnswer) {
    /* From the answer, from zero and from a corner, the same answer to rounding. */
    typedef TypeParam T;
    Sequence random;
    for (int n = 0; n < 50; n++) {
        const Case<T, 5> c = make_case<T, 5>(random);
        Eigen::Array<bool, 5, 1> held;
        Eigen::LLT<NMPC::Matrix<T, 5, 5>> factor;
        NMPC::Vector<T, 5> from_zero = NMPC::Vector<T, 5>::Zero();
        ASSERT_GT((NMPC::solve_box_qp<T, 5>(c.H, c.g, c.lower, c.upper, from_zero, held, factor, 64)), 0);
        NMPC::Vector<T, 5> from_corner = c.upper;
        ASSERT_GT((NMPC::solve_box_qp<T, 5>(c.H, c.g, c.lower, c.upper, from_corner, held, factor, 64)), 0);
        NMPC::Vector<T, 5> again = from_zero;
        ASSERT_GT((NMPC::solve_box_qp<T, 5>(c.H, c.g, c.lower, c.upper, again, held, factor, 64)), 0);
        const double tolerance = std::is_same<T, float>::value ? 1e-5 : 1e-13;
        for (int i = 0; i < 5; i++) {
            EXPECT_NEAR(double(from_corner(i)), double(from_zero(i)), tolerance) << "case " << n;
            EXPECT_NEAR(double(again(i)), double(from_zero(i)), tolerance) << "case " << n;
        }
    }
}

TYPED_TEST(BoxQPTest, AWholeStepFromFarAwayReachesItsTarget) {
    /*
    The whole Newton step sets the free variables to the target. Adding the step to a start far from it
    lost it to cancellation: in single precision, from 1e8 the answer to H = 1, g = -1 came out 0, with a
    gradient of -1 left, and the method stopped there. Double precision loses it from 1e17.
    */
    typedef TypeParam T;
    for (const double start : {1e8, -1e8, 1e12, 1e17}) {
        NMPC::Matrix<T, 1, 1> H;
        H << T(1);
        NMPC::Vector<T, 1> g, lower, upper, d;
        g << T(-1);
        lower << T(-1e30);
        upper << T(1e30);
        d << T(start);
        Eigen::Array<bool, 1, 1> held;
        Eigen::LLT<NMPC::Matrix<T, 1, 1>> factor;
        bool finished = false;
        ASSERT_GT((NMPC::solve_box_qp<T, 1>(H, g, lower, upper, d, held, factor, 32,
            T(16) * std::numeric_limits<T>::epsilon(), &finished)), 0) << start;
        EXPECT_EQ(d(0), T(1)) << start;
        EXPECT_TRUE(finished) << start;
    }
    /* Two variables, one ending on a bound, from far away. */
    NMPC::Matrix<T, 2, 2> H;
    H << T(2), T(0.5), T(0.5), T(1);
    const NMPC::Vector<T, 2> g(T(-1), T(4)), lower(T(-1e9), T(-1)), upper(T(1e9), T(1));
    NMPC::Vector<T, 2> d(T(1e8), T(0));
    Eigen::Array<bool, 2, 1> held;
    Eigen::LLT<NMPC::Matrix<T, 2, 2>> factor;
    ASSERT_GT((NMPC::solve_box_qp<T, 2>(H, g, lower, upper, d, held, factor, 32)), 0);
    EXPECT_EQ(d(1), T(-1));
    EXPECT_TRUE(held(1));
    EXPECT_NEAR(double(d(0)), (1.0 + 0.5) / 2.0, 1e-6);
}

TYPED_TEST(BoxQPTest, AFarBoundCrossedByTheTargetIsHeld) {
    /*
    A target past a bound always stops the step on that bound and holds it, even when the fraction of the
    step to the bound rounds to one. From -1e8 to the bound 0 with a target of 1, (0 - -1e8) / (1 - -1e8)
    is 1 in single precision: the step was taken as whole, d came out 1, outside the bounds, and the
    variable was left free. Double precision does the same from -1e17.
    */
    typedef TypeParam T;
    for (const double start : {-1e8, -1e12, -1e17, -1e30}) {
        NMPC::Matrix<T, 1, 1> H;
        H << T(1);
        NMPC::Vector<T, 1> g, lower, upper, d;
        g << T(-1);
        lower << T(start);
        upper << T(0);
        d << T(start);
        Eigen::Array<bool, 1, 1> held;
        Eigen::LLT<NMPC::Matrix<T, 1, 1>> factor;
        bool finished = false;
        ASSERT_GT((NMPC::solve_box_qp<T, 1>(H, g, lower, upper, d, held, factor, 32,
            T(16) * std::numeric_limits<T>::epsilon(), &finished)), 0) << start;
        EXPECT_EQ(d(0), T(0)) << start;
        EXPECT_TRUE(held(0)) << start;
        EXPECT_TRUE(finished) << start;
    }
    /* The same at the lower bound, and with a second variable that ends free. */
    for (const double start : {1e8, 1e17}) {
        NMPC::Matrix<T, 2, 2> H;
        H << T(1), T(0), T(0), T(2);
        const NMPC::Vector<T, 2> g(T(1), T(-1)), lower(T(0), T(-10)), upper(T(start), T(10));
        NMPC::Vector<T, 2> d(T(start), T(0));
        Eigen::Array<bool, 2, 1> held;
        Eigen::LLT<NMPC::Matrix<T, 2, 2>> factor;
        bool finished = false;
        ASSERT_GT((NMPC::solve_box_qp<T, 2>(H, g, lower, upper, d, held, factor, 32,
            T(16) * std::numeric_limits<T>::epsilon(), &finished)), 0) << start;
        EXPECT_EQ(d(0), T(0)) << start;
        EXPECT_TRUE(held(0)) << start;
        EXPECT_FALSE(held(1)) << start;
        EXPECT_NEAR(double(d(1)), 0.5, 1e-6) << start;
        EXPECT_TRUE(finished) << start;
    }
}

TYPED_TEST(BoxQPTest, ReportsWhetherItFinished) {
    typedef TypeParam T;
    Sequence random;
    const Case<T, 5> c = make_case<T, 5>(random);
    Eigen::Array<bool, 5, 1> held;
    Eigen::LLT<NMPC::Matrix<T, 5, 5>> factor;
    const T tolerance = T(16) * std::numeric_limits<T>::epsilon();
    bool finished = true;
    NMPC::Vector<T, 5> d = c.upper;
    EXPECT_EQ((NMPC::solve_box_qp<T, 5>(c.H, c.g, c.lower, c.upper, d, held, factor, 0, tolerance, &finished)), 0);
    EXPECT_FALSE(finished);
    d = c.upper;
    ASSERT_GT((NMPC::solve_box_qp<T, 5>(c.H, c.g, c.lower, c.upper, d, held, factor, 64, tolerance, &finished)), 0);
    EXPECT_TRUE(finished);
}
