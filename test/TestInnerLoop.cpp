#include "Typed.h"
#include "NMPC/InnerLoop.h"
#include "NMPC/Reference.h"
#include <vector>

/*
The reference table and clock, and the inner loop's building blocks.
*/

template <typename T>
class InnerLoopTest : public ::testing::Test {};

typedef ::testing::Types<float, double> Precisions;
TYPED_TEST_SUITE(InnerLoopTest, Precisions);

TYPED_TEST(InnerLoopTest, TableInterpolatesAndHoldsItsEnds) {
    typedef TypeParam T;
    /* Three columns: t, 2t and 1 - t, at 0.5 s. */
    std::vector<T> rows;
    for (int i = 0; i < 5; i++) {
        const T t = T(0.5) * T(i);
        rows.push_back(t);
        rows.push_back(T(2) * t);
        rows.push_back(T(1) - t);
    }
    NMPC::PathTable<T, 3> table;
    table.hold(rows.data(), 5, T(0.5));
    EXPECT_NEAR(table.duration(), T(2), T(1e-6));
    const auto s = table.span(T(1.2));
    const NMPC::Vector<T, 3> r(s.r0 + s.a * (s.r1 - s.r0));
    EXPECT_NEAR(r(0), T(1.2), T(1e-5));
    EXPECT_NEAR(r(1), T(2.4), T(1e-5));
    EXPECT_NEAR(r(2), T(-0.2), T(1e-5));
    /* Before the first row and after the last the ends hold. */
    const auto before = table.span(T(-1)), after = table.span(T(9));
    EXPECT_EQ(before.r0(0) + before.a * (before.r1(0) - before.r0(0)), T(0));
    EXPECT_EQ(after.r0(0) + after.a * (after.r1(0) - after.r0(0)), T(2));
    EXPECT_EQ(table.misses, 0);
}

TYPED_TEST(InnerLoopTest, WindowCountsLookupsOutsideIt) {
    typedef TypeParam T;
    std::vector<T> all;
    for (int i = 0; i < 10; i++) {
        all.push_back(T(i));
    }
    /* Rows 4 to 6 of a table of ten held. */
    NMPC::PathTable<T, 1> table;
    table.rows = all.data() + 4;
    table.count = 10;
    table.first = 4;
    table.held = 3;
    table.dt = T(1);
    const auto inside = table.span(T(4.5));
    EXPECT_EQ(inside.r0(0), T(4));
    EXPECT_EQ(table.misses, 0);
    const auto outside = table.span(T(8.5));
    EXPECT_EQ(table.misses, 1);
    /* The nearest rows held stand in. */
    EXPECT_EQ(outside.r0(0), T(5));
    EXPECT_EQ(outside.r1(0), T(6));
}

TYPED_TEST(InnerLoopTest, ClockKeepsTimeOverAFlight) {
    typedef TypeParam T;
    /* 30 s at 1 kHz with a rate of 0.9: the clock must end at 27 s to well under a millisecond. */
    NMPC::Clock<T> clock;
    for (long i = 0; i < 30000; i++) {
        clock.advance(T(0.001) * T(0.9));
    }
    EXPECT_NEAR(double(clock.value), 27.0, 1e-5);
    if (sizeof(T) < sizeof(double)) {
        /* Without compensation single precision drifts by milliseconds. */
        NMPC::Clock<T> plain;
        plain.compensated = false;
        for (long i = 0; i < 30000; i++) {
            plain.advance(T(0.001) * T(0.9));
        }
        EXPECT_GT(std::abs(double(plain.value) - 27.0), 1e-4);
    }
}

TYPED_TEST(InnerLoopTest, ButterworthPassesSteadyAndStopsFast) {
    typedef TypeParam T;
    NMPC::Butterworth2<T> f;
    f.design(T(30), T(1000));
    T y = T(0);
    for (int i = 0; i < 2000; i++) {
        y = f(T(1));
    }
    EXPECT_NEAR(y, T(1), T(1e-4));
    /* A tone at 300 Hz, ten times the cutoff, comes through at about a hundredth. */
    f.reset();
    T peak = T(0);
    for (int i = 0; i < 2000; i++) {
        const T out = f(T(std::sin(2 * 3.14159265358979 * 300 * 0.001 * i)));
        if (i > 1000) {
            peak = std::max(peak, std::abs(out));
        }
    }
    EXPECT_LT(peak, T(0.02));
}

TYPED_TEST(InnerLoopTest, IncrementalAllocatesWithinItsLimits) {
    typedef TypeParam T;
    NMPC::Incremental<T, 2, 2> c;
    c.design(T(30), T(1000));
    c.gain << T(2), T(-1);
    c.allocation << T(1), T(1), T(1), T(-1);
    c.limit << T(10), T(0.5);
    NMPC::Vector<T, 2> held = NMPC::Vector<T, 2>::Zero();
    for (int i = 0; i < 2000; i++) {
        held = c.update(NMPC::Vector<T, 2>(T(1), T(1)));
    }
    /* Wanted 2 and -1 once the filters settle: increments 1 and 3, the second held to its limit. */
    EXPECT_NEAR(held(0), T(1), T(1e-3));
    EXPECT_EQ(held(1), T(0.5));
    c.reset();
    EXPECT_EQ(c.held(0), T(0));
}

TYPED_TEST(InnerLoopTest, ClampAndBound) {
    typedef TypeParam T;
    const NMPC::Vector<T, 3> lower(T(0), T(0), T(0)), upper(T(1), T(1), T(1));
    const NMPC::Vector<T, 3> u(T(-1), T(0.5), T(2));
    const NMPC::Vector<T, 3> c = NMPC::clamp(u, lower, upper);
    EXPECT_EQ(c(0), T(0));
    EXPECT_EQ(c(1), T(0.5));
    EXPECT_EQ(c(2), T(1));
    EXPECT_TRUE(NMPC::on_bound(u, lower, upper));
    const NMPC::Vector<T, 3> middle = NMPC::Vector<T, 3>::Constant(T(0.5));
    EXPECT_FALSE(NMPC::on_bound(middle, lower, upper));
}

TYPED_TEST(InnerLoopTest, CascadeTurnsBackTowardThePlan) {
    typedef TypeParam T;
    NMPC::AttitudeCascade<T> c;
    c.design(T(0.001));
    const NMPC::Vector<T, 3> inertia(T(0.02), T(0.02), T(0.04)), zero = NMPC::Vector<T, 3>::Zero();
    const NMPC::Quaternion<T> level = NMPC::Quaternion<T>::Identity();

    /* On the plan: no moment, and the integral stays empty. */
    EXPECT_TRUE(c.moment(level, zero, level, zero, inertia).isZero());
    c.integrate(false, inertia);
    EXPECT_TRUE(c.integral.isZero());

    /* Rolled right of the plan: a moment to roll left, and only about x. */
    const T angle = T(0.1);
    const NMPC::Quaternion<T> rolled(Eigen::AngleAxis<T>(angle, NMPC::Vector<T, 3>::UnitX()));
    c.reset();
    const NMPC::Vector<T, 3> m = c.moment(rolled, zero, level, zero, inertia);
    EXPECT_LT(m.x(), T(0));
    EXPECT_NEAR(m.y(), T(0), T(1e-6));
    EXPECT_NEAR(m.z(), T(0), T(1e-6));
    /* The rate error is the attitude loop's rate command. */
    const T error = -c.k_attitude.x() * angle;
    EXPECT_NEAR(c.error.x(), error, T(1e-5));

    /* The integral grows toward its limit, and stops while the command is clipped. */
    for (int i = 0; i < 20000; i++) {
        c.moment(rolled, zero, level, zero, inertia);
        c.integrate(false, inertia);
    }
    EXPECT_EQ(c.integral.x(), -c.integral_limit.x());
    c.reset();
    c.moment(rolled, zero, level, zero, inertia);
    c.integrate(true, inertia);
    EXPECT_TRUE(c.integral.isZero());
}
