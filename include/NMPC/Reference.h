/*
Part of NMPC. Distributed under the MIT licence: see LICENSE in the root of
the repository.
*/

#ifndef NMPC_REFERENCE_H
#define NMPC_REFERENCE_H

#include "NMPC/Types.h"
#include <algorithm>

namespace NMPC {

/*
A reference as a table: rows of a fixed number of columns at a fixed spacing
in time, with no heap. The table may hold only part of a flight's rows (a
window), in which case `first` and `held` say which; a lookup outside the
window holds the window's nearest end and is counted in `misses`, so that an
application can tell a window that was moved too late.

The table knows nothing of what the columns mean. span() finds the two rows
either side of a time and the fraction between them; the vehicle's path type
(for the multirotor, SFWA::Guidance::MultirotorPath) turns them into a state
and a command, interpolating each part in its own way. Before the first row
and after the last, the ends hold.
*/
template <class T, int Columns>
struct PathTable {
    using Time = T;
    using Row = Vector<T, Columns>;
    enum { Width = Columns };

    const T* rows = nullptr;    /* the rows held, Columns numbers each */
    int count = 0;              /* rows of the whole table */
    int first = 0;              /* the first row held, of the whole table */
    int held = 0;               /* rows held */
    T dt = T(0.002);            /* spacing of the rows, s */
    mutable int misses = 0;     /* lookups that the rows held did not cover */

    struct Span {
        Eigen::Map<const Row> r0, r1;
        T a;                    /* 0 at r0, 1 at r1 */
    };

    /* The whole table at once. */
    void hold(const T* all, int rows_in_all, T spacing) {
        rows = all;
        count = held = rows_in_all;
        first = 0;
        dt = spacing;
    }

    T duration() const { return T(count - 1) * dt; }

    /*
    A lookup outside the window is moved to the nearest end of the window, the time and not only the rows,
    so that the reference holds there as it does at the ends of the table instead of jumping back and
    ramping through the last rows again.
    */
    Span span(T t) const {
        NMPC_ASSERT(count >= 2 && held >= 2 && first >= 0 && first + held <= count);
        T f = std::min(std::max(t / dt, T(0)), T(count - 1));
        const T start = T(first), end = T(first + held - 1);
        if (f < start || f > end) {
            misses++;
            f = std::min(std::max(f, start), end);
        }
        const int i = std::min(int(f), first + held - 2);
        const T a = f - T(i);
        const int local = i - first;
        return Span{Eigen::Map<const Row>(rows + Columns * local), Eigen::Map<const Row>(rows + Columns * (local + 1)), a};
    }
};

/*
The reference's clock, summed every tick. A millisecond is a few hundred units
of the last place of a single-precision clock of some seconds, so a plain sum
loses up to half a microsecond a tick, always in the same direction: 2 to 5 ms
over a flight of 20 to 30 seconds. The sum carries what it rounds off (Kahan).
In double precision the plain sum is exact enough, and is what the review
prototype did, so compensation is on by default for types narrower than double.
*/
template <class T>
struct Clock {
    T value = T(0);
    T carry = T(0);             /* what the last sums rounded off */
    bool compensated = sizeof(T) < sizeof(double);

    void reset(T start = T(0)) {
        value = start;
        carry = T(0);
    }

    void advance(T step) {
        if (compensated) {
            const T y = step - carry, next = value + y;
            carry = (next - value) - y;
            value = next;
        } else {
            value += step;
        }
    }
};

}

#endif
