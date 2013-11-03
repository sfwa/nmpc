/*
Copyright (C) 2013 Daniel Dyer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <cmath>
#include <stdint.h>
#include <acado_toolkit.hpp>
#include <include/acado_gnuplot/gnuplot_window.hpp>

#include "nmpc/state.h"
#include "nmpc/dynamics.h"
#include "nmpc/integrator.h"
#include "nmpc/types.h"
#include "debug.h"

#include <iostream>

#define SIM_TIMESTEP	0.01
#define SIM_LENGTH		25
#define HORIZON_LENGTH	10

void dynamics(double *x, double *f, void *user_data) {
	int i;

	/* Create new state vector using input array. */
	State s;
	s.position() << x[0], x[1], x[2];
	s.velocity() << x[3], x[4], x[5];
	s.acceleration() << 0, 0, 0;
	s.attitude() << x[9], x[10], x[11], x[12];
	Quaternionr q_1(s.attitude());
	q_1.normalize();
	s.attitude() << q_1.vec(), q_1.w();
	s.angular_velocity() << x[13], x[14], x[15];
	s.angular_acceleration() << 0, 0, 0;
	s.wind_velocity() << 0, 0, 0;

	for(i=0; i<19; i++) { AssertNotNaN(x[i]); }

	ControlVector c(4);
	c << x[19], x[20], x[21], 0;

	for(i=0; i<4; i++) { AssertNotNaN(c[i]); }

	/* Evaluate the dynamics model. */
	FixedWingFlightDynamicsModel *d =
		static_cast<FixedWingFlightDynamicsModel *>(user_data);
	AccelerationVector a = d->evaluate(s, c);

	for(i=0; i<6; i++) {
		AssertNotNaN(a[i]);
		if(abs(a[i]) > 80) {
			std::cout << "control:\n" << c << std::endl << std::endl;
			std::cout << "state:\n" << s << std::endl << std::endl;
			std::cout << "accel:\n" << a << std::endl << std::endl;
		}
	}

	/* Evaluate the process model using state and dynamics outputs. */
	s.acceleration() << a[0], a[1], a[2];
	s.angular_acceleration() << a[3], a[4], a[5];

	IntegratorRK4 integrator;
	State dot = integrator.integrate(s, SIM_TIMESTEP);
	Quaternionr q_2(dot.attitude());
	q_2.normalize();
	dot.attitude() << q_2.vec(), q_2.w();

	for(i=0; i<19; i++) { AssertNotNaN(dot[i]); }

	/* Copy results to output vector. */
	for(i=0; i<3; ++i) f[i] = dot.position()[i];
	for(i=0; i<3; ++i) f[i+3] = dot.velocity()[i];
	for(i=0; i<3; ++i) f[i+6] = dot.acceleration()[i];
	for(i=0; i<4; ++i) f[i+9] = dot.attitude()[i];
	for(i=0; i<3; ++i) f[i+13] = dot.angular_velocity()[i];
	for(i=0; i<3; ++i) f[i+16] = dot.angular_acceleration()[i];
}

int main()
{
    USING_NAMESPACE_ACADO

	int i;

	/* State vector. */
	DifferentialState state_vector(19);

	/* Control vector. */
	Control motor;
	Control elevon(2);

	/* Dynamics model. */
	FixedWingFlightDynamicsModel dynamics_model;

	/* Dynamics model parameters. */
	dynamics_model.set_mass(3.8);
	dynamics_model.set_inertia_tensor((Matrix3x3r() <<
		2.59e-1, 0, -0.334e-1,
		0, 1.47e-1, 0,
		-0.334e-1, 0, 4.05e-1).finished());
	dynamics_model.set_prop_coeffs(0.025, 0.00250);
	dynamics_model.set_lift_coeffs((Vector5r() <<
		-3.7, -5.4, 1.3, 1.7, 0.05).finished());
	dynamics_model.set_drag_coeffs((Vector5r() <<
		0.0, 0.0, 0.2, 0.0, 0.05).finished());
	dynamics_model.set_side_coeffs((Vector4r() <<
		0, 2.35e-01, -1.87e-03, 4.53e-04).finished(),
		(ControlVector(4) << 0.0, 1.1e-02, -1.1e-02, 0.0).finished());
	dynamics_model.set_pitch_moment_coeffs(
		(Vector2r() << -0.01, -0.0018).finished(),
		(ControlVector(4) << 0.0, -0.001, -0.001, 0.0).finished());
	dynamics_model.set_roll_moment_coeffs(
		(Vector1r() << -0.002).finished(),
		(ControlVector(4) << 0.0, -0.003, 0.003, 0.0).finished());
	dynamics_model.set_yaw_moment_coeffs(
		(Vector2r() << 0, -0.005).finished(),
		(ControlVector(4) << 0.0, 0.0, 0.0, 0.0).finished());
	dynamics_model.set_motor_index(0);

	/* Combine the above vectors into one to pass to the dynamics function. */
	IntermediateState is(22);

	for(i=0; i<19; ++i) is(i) = state_vector(i);
	is(19) = motor;
	for(i=0; i<2; ++i) is(i+20) = elevon(i);

	/* Define the differential equation. */
	CFunction flight_model(19, dynamics);
	flight_model.setUserData(&dynamics_model);
	DiscretizedDifferentialEquation f(SIM_TIMESTEP);

	/* Define the optimal control problem. */
	Function h;

	h << state_vector(0);
	h << state_vector(1);
	h << state_vector(2);
	h << state_vector(3);
	h << state_vector(4);
	h << state_vector(5);
	h << state_vector(9);
	h << state_vector(10);
	h << state_vector(11);
	h << state_vector(12);
	h << state_vector(13);
	h << state_vector(14);
	h << state_vector(15);

	/* Least-squares weighting matrix. */
	Matrix Q = zeros(13, 13);
	Q(0, 0) = 0.1;
	Q(1, 1) = 0.1;
	Q(2, 2) = 0.1;
	Q(3, 3) = 1.0;
	Q(4, 4) = 1.0;
	Q(5, 5) = 1.0;
	Q(6, 6) = 5.0;
	Q(7, 7) = 5.0;
	Q(8, 8) = 5.0;
	Q(9, 9) = 5.0;
	Q(10, 10) = 1.0;
	Q(11, 11) = 1.0;
	Q(12, 12) = 1.0;

	OCP ocp(0.0, HORIZON_LENGTH, 50);

	Vector r(13);
	r.setZero();

	ocp.minimizeLSQ(Q, h, r);

	ocp.subjectTo(f << flight_model(is));

	/* Flight envelope constraints. */

	/* Acceleration. */
	ocp.subjectTo(-80 <= state_vector(6) <= 80);
	ocp.subjectTo(-10 <= state_vector(7) <= 10);
	ocp.subjectTo(-80 <= state_vector(8) <= 80);

	/* Angular velocity. */
	ocp.subjectTo(-M_PI <= state_vector(13) <= M_PI);
	ocp.subjectTo(-M_PI <= state_vector(14) <= M_PI);
	ocp.subjectTo(-M_PI_4 <= state_vector(15) <= M_PI_4);

	/* Control constraints. */
	ocp.subjectTo(0.0 <= motor <= 30000.0);
	ocp.subjectTo(-1.0 <= elevon <= 1.0);

	/* Set up the simulation process. */
	OutputFcn identity;
	DynamicSystem dynamicSystem(f, identity);
	Process process(dynamicSystem, INT_RK78);

	/* Set up the realtime algorithm. */
	RealTimeAlgorithm alg(ocp, SIM_TIMESTEP);
	alg.set(INTEGRATOR_TYPE, INT_RK78);
	//alg.set(INTEGRATOR_PRINTLEVEL, HIGH);

	/* Set up controller and load reference trajectory. */
	VariablesGrid ref_data("data/reference.txt");
	StaticReferenceTrajectory reference(ref_data);
	Controller controller(alg, reference);

	/* Initialise and run the simulation. */
	SimulationEnvironment sim(0.0, SIM_LENGTH, process, controller);

	Vector r0(13);
	r0 = ref_data.getFirstVector();
	Vector x0(19);
	x0.setZero();
	x0(0) = r0(0);
	x0(1) = r0(1);
	x0(2) = r0(2);
	x0(3) = r0(3);
	x0(4) = r0(4);
	x0(5) = r0(5);
	x0(9) = r0(6);
	x0(10) = r0(7);
	x0(11) = r0(8);
	x0(12) = r0(9);
	x0(13) = r0(10);
	x0(14) = r0(11);
	x0(15) = r0(12);
	sim.init(x0);
	sim.run();

	/* Plot the results. */
	VariablesGrid diffStates;
	sim.getProcessDifferentialStates(diffStates);

	VariablesGrid feedbackControl;
	sim.getFeedbackControl(feedbackControl);

	OutputFcn plot_fn;
	plot_fn << state_vector(2);
	plot_fn << state_vector(3);
	plot_fn << state_vector(4);
	plot_fn << state_vector(5);
	plot_fn << state_vector(10);
	plot_fn << state_vector(11);
	plot_fn << state_vector(12);
	VariablesGrid plot_out;
	plot_fn.evaluate(&diffStates, NULL, NULL, NULL, NULL, &plot_out);

	GnuplotWindow window;
	window.addSubplot(plot_out(0), "Altitude [m]");
	window.addSubplot(plot_out(1), "Velocity X [m/s]");
	window.addSubplot(plot_out(2), "Velocity Y [m/s]");
	window.addSubplot(plot_out(3), "Velocity Z [m/s]");
	window.addSubplot(plot_out(4), "Angular Rate P [rad/s]");
	window.addSubplot(plot_out(5), "Angular Rate Q [rad/s]");
	window.addSubplot(plot_out(6), "Angular Rate R [rad/s]");
	window.addSubplot(feedbackControl(0), "Throttle");
	window.addSubplot(feedbackControl(1), "Elevon 1 Deflection [rad]");
	window.addSubplot(feedbackControl(2), "Elevon 2 Deflection [rad]");
	window.plot();

	OutputFcn xplane_fn;
	xplane_fn << state_vector(0);
	xplane_fn << state_vector(1);
	xplane_fn << state_vector(2);
	xplane_fn << state_vector(9);
	xplane_fn << state_vector(10);
	xplane_fn << state_vector(11);
	xplane_fn << state_vector(12);
	VariablesGrid xplane_out;
	xplane_fn.evaluate(&diffStates, NULL, NULL, NULL, NULL, &xplane_out);
	Grid xplane_grid(0.0, SIM_LENGTH-SIM_TIMESTEP, SIM_LENGTH / SIM_TIMESTEP);
	xplane_out.coarsenGrid(xplane_grid);
	xplane_out.print("", "", "", DEFAULT_WIDTH, DEFAULT_PRECISION, "\t", "\n");

    return 0;
}



