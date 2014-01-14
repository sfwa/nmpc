import os
from ctypes import *


# Taken from c/cnmpc.h
NMPC_PRECISION_FLOAT = 0
NMPC_PRECISION_DOUBLE = 1

state = None

# Internal globals, set during init
_cnmpc = None
_REAL_T = None
_CONTROL_DIM = None


class _State(Structure):
    def __repr__(self):
        fields = {
            "position": tuple(self.position),
            "velocity": tuple(self.velocity),
            "attitude": tuple(self.attitude),
            "angular_velocity": tuple(self.angular_velocity),
            "wind_velocity": tuple(self.wind_velocity)
        }
        return str(fields)


# Public interface
def integrate(dt, control=None):
    global _cnmpc, state

    if not _cnmpc:
        raise RuntimeError("Please call nmpc.init()")

    if control is None:
        control = (0.0, ) * _CONTROL_DIM
    elif len(control) != _CONTROL_DIM:
        raise ValueError("Control vector must contain %d elements" %
                         _CONTROL_DIM)

    _cnmpc.nmpc_fixedwingdynamics_set_state(state)
    _cnmpc.nmpc_fixedwingdynamics_integrate(
        dt, (_REAL_T * _CONTROL_DIM)(*control))
    _cnmpc.nmpc_fixedwingdynamics_get_state(state)


def configure_airframe(mass=None, inertia_tensor=None, prop_coeffs=None,
        drag_coeffs=None, lift_coeffs=None, side_coeffs=None,
        pitch_moment_coeffs=None, roll_moment_coeffs=None,
        yaw_moment_coeffs=None):
    _cnmpc.nmpc_fixedwingdynamics_set_mass(mass)
    _cnmpc.nmpc_fixedwingdynamics_set_inertia_tensor(
        (_REAL_T * 9)(*inertia_tensor))
    _cnmpc.nmpc_fixedwingdynamics_set_prop_coeffs(prop_coeffs[0],
        prop_coeffs[1])
    _cnmpc.nmpc_fixedwingdynamics_set_lift_coeffs((_REAL_T * 5)(*lift_coeffs))
    _cnmpc.nmpc_fixedwingdynamics_set_drag_coeffs((_REAL_T * 5)(*drag_coeffs))
    _cnmpc.nmpc_fixedwingdynamics_set_side_coeffs(
        (_REAL_T * 4)(*side_coeffs[0:-_CONTROL_DIM]),
        (_REAL_T * _CONTROL_DIM)(*side_coeffs[-_CONTROL_DIM:]))
    _cnmpc.nmpc_fixedwingdynamics_set_pitch_moment_coeffs(
        (_REAL_T * 2)(*pitch_moment_coeffs[0:-_CONTROL_DIM]),
        (_REAL_T * _CONTROL_DIM)(*pitch_moment_coeffs[-_CONTROL_DIM:]))
    _cnmpc.nmpc_fixedwingdynamics_set_roll_moment_coeffs(
        (_REAL_T * 1)(*roll_moment_coeffs[0:-_CONTROL_DIM]),
        (_REAL_T * _CONTROL_DIM)(*roll_moment_coeffs[-_CONTROL_DIM:]))
    _cnmpc.nmpc_fixedwingdynamics_set_yaw_moment_coeffs(
        (_REAL_T * 2)(*yaw_moment_coeffs[0:-_CONTROL_DIM]),
        (_REAL_T * _CONTROL_DIM)(*yaw_moment_coeffs[-_CONTROL_DIM:]))


def init(implementation="c"):
    global _cnmpc, _REAL_T, _CONTROL_DIM, state

    # Load the requested library and determine configuration parameters
    if implementation == "c":
        lib = os.path.join(os.path.dirname(__file__), "c", "libcnmpc.dylib")
    else:
        raise NameError(
            "Unknown NMPC implementation: %s (options are 'c')" %
            implementation)

    _cnmpc = cdll.LoadLibrary(lib)

    _cnmpc.nmpc_init.argtypes = []
    _cnmpc.nmpc_init.restype = None

    _cnmpc.nmpc_config_get_precision.argtypes = []
    _cnmpc.nmpc_config_get_precision.restype = c_long

    _cnmpc.nmpc_config_get_state_dim.argtypes = []
    _cnmpc.nmpc_config_get_state_dim.restype = c_long

    _cnmpc.nmpc_config_get_control_dim.argtypes = []
    _cnmpc.nmpc_config_get_control_dim.restype = c_long

    _PRECISION = _cnmpc.nmpc_config_get_precision()
    _REAL_T = c_double if _PRECISION == NMPC_PRECISION_DOUBLE else c_float
    _CONTROL_DIM = _cnmpc.nmpc_config_get_control_dim()
    _STATE_DIM = _cnmpc.nmpc_config_get_state_dim()

    _State._fields_ = [
        ("position", _REAL_T * 3),
        ("velocity", _REAL_T * 3),
        ("attitude", _REAL_T * 4),
        ("angular_velocity", _REAL_T * 3),
        ("wind_velocity", _REAL_T * 3)
    ]

    _cnmpc.nmpc_preparation_step.argtypes = []
    _cnmpc.nmpc_preparation_step.restype = None

    _cnmpc.nmpc_feedback_step.argtype = [POINTER(_State)]
    _cnmpc.nmpc_feedback_step.restype = None

    _cnmpc.nmpc_set_state_weights.argtype = [
        POINTER(_REAL_T * (_STATE_DIM-1))]
    _cnmpc.nmpc_set_state_weights.restype = None

    _cnmpc.nmpc_set_control_weights.argtype = [
        POINTER(_REAL_T * (_CONTROL_DIM))]
    _cnmpc.nmpc_set_control_weights.restype = None

    _cnmpc.nmpc_set_terminal_weights.argtype = [
        POINTER(_REAL_T * (_STATE_DIM-1))]
    _cnmpc.nmpc_set_terminal_weights.restype = None

    _cnmpc.nmpc_set_lower_control_bound.argtype = [
        POINTER(_REAL_T * (_CONTROL_DIM))]
    _cnmpc.nmpc_set_lower_control_bound.restype = None

    _cnmpc.nmpc_set_upper_control_bound.argtype = [
        POINTER(_REAL_T * (_CONTROL_DIM))]
    _cnmpc.nmpc_set_upper_control_bound.restype = None

    _cnmpc.nmpc_set_reference_point.argtype = [
        POINTER(_REAL_T * (_STATE_DIM+_CONTROL_DIM)),
        c_uint]
    _cnmpc.nmpc_set_reference_point.restype = None

    # Set up the function prototypes
    _cnmpc.nmpc_fixedwingdynamics_set_position.argtypes =
        [_REAL_T, _REAL_T, _REAL_T]
    _cnmpc.nmpc_fixedwingdynamics_set_position.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_velocity.argtypes =
        [_REAL_T, _REAL_T, _REAL_T]
    _cnmpc.nmpc_fixedwingdynamics_set_velocity.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_attitude.argtypes =
        [_REAL_T, _REAL_T, _REAL_T, _REAL_T]
    _cnmpc.nmpc_fixedwingdynamics_set_attitude.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_angular_velocity.argtypes =
        [_REAL_T, _REAL_T, _REAL_T]
    _cnmpc.nmpc_fixedwingdynamics_set_angular_velocity.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_wind_velocity.argtypes =
        [_REAL_T, _REAL_T, _REAL_T]
    _cnmpc.nmpc_fixedwingdynamics_set_wind_velocity.restype = None

    _cnmpc.nmpc_fixedwingdynamics_get_state.argtypes =
        [POINTER(_State)]
    _cnmpc.nmpc_fixedwingdynamics_get_state.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_state.argtypes =
        [POINTER(_State)]
    _cnmpc.nmpc_fixedwingdynamics_set_state.restype = None

    _cnmpc.nmpc_fixedwingdynamics_integrate.argtypes =
        [c_float, POINTER(_REAL_T * _CONTROL_DIM)]
    _cnmpc.nmpc_fixedwingdynamics_integrate.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_mass.argtypes = [_REAL_T]
    _cnmpc.nmpc_fixedwingdynamics_set_mass.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_inertia_tensor.argtypes = [
        POINTER(_REAL_T * 9)]
    _cnmpc.nmpc_fixedwingdynamics_set_inertia_tensor.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_prop_coeffs.argtypes = [
        _REAL_T, _REAL_T]
    _cnmpc.nmpc_fixedwingdynamics_set_prop_coeffs.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_drag_coeffs.argtypes = [
        POINTER(_REAL_T * 5)]
    _cnmpc.nmpc_fixedwingdynamics_set_drag_coeffs.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_lift_coeffs.argtypes = [
        POINTER(_REAL_T * 5)]
    _cnmpc.nmpc_fixedwingdynamics_set_lift_coeffs.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_side_coeffs.argtypes = [
        POINTER(_REAL_T * 4),
        POINTER(_REAL_T * _CONTROL_DIM)]
    _cnmpc.nmpc_fixedwingdynamics_set_side_coeffs.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_pitch_moment_coeffs.argtypes = [
        POINTER(_REAL_T * 2), POINTER(_REAL_T * _CONTROL_DIM)]
    _cnmpc.nmpc_fixedwingdynamics_set_pitch_moment_coeffs.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_roll_moment_coeffs.argtypes = [
        POINTER(_REAL_T * 1), POINTER(_REAL_T * _CONTROL_DIM)]
    _cnmpc.nmpc_fixedwingdynamics_set_roll_moment_coeffs.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_yaw_moment_coeffs.argtypes = [
        POINTER(_REAL_T * 2), POINTER(_REAL_T * _CONTROL_DIM)]
    _cnmpc.nmpc_fixedwingdynamics_set_yaw_moment_coeffs.restype = None

    # Set up the state
    state = _State()
    _cnmpc.nmpc_fixedwingdynamics_get_state(state)
