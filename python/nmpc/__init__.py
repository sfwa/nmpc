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
_STATE_DIM = None

# Externally accessible globals
HORIZON_LENGTH = None
STEP_LENGTH = None

class _State(Structure):
    def __repr__(self):
        fields = {
            "position": tuple(self.position),
            "velocity": tuple(self.velocity),
            "attitude": tuple(self.attitude),
            "angular_velocity": tuple(self.angular_velocity)
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


def setup(state_weights, control_weights, terminal_weights,
        upper_control_bound, lower_control_bound):
    _cnmpc.nmpc_set_state_weights(
        (_REAL_T * (_STATE_DIM-1))(*state_weights))
    _cnmpc.nmpc_set_control_weights(
        (_REAL_T * (_CONTROL_DIM))(*control_weights))
    _cnmpc.nmpc_set_terminal_weights(
        (_REAL_T * (_STATE_DIM-1))(*terminal_weights))
    _cnmpc.nmpc_set_lower_control_bound(
        (_REAL_T * (_CONTROL_DIM))(*lower_control_bound))
    _cnmpc.nmpc_set_upper_control_bound(
        (_REAL_T * (_CONTROL_DIM))(*upper_control_bound))

def set_reference(point, index):
    _cnmpc.nmpc_set_reference_point(
        (_REAL_T * (_CONTROL_DIM+_STATE_DIM))(*point),
        index)

def initialise_horizon():
    _cnmpc.nmpc_init()

def prepare():
    _cnmpc.nmpc_preparation_step()

def solve(measurement):
    _cnmpc.nmpc_feedback_step((_REAL_T * (_STATE_DIM))(*measurement))

def get_controls():
    new_controls = (_REAL_T * _CONTROL_DIM)(*([0] * _CONTROL_DIM))
    _cnmpc.nmpc_get_controls(pointer(new_controls))

    return list(new_controls)

def set_wind_velocity(wind_velocity):
    _cnmpc.nmpc_set_wind_velocity(
        _REAL_T(wind_velocity[0]),
        _REAL_T(wind_velocity[1]),
        _REAL_T(wind_velocity[2]))

def update_horizon(new_reference):
    _cnmpc.nmpc_update_horizon(
        (_REAL_T * (_STATE_DIM+_CONTROL_DIM))(*new_reference))

def init(implementation="c"):
    global _cnmpc, _REAL_T, _STATE_DIM, _CONTROL_DIM, state
    global HORIZON_LENGTH, STEP_LENGTH

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

    _cnmpc.nmpc_config_get_horizon_length.argtypes = []
    _cnmpc.nmpc_config_get_horizon_length.restype = c_long

    _cnmpc.nmpc_config_get_step_length.argtypes = []
    _cnmpc.nmpc_config_get_step_length.restype = _REAL_T

    HORIZON_LENGTH = _cnmpc.nmpc_config_get_horizon_length()
    STEP_LENGTH = _cnmpc.nmpc_config_get_step_length()

    _State._fields_ = [
        ("position", _REAL_T * 3),
        ("velocity", _REAL_T * 3),
        ("attitude", _REAL_T * 4),
        ("angular_velocity", _REAL_T * 3)
    ]

    _cnmpc.nmpc_preparation_step.argtypes = []
    _cnmpc.nmpc_preparation_step.restype = None

    _cnmpc.nmpc_feedback_step.argtype = [POINTER(_State)]
    _cnmpc.nmpc_feedback_step.restype = None

    _cnmpc.nmpc_get_controls.argtype = [POINTER(_REAL_T * _CONTROL_DIM)]
    _cnmpc.nmpc_get_controls.restype = None

    _cnmpc.nmpc_update_horizon.argtype = [
        POINTER(_REAL_T * (_STATE_DIM+_CONTROL_DIM))]
    _cnmpc.nmpc_update_horizon.restype = None

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

    _cnmpc.nmpc_set_wind_velocity.argtype = [
        _REAL_T, _REAL_T, _REAL_T]
    _cnmpc.nmpc_set_wind_velocity.restype = None

    # Set up the function prototypes
    _cnmpc.nmpc_fixedwingdynamics_set_position.argtypes = [
        _REAL_T, _REAL_T, _REAL_T]
    _cnmpc.nmpc_fixedwingdynamics_set_position.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_velocity.argtypes = [
        _REAL_T, _REAL_T, _REAL_T]
    _cnmpc.nmpc_fixedwingdynamics_set_velocity.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_attitude.argtypes = [
        _REAL_T, _REAL_T, _REAL_T, _REAL_T]
    _cnmpc.nmpc_fixedwingdynamics_set_attitude.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_angular_velocity.argtypes = [
        _REAL_T, _REAL_T, _REAL_T]
    _cnmpc.nmpc_fixedwingdynamics_set_angular_velocity.restype = None

    _cnmpc.nmpc_fixedwingdynamics_get_state.argtypes = [
        POINTER(_State)]
    _cnmpc.nmpc_fixedwingdynamics_get_state.restype = None

    _cnmpc.nmpc_fixedwingdynamics_set_state.argtypes = [
        POINTER(_State)]
    _cnmpc.nmpc_fixedwingdynamics_set_state.restype = None

    _cnmpc.nmpc_fixedwingdynamics_integrate.argtypes = [
        c_float, POINTER(_REAL_T * _CONTROL_DIM)]
    _cnmpc.nmpc_fixedwingdynamics_integrate.restype = None

    # Set up the state
    state = _State()
    _cnmpc.nmpc_fixedwingdynamics_get_state(state)
