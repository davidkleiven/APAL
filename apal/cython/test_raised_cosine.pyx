cdef extern from "test_raised_cosine_filter.hpp":
    object test_omega_min()
    object test_omega_max()
    object test_eval_low_freq()
    object test_eval_at_cut()
    object test_eval_large_freq()

def pytest_omega_min():
    return test_omega_min()

def pytest_omega_max():
    return test_omega_max()

def pytest_eval_low_freq():
    return test_eval_low_freq()

def pytest_eval_at_cut():
    return test_eval_at_cut()

def pytest_eval_large_freq():
    return test_eval_large_freq()