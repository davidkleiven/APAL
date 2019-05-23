cdef extern from "test_build_matrices.hpp":
    object biharmonic_matrix()

def pytest_biharmonic_matrix():
    return biharmonic_matrix()