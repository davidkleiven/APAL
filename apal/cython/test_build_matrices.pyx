cdef extern from "test_build_matrices.hpp":
    object biharmonic_matrix()
    object laplacian_matrix3D()
    object test_add_matrices_same_elements()
    object test_add_matrices_mixed_elements()

def pytest_biharmonic_matrix():
    return biharmonic_matrix()

def pytest_laplacian_matrix3D():
    return laplacian_matrix3D()

def pytest_add_same_element():
    return test_add_matrices_same_elements()

def pytest_add_mixed_element():
    return test_add_matrices_mixed_elements()