cdef extern from "test_build_matrices.hpp":
    object biharmonic_matrix() except+
    object laplacian_matrix3D() except+
    object get_small_biharmonic() except+

def pytest_biharmonic_matrix():
    return biharmonic_matrix()

def pytest_laplacian_matrix3D():
    return laplacian_matrix3D()

def pytest_get_small_biharmonic():
    return get_small_biharmonic()