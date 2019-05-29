from libcpp.vector cimport vector


cdef extern from "test_fftw.hpp":
    int fft1D(vector[double] &vec, int direction)
    int fft2D(vector[double] &vec, int direction)
    int fft3D(vector[double] &vec, int direction)


def pyfft1D(array, direction):
    cdef vector[double] vec
    for i in range(len(array)):
        vec.push_back(array[i])

    ret_code = fft1D(vec, direction)

    if ret_code != 0:
        return None

    for i in range(len(array)):
        array[i] = vec[i]
    return array


def pyfft2D(array, direction):
    cdef vector[double] vec
    shape = array.shape
    array = array.ravel()

    for i in range(len(array)):
        vec.push_back(array[i])
    
    ret_code = fft2D(vec, direction)

    if ret_code != 0:
        return None

    for i in range(len(array)):
        array[i] = vec[i]
    array = array.reshape(shape)
    return array

def pyfft3D(array, direction):
    cdef vector[double] vec
    shape = array.shape
    array = array.ravel()

    for i in range(len(array)):
        vec.push_back(array[i])
    
    ret_code = fft3D(vec, direction)

    if ret_code != 0:
        return None

    for i in range(len(array)):
        array[i] = vec[i]
    array = array.reshape(shape)
    return array