from libcpp.vector cimport vector


cdef extern from "test_fftw.hpp":
    void fft1D(vector[double] &vec, int direction)
    void fft2D(vector[double] &vec, int direction)
    void fft3D(vector[double] &vec, int direction)


def pyfft1D(array, direction):
    cdef vector[double] vec
    for i in range(len(array)):
        vec.push_back(array[i])

    fft1D(vec, direction)

    for i in range(len(array)):
        array[i] = vec[i]
    return array


def pyfft2D(array, direction):
    cdef vector[double] vec
    shape = array.shape
    array = array.ravel()

    for i in range(len(array)):
        vec.push_back(array[i])
    
    fft2D(vec, direction)

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
    
    fft3D(vec, direction)

    for i in range(len(array)):
        array[i] = vec[i]
    array = array.reshape(shape)
    return array