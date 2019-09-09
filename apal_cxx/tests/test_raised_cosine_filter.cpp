#include "raised_cosine.hpp"
#include <cmath>

const double OMEGA_CUT = 1.5;
const double ROLL_OFF = 0.4;
const RaisedCosine cos_filter(OMEGA_CUT, ROLL_OFF);

PyObject* test_omega_min(){
    double expected = 0.6*1.5;
    if (abs(cos_filter.omega_min() - expected) < 1E-6){
        Py_RETURN_TRUE;
    }
    Py_RETURN_FALSE;
}

PyObject* test_omega_max(){
    double expected = 1.4*1.5;
    if (abs(cos_filter.omega_max() - expected) < 1E-6){
        Py_RETURN_TRUE;
    }
    Py_RETURN_FALSE;
}

PyObject* test_eval_low_freq(){
    double value = cos_filter.evaluate(0.1);
    if (abs(value - 1.0) < 1E-6){
        Py_RETURN_TRUE;
    }
    Py_RETURN_FALSE;
}

PyObject* test_eval_at_cut(){
    double value = cos_filter.evaluate(OMEGA_CUT);
    if (abs(value - 0.5) < 1E-6){
        Py_RETURN_TRUE;
    }
    Py_RETURN_FALSE;
}

PyObject* test_eval_large_freq(){
    double value = cos_filter.evaluate(10.0);
    if (abs(value) < 1E-6){
        Py_RETURN_TRUE;
    }
    Py_RETURN_FALSE;
}