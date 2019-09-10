#include "test_vandeven.hpp"
#include "vandeven.hpp"
#include <iostream>

using namespace std;

PyObject* test_order1(){
    Vandeven vandeven(1);
    double dx = 1.0/9.0;
    for (unsigned int i=0;i<10;i++){
        double value = vandeven.evaluate(i*dx*PI/2);
        double expect = 1.0 - i*dx;
        if (abs(value - expect) > 1E-6){
            Py_RETURN_FALSE;
        }
    }

    Py_RETURN_TRUE;
}

PyObject* test_order2(){
    Vandeven vandeven(2);
    double dx = 1.0/9.0;
    for (unsigned int i=0;i<10;i++){
        double x = i*dx;
        double value = vandeven.evaluate(x*PI/2);
        double expect = 1.0  - 6*(x*x/2 - x*x*x/3);

        if (abs(value - expect) > 1E-4){
            Py_RETURN_FALSE;
        }
    }
    Py_RETURN_TRUE;
}