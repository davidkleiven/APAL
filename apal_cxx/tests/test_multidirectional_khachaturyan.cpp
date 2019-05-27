#include "test_multidirectional_khachaturyan.hpp"
#include "multidirectional_khachaturyan.hpp"
#include "tools.hpp"
#include <iostream>

using namespace std;

const unsigned int L = 128;

PyObject *field2npy(MMSP::grid<2, MMSP::vector<fftw_complex> > &gr, unsigned int field){
    npy_intp dims[2] = {L, L};
    PyObject* npy_array = PyArray_SimpleNew(2, dims, NPY_DOUBLE);

    for (unsigned int node=0;node<MMSP::nodes(gr);node++){
        MMSP::vector<int> pos = gr.position(node);
        double* val = static_cast<double*>(PyArray_GETPTR2(npy_array, pos[0], pos[1]));
        *val = real(gr(node)[field]);
    }
    return npy_array;
}

PyObject* logger2dict_func_deriv(MultidirectionalKhachDataLogger<2> logger){
    PyObject* dict = PyDict_New();

    PyDict_SetItemString(dict, "shape_squared_in", field2npy(*logger.shape_squared_in, 0));
    PyDict_SetItemString(dict, "ft_shape_real", field2npy(*logger.fourier_shape_squared, 0));
    PyDict_SetItemString(dict, "b_tensor_dot_ft_squared", field2npy(*logger.b_tensor_dot_ft_squared, 0));
    PyDict_SetItemString(dict, "misfit_energy_contrib", field2npy(*logger.misfit_energy_contrib, 0));
    return dict;
}

PyObject* test_functional_derivative(PyObject *elastic, PyObject *misfit, const vector<double> &values){
    Khachaturyan khach(2, elastic, misfit);

    MultidirectionalKhachDataLogger<2> logger;
    MultidirectionalKhachaturyan multi(0.8);
    multi.set_logger(logger);

    multi.add_model(khach, 0);

    
    MMSP::grid<2, MMSP::vector<fftw_complex> > gr(1, 0, L, 0, L);

    if (values.size() != MMSP::nodes(gr)){
        throw invalid_argument("Invalide size of passed numpy arra !");
    }
    // Insert field values
    for (int node=0;node<MMSP::nodes(gr);node++){
        real(gr(node)[0]) = values[node];
        imag(gr(node)[0]) = 0.0;
    }

    MMSP::grid<2, MMSP::vector<fftw_complex> > grid_out(gr);

    vector<int> shape_fields;
    shape_fields.push_back(0);

    multi.functional_derivative(gr, grid_out, shape_fields);

    PyObject* dict = logger2dict_func_deriv(logger);
    clean_up(logger);

    PyDict_SetItemString(dict, "func_deriv", field2npy(grid_out, 0));
    return dict;
}

PyObject *test_contract_tensors(PyObject* tensor1, PyObject *tensor2){
    array< array<double, 3>, 3> t1, t2;

    for (unsigned int i=0;i<3;i++){
        PyObject *l1 = PyList_GetItem(tensor1, i);
        PyObject *l2 = PyList_GetItem(tensor2, i);
        for (unsigned int j=0;j<3;j++){
            t1[i][j] = PyFloat_AsDouble(PyList_GetItem(l1, j));
            t2[i][j] = PyFloat_AsDouble(PyList_GetItem(l2, j));
        }
    }

    double contract = contract_tensors(t1, t2);
    return PyFloat_FromDouble(contract);
}

void py2mat3x3(PyObject *obj, mat3x3 &mat){
    for (unsigned int i=0;i<3;i++){
        PyObject *list = PyList_GetItem(obj, i);
        for (unsigned int j=0;j<3;j++){
            mat[i][j] = PyFloat_AsDouble(PyList_GetItem(list, j));
        }
    }
}

PyObject *test_B_tensor_element(const vector<double> &vec, PyObject *pygf, PyObject *tensor1, PyObject *tensor2){
    mat3x3 t1, t2, gf;
    py2mat3x3(pygf, gf);
    py2mat3x3(tensor1, t1);
    py2mat3x3(tensor2, t2);

    MMSP::vector<double> mmsp_vec(3);
    memcpy(&(mmsp_vec[0]), &(vec[0]), 3*sizeof(double));

    double element = B_tensor_element(mmsp_vec, gf, t1, t2);
    return PyFloat_FromDouble(element);
}


PyObject* test_strain_energy_sphere(PyObject *elastic, PyObject *misfit){

    unsigned int size = 256;
    Khachaturyan khach(3, elastic, misfit);

    MultidirectionalKhachaturyan multi(0.8);
    multi.add_model(khach, 0);

    MMSP::grid<3, MMSP::vector<fftw_complex> > gr(1, 0, size, 0, size, 0, size);
    MMSP::grid<3, MMSP::vector<fftw_complex> > gr_squared_norm(1, 0, size, 0, size, 0, size);

    // Fill array
    double volume = 0.0;
    unsigned int R = 20;
    int center = size/2;
    for (int node=0;node<MMSP::nodes(gr);node++){
        MMSP::vector<int> pos = gr.position(node);
        double r_sq = 0.0;
        for (unsigned int i=0;i<3;i++){
            r_sq += pow(pos[i] - center, 2);
        }
        double r = sqrt(r_sq);

        real(gr(node)[0]) = r <= R ? 0.8 : 0.0;
        imag(gr(node)[0]) = 0.0;

        real(gr_squared_norm(node)[0]) = r <= R ? 1.0 : 0.0;
        imag(gr_squared_norm(node)[0]) = 0.0;
        volume += r <= R ? 1 : 0;
    }

    // Calculate the functional derivative which triggers an energy calculation
    MMSP::grid<3, MMSP::vector<fftw_complex> > grid_out(gr);

    vector<int> shape_fields;
    shape_fields.push_back(0);

    double misfit_contrib = multi.misfit_contribution(gr_squared_norm, shape_fields);

    multi.functional_derivative(gr, grid_out, shape_fields);

    PyObject *result_dict = PyDict_New();
    PyDict_SetItemString(result_dict, "misfit_contrib", PyFloat_FromDouble(misfit_contrib));
    PyDict_SetItemString(result_dict, "radius", PyFloat_FromDouble(R));
    PyDict_SetItemString(result_dict, "volume", PyFloat_FromDouble(volume));
    double energy = multi.get_last_strain_energy();
    PyDict_SetItemString(result_dict, "energy", PyFloat_FromDouble(energy));

    return result_dict;
}