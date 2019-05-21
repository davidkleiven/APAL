#include "additional_tools.hpp"
#include <stdexcept>
#include <sstream>

using namespace std;

int kronecker(int i, int j)
{
  if (i==j) return 1;
  return 0;
};

PyObject* string2py(const string &str)
{
  #if PY_MAJOR_VERSION >= 3
    // Python 3
    return PyUnicode_FromString(str.c_str());
  #else
    // Python 2
    return PyUnicode_FromString(str.c_str());
  #endif
}

string py2string(PyObject *str)
{
  #if PY_MAJOR_VERSION >= 3
    // Python 3
    return PyUnicode_AsUTF8(str);
  #else
    // Python 2
    return PyString_AsString(str);
  #endif
}

PyObject *int2py(int integer)
{
  #if PY_MAJOR_VERSION >= 3
    return PyLong_FromLong(integer);
  #else
    return PyInt_FromLong(integer);
  #endif
}

int py2int(PyObject *integer)
{
  #if PY_MAJOR_VERSION >= 3
    return PyLong_AsLong(integer);
  #else
    return PyInt_AsLong(integer);
  #endif
}

PyObject* get_attr(PyObject* obj, const char* name)
{
  PyObject* attr = PyObject_GetAttrString(obj, name);
  if (attr == nullptr)
  {
    stringstream ss;
    ss << "Python object has not attribute " << name;
    throw invalid_argument(ss.str());
  }
  return attr;
}

unsigned int list_size(PyObject *list)
{
  if (!PyList_Check(list))
  {
    throw invalid_argument("Python object is not a list. Cannot retrieve the length!");
  }
  return PyList_Size(list);
}
