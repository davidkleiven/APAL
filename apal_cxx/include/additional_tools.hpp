#ifndef ADDITIONAL_TOOLS_H
#define ADDITIONAL_TOOLS_H
#include <vector>
#include <iostream>
#include <map>
#include <set>
#include <array>
#include <Python.h>

//class SymbolChange;

template<class key,class value>
std::ostream& operator <<(std::ostream &out, const std::map<key,value> &map );

template<class T>
std::ostream& operator <<( std::ostream &out, const std::vector<T> &vec );

template<class T>
std::ostream& operator <<(std::ostream &out, const std::set<T> &set);

template<class T>
std::vector<T>& cyclic_permute( std::vector<T> &vec );

template<class T, unsigned int N>
std::ostream& operator <<(std::ostream &out, const std::array<T, N> &array);

template<class T>
void keys( std::map<std::string,T> &, std::vector<std::string> &keys );

template<class T>
void set2vector( const std::set<T> &s, std::vector<T> &v );

int kronecker(int i, int j);

PyObject* string2py(const std::string &string);
std::string py2string(PyObject *str);

PyObject* int2py(int integer);
int py2int(PyObject *integer);

PyObject *get_attr(PyObject *obj, const char* name);

/** Return the length of a python list */
unsigned int list_size(PyObject *list);

/** Return true if element is in vector*/
template<class T>
bool is_in_vector(const T &value, const std::vector<T> &vec);

template<class T>
void insert_in_set(const std::vector<T> &vec, std::set<T> &unique);

#include "additional_tools.tpp"
#endif
