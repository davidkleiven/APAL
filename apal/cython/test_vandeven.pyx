cdef extern from "test_vandeven.hpp":
    object test_order1()
    object test_order2()

def pytest_test_order1():
    return test_order1()

def pytest_test_order2():
    return test_order2()
