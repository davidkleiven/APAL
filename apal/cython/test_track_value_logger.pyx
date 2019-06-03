from libcpp cimport bool as bool_t

cdef extern from "test_track_value_logger.hpp":
    bool_t read_from_file()
    void init_keys_from_entry()
    void track_values_append()

def pyread_from_file():
    return read_from_file()

def pyinit_keys_from_entry():
    init_keys_from_entry()

def pytrack_values_append():
    track_values_append()