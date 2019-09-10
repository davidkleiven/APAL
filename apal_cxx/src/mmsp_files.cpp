// Hack: We need to include all cpp files here
// because the header files of the MMSP project
// contains implementations
#include "tools.cpp"
#include "phase_field_simulation.cpp"
#include "cahn_hilliard_phasefield.cpp"
#include "chgl.cpp"
#include "multidirectional_khachaturyan.cpp"
#include "fftw_mmsp.cpp"
#include "chgl_realspace.cpp"
#include "conjugate_gradient.cpp"
#include "two_phase_landau.cpp"
#include "chc_noise.cpp"
#include "elasticity_util.cpp"
#include "origin_singularity_integration.cpp"
#include "raised_cosine.cpp"
#include "gaussian_filter.cpp"
#include "vandeven.cpp"

// Test files
#include "test_multidirectional_khachaturyan.cpp"
#include "test_build_matrices.cpp"
#include "test_fftw.cpp"
#include "test_track_value_logger.cpp"
#include "test_raised_cosine_filter.cpp"
#include "test_vandeven.cpp"