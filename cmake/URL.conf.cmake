
# Here we declare the different PATH, TAG and PATCH to get the HDB5_IO dependencies


# MathUtils
set(mathutils_URL git@frydom-ce.org:ce/mathutils.git)
set(mathutils_TAG v1.14.1 CACHE STRING "MathUtils version")
set(MATHUTILS_BUILD_TESTS OFF CACHE BOOL "")
set(MATHUTILS_BUILD_BOOST_TESTS OFF CACHE BOOL "")


# HDF5
set(hdf5_URL https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5-1_13_1.tar.gz)
set(hdf5_TAG 1.13.1)


# HighFive
set(highfive_URL https://github.com/BlueBrain/HighFive.git)
set(highfive_TAG v2.4.1)
set(highfive_PATCH highfive.patch CACHE STRING "highfive version")


# GoogleTest
set(googletest_URL https://github.com/google/googletest.git)
set(googletest_TAG release-1.10.0 CACHE STRING "googletest version")


# meshoui
set(meshoui_URL git@frydom-ce.org:ce/meshoui.git)
set(meshoui_TAG v1.3.1 CACHE STRING "meshoui version")

# VTK
set(VTK_TAG 8.2 CACHE STRING "VTK version")
