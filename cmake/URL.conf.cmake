
# Here we declare the different PATH, TAG and PATCH to get the HDB5_IO dependencies


# MathUtils
set(mathutils_URL git@frydom-ce.org:ce/mathutils.git)
set(mathutils_TAG v1.5 CACHE STRING "MathUtils version")
set(MATHUTILS_BUILD_TESTS OFF CACHE BOOL "")
set(MATHUTILS_BUILD_BOOST_TESTS OFF CACHE BOOL "")


# HighFive
set(highfive_URL https://github.com/BlueBrain/HighFive.git)
set(highfive_TAG v2.2.2)
set(highfive_PATCH highfive.patch CACHE STRING "highfive version")


# GoogleTest
set(googletest_URL https://github.com/google/googletest.git)
set(googletest_TAG release-1.10.0 CACHE STRING "googletest version")


# MeshOui
set(meshoui_URL git@frydom-ce.org:ce/meshoui.git)
set(meshoui_TAG v1.0.2 CACHE STRING "meshoui version")


# HDF5
set(HDF5_URL https://github.com/HDFGroup/hdf5.git)
set(HDF5_TAG hdf5-1_10_6 CACHE STRING "HDF5 version")
set(HDF5_PATCH hdf5.patch)