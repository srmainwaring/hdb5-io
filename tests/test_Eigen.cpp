//
// Created by lletourn on 02/03/20.
//
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <typeinfo>
#include <vector>

#include <Eigen/Eigen>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5Utility.hpp>

#define BOOST_TEST_MAIN HighFiveTestBase
#include <boost/test/unit_test.hpp>

#include "tests_high_five.hpp"

using namespace HighFive;

int main() {

  const std::string FILE_NAME("test_eigen.h5");
  const std::string DS_NAME = "ds";

// Create a new file using the default property lists.
  File file(FILE_NAME, File::ReadWrite | File::Create | File::Truncate);

  auto test = [&DS_NAME, &file](const std::string &test_flavor, const auto &vec_input, auto &vec_output) {
    file.createDataSet(DS_NAME + test_flavor, vec_input).write(vec_input);
    file.getDataSet(DS_NAME + test_flavor).read(vec_output);
    BOOST_CHECK(vec_input == vec_output);
  };

  std::string DS_NAME_FLAVOR;


// std::vector<of vector <of POD>>
  {
    DS_NAME_FLAVOR = "VectorOfVectorOfPOD";
    std::vector <std::vector<float>> vec_in{{5.0f, 6.0f, 7.0f},
                                            {5.1f, 6.1f, 7.1f},
                                            {5.2f, 6.2f, 7.2f}};
    std::vector <std::vector<float>> vec_out;

    test(DS_NAME_FLAVOR, vec_in, vec_out);
  }

//std::vector<Eigen::Vector3d>
  {
    DS_NAME_FLAVOR = "VectorOfEigenVector3d";
    std::vector <Eigen::Vector3d> vec_in{{5.0, 6.0, 7.0},
                                         {7.0, 8.0, 9.0}};
    std::vector <Eigen::Vector3d> vec_out;
    test(DS_NAME_FLAVOR, vec_in, vec_out);
  }

// Eigen Vector2d
  {
    DS_NAME_FLAVOR = "EigenVector2d";
    Eigen::Vector2d vec_in{5.0, 6.0};
    Eigen::Vector2d vec_out;

    test(DS_NAME_FLAVOR, vec_in, vec_out);
  }

// Eigen Matrix
  {
    DS_NAME_FLAVOR = "EigenMatrix";
    Eigen::Matrix<double, 3, 3> vec_in;
    vec_in << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    Eigen::Matrix<double, 3, 3> vec_out;

    test(DS_NAME_FLAVOR, vec_in, vec_out);
  }

// Eigen MatrixXd
  {
    DS_NAME_FLAVOR = "EigenMatrixXd";
    Eigen::MatrixXd vec_in = 100. * Eigen::MatrixXd::Random(20, 5);
    Eigen::MatrixXd vec_out(20, 5);

    test(DS_NAME_FLAVOR, vec_in, vec_out);
  }

// std::vector<of EigenMatrixXd>
  {
    DS_NAME_FLAVOR = "VectorEigenMatrixXd";

    Eigen::MatrixXd m1 = 100. * Eigen::MatrixXd::Random(20, 5);
    Eigen::MatrixXd m2 = 100. * Eigen::MatrixXd::Random(20, 5);
    std::vector <Eigen::MatrixXd> vec_in;
    vec_in.push_back(m1);
    vec_in.push_back(m2);
    std::vector <Eigen::MatrixXd> vec_out(2, Eigen::MatrixXd::Zero(20, 5));

    test(DS_NAME_FLAVOR, vec_in, vec_out);
  }

// std::vector<of EigenMatrixXd> - exception
  {
    DS_NAME_FLAVOR = "VectorEigenMatrixXdExc";

    Eigen::MatrixXd m1 = 100. * Eigen::MatrixXd::Random(20, 5);
    Eigen::MatrixXd m2 = 100. * Eigen::MatrixXd::Random(20, 5);
    std::vector <Eigen::MatrixXd> vec_in;
    vec_in.push_back(m1);
    vec_in.push_back(m2);
    file.createDataSet(DS_NAME + DS_NAME_FLAVOR, vec_in).write(vec_in);

    std::vector <Eigen::MatrixXd> vec_out_exception;
    SilenceHDF5 silencer;
    BOOST_CHECK_THROW(file.getDataSet(DS_NAME + DS_NAME_FLAVOR).read(vec_out_exception), HighFive::DataSetException);
  }
}