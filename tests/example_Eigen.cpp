#include <iostream>
#include <Eigen/Eigen>
#include <highfive/H5Easy.hpp>

int main()
{
  {
    Eigen::Vector3f vec(1,2,3);
    H5Easy::File file("example.h5", H5Easy::File::Overwrite);

    H5Easy::dump(file, "/path/to/A", vec);

    auto a = H5Easy::load<Eigen::Vector3f>(file, "/path/to/A");
    std::cout << a << std::endl;
  }
  return 0;
}



//#include <highfive/H5Easy.hpp>
//#include <eigen3/Eigen/Eigen>
//
//int main()
//{
//  H5Easy::File file("example.h5", H5Easy::File::Overwrite);
//
//  Eigen::MatrixXd vec_in = Eigen::MatrixXd::Random(10,5);
//
//  H5Easy::dump(file, "Eigen", vec_in);
//
//  Eigen::MatrixXd vec_out = H5Easy::load<Eigen::Vector3d>(file, "Eigen");
//
//  std::cout<<vec_out<<std::endl;
//}