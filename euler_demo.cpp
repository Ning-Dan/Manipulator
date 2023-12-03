#include <iostream>
#include "kdl/chain.hpp"
#include "kdl/chainiksolver.hpp"
#include "kdl/chainiksolverpos_lma.hpp"
#include "kdl/chainfksolverpos_recursive.hpp"

#include "eigen3/Eigen/Dense"

using namespace KDL;
#define pi 3.1415926

int main() {
    // 定义旋转矩阵 R
    Eigen::Matrix3d R;
    R <<     0.732, -0.661, 0.168,
         0.661, 0.732, -0.168,
        -0.168, 0.168, 0.970;

    // 将旋转矩阵转换为欧拉角
    Eigen::Vector3d eul = R.eulerAngles(2, 1, 0); // ZYX顺序
    // 显示欧拉角
    // eul[0] = eul[0]-pi;
    // eul[1] = -eul[1]-pi;
    // eul[2] = eul[2]+pi;
    std::cout << "欧拉角 (ZYX):\n" << eul.transpose()*180/pi << std::endl;
    return 0;
}

