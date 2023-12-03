#include "KineSolver.h"

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <memory>

using namespace std;

double get_random_num()
{
    return 2 * M_PI * (rand() / (std::pow(2, 31))) - M_PI;
}

int main()
{
    // set DH table
    std::array<double, 6> a = {0, 0.78, 0.2, 0, 0, 0};
    std::array<double, 6> alpha = {M_PI / 2, 0, M_PI / 2, -M_PI / 2, M_PI / 2, 0};
    std::array<double, 6> d = {0, 0, 0, 1.08, 0, 0.1};
    std::array<double, 6> theta = {0, M_PI / 2, 0, 0, 0, 0};

    // set jnts limit
    std::array<double, 6> upper_limit{170 * DEG_2_RAD, 160 * DEG_2_RAD, 311.9 * DEG_2_RAD, 200 * DEG_2_RAD, 140 * DEG_2_RAD, 270 * DEG_2_RAD};
    std::array<double, 6> lower_limit{-170 * DEG_2_RAD, -90 * DEG_2_RAD, -135.1 * DEG_2_RAD, -200 * DEG_2_RAD, -140 * DEG_2_RAD, -270 * DEG_2_RAD};

    // init kine solver
    std::unique_ptr<KineSolver> kine_ptr = std::make_unique<KineSolver>(a, alpha, d, theta);
    kine_ptr->SetJntLimit(lower_limit, upper_limit);

    // utest
    for (int i = 0; i < 1000; ++i)
    {
        Eigen::VectorXd jnts(6);
        std::cout << "+++++++++++" << std::endl;
        // set random jnts
        jnts << get_random_num(), get_random_num(), get_random_num(), get_random_num(), get_random_num(), get_random_num();
        // jnts << -0.1, 0.2, -0.3, 0.4, -0.5, 0.6;
        std::cout << "jnts:" << jnts.transpose() << std::endl;
        // solve fk
        Eigen::Isometry3d tf = kine_ptr->ComputeFK(jnts);

        std::vector<Eigen::VectorXd> res_vec;
        // solve ik
        kine_ptr->ComputeAllIK(tf, &res_vec);
        // check whether 8 results were returned,if not,break this loop and print waring
        if (res_vec.size() < 8)
        {
            std::cout << "[warning] Ik failed:" << res_vec.size() << "  i:" << i << std::endl;
            break;
        }
        // Eigen::VectorXd init_joint(6), solved_jnts;
        // init_joint << 0.2, 0.2, 0.2, 0.2, 0.2, 0.2;
        // kine_ptr->ComputeIK(tf, init_joint, &solved_jnts);
        // std::cout << "solved_jnts:" << solved_jnts.transpose() << std::endl;
        // Eigen::Isometry3d tf1 = kine_ptr->ComputeFK(solved_jnts);
        // std::cout << "tf1:" << tf.matrix() << std::endl;
        // std::cout << "tf2:" << tf1.matrix() << std::endl;
        std::cout << "__________________" << std::endl;
    }

    return 0;
}