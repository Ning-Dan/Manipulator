#include "KineSolver.h"

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <memory>

using namespace std;

double get_random_num()
{
    return M_PI * (rand() / (std::pow(2, 31))) - M_PI / 2;
}

void test_ik()
{
    // std::array<double, 7> a = {0.15, 0.78, 0.2, 0, 0, 0};
    std::array<double, 7> a = {0, 0.15, 0.78, 0.2, 0, 0, 0};
    std::array<double, 7> alpha = {0, -M_PI / 2, 0, M_PI / 2, -M_PI / 2, -M_PI / 2, -M_PI};
    std::array<double, 7> d = {0, 0, 0, -1.08, 0, -0.1, 0};
    std::array<double, 7> theta = {0, -M_PI / 2, 0, 0, M_PI, 0, -M_PI};
    // set jnts limit
    std::array<double, 6> upper_limit{170 * DEG_2_RAD, 160 * DEG_2_RAD, 311.9 * DEG_2_RAD, 200 * DEG_2_RAD, 140 * DEG_2_RAD, 270 * DEG_2_RAD};
    std::array<double, 6> lower_limit{-170 * DEG_2_RAD, -90 * DEG_2_RAD, -135.1 * DEG_2_RAD, -200 * DEG_2_RAD, -140 * DEG_2_RAD, -270 * DEG_2_RAD};

    // init kine solver
    std::unique_ptr<KineSolver> kine_ptr = std::make_unique<KineSolver>(a, alpha, d, theta);
    kine_ptr->SetJntLimit(lower_limit, upper_limit);

    Eigen::VectorXd res_jnts(6);
    for (int i = 0; i < 1000; ++i)
    {
        Eigen::VectorXd jnts(6);
        std::cout << "+++++++++++" << std::endl;
        jnts << get_random_num(), get_random_num(), get_random_num(), get_random_num(), get_random_num(), get_random_num();
        std::cout << "jnts:" << jnts.transpose() * 180 / M_PI << std::endl;

        Eigen::Isometry3d tf = kine_ptr->ComputeFK(jnts, false);
        // std::cout << "tf:" << tf.matrix() << std::endl;

        IKStatus res = kine_ptr->ComputeIK(tf, jnts, &res_jnts);
        std::cout << "res_jnts:" << res_jnts.transpose() * 180 / M_PI << std::endl;
        if ((jnts - res_jnts).norm() > 1e-3)
        {
            std::cout << "[warning] Ik failed at:" << i << std::endl;
            break;
        }

        std::cout << "__________________  " << i << std::endl;
    }
}

void test_kine_class()
{
    // set DH table
    std::array<double, 7> a = {0, 0.78, 0, 0, 0, 0, 0};
    std::array<double, 7> alpha = {M_PI / 2, 0, M_PI / 2, -M_PI / 2, M_PI / 2, 0, -M_PI};
    std::array<double, 7> d = {0, 0, 0, 1.08, 0, 0.1, 0};
    std::array<double, 7> theta = {0, M_PI / 2, 0, 0, 0, 0, -M_PI};

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

        std::cout << "jnts:" << jnts.transpose() << std::endl;
        // solve fk
        Eigen::Isometry3d tf = kine_ptr->ComputeFK(jnts);
        std::vector<Eigen::VectorXd> res_vec;
        // solve ik
        // kine_ptr->ComputeAllIK(tf, &res_vec);

        // std::cout << "res_vec:" << res_vec.size() << std::endl;

        // // check whether 8 results were returned,if not,break this loop and print waring
        // if (res_vec.size() < 4)
        // {
        //     std::cout << "[warning] Ik failed:" << res_vec.size() << "  i:" << i << std::endl;
        //     break;
        // }
        Eigen::VectorXd res_jnts(6);

        kine_ptr->ComputeIK(tf, jnts, &res_jnts);
        std::cout << "res_jnts:" << res_jnts.transpose() << std::endl;

        if ((jnts - res_jnts).norm() > 1e-3)
        {
            std::cout << "[warning] Ik failed:" << res_vec.size() << "  i:" << i << std::endl;
            break;
        }

        std::cout << "__________________  " << i << std::endl;
    }
}

int main()
{
    test_ik();
}