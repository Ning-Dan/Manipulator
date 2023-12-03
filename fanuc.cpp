#include <iostream>
#include "kdl/chain.hpp"
#include "kdl/chainiksolver.hpp"
#include "kdl/chainiksolverpos_lma.hpp"
#include "kdl/chainiksolverpos_nr.hpp"
#include "kdl/chainfksolverpos_recursive.hpp"
#include "kdl/utilities/utility.h"

#include "eigen3/Eigen/Dense"

using namespace KDL;

void NewSegment(double a, double alpha, double d, double theta, Chain &chain);

Chain faunc()
{
    Chain robot;
    NewSegment(0, 0, 0.0, 0, robot);
    NewSegment(150.0, -PI / 2, 0.0, -PI / 2, robot);
    NewSegment(780.0, 0, 0.0, 0.0, robot);
    NewSegment(200.0, PI / 2, -1080, 0.0, robot);
    NewSegment(0.0, -PI / 2, 0.0, PI, robot);
    NewSegment(0.0, -PI / 2, -100, 0.0, robot);
    robot.addSegment(Segment(Joint(Joint::None),
                             Frame::DH_Craig1989(0, -PI, 0, -PI)));
    return robot;
}

void NewSegment(double a, double alpha, double d, double theta, Chain &chain)
{
    chain.addSegment(Segment(Joint(Joint::None),
                             Frame::DH_Craig1989(a, alpha, d, theta)));

    chain.addSegment(KDL::Segment(KDL::Joint(KDL::Joint::RotZ)));
}

void print_frame(Frame pos_goal)
{
    Eigen::Isometry3d tf;
    tf.translation() << pos_goal.p.data[0], pos_goal.p.data[1], pos_goal.p.data[2];
    tf.linear() << pos_goal.M.data[0], pos_goal.M.data[1], pos_goal.M.data[2],
        pos_goal.M.data[3], pos_goal.M.data[4], pos_goal.M.data[5],
        pos_goal.M.data[6], pos_goal.M.data[7], pos_goal.M.data[8];

    std::cout << "tf translation:" << tf.translation().transpose() << std::endl;
    auto euler = tf.linear().eulerAngles(2, 1, 0);
    std::cout << "euler angle:" << euler * 180 / PI << std::endl;
}

int main()
{
    auto robot = faunc();
    int n = robot.getNrOfJoints();
    KDL::ChainFkSolverPos_recursive fk_solver(robot); // init fk
    KDL::JntArray q_given(n);
    q_given.data << 10, 20, -30, 40, -50, 60;
    q_given.data = q_given.data * M_PI / 180.0;

    KDL::Frame pos_goal;
    fk_solver.JntToCart(q_given, pos_goal);

    print_frame(pos_goal);

    KDL::JntArray q_init(n);
    q_init.data.setZero();

    int retval;
    KDL::JntArray q_sol(n);

    Frame p_out = Frame::Identity();

    Eigen::Matrix<double, 6, 1> L;
    L(0) = 1;
    L(1) = 1;
    L(2) = 1;
    L(3) = 0.01;
    L(4) = 0.01;
    L(5) = 0.01;
    ChainIkSolverPos_LMA ik_solver(robot,L,1e-10,8000);
    retval = ik_solver.CartToJnt(q_init, pos_goal, q_sol);

    std::cout << "ik return"
              << ":" << retval << std::endl;
    std::cout << "q_given"
              << " :" << q_given.data.transpose() / M_PI * 180.0 << std::endl;
    std::cout << "q_sol  "
              << " :" << q_sol.data.transpose() / M_PI * 180.0 << std::endl;

    fk_solver.JntToCart(q_sol, pos_goal);
    print_frame(pos_goal);
}
