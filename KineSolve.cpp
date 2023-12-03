#include "KineSolver.h"

#include <cmath>
#include <iostream>

KineSolver::KineSolver(const std::array<double, 6> &a, const std::array<double, 6> &alpha, const std::array<double, 6> &d, const std::array<double, 6> &theta)
    : a_(a), alpha_(alpha), d_(d), theta_(theta)
{
}

KineSolver::~KineSolver()
{
}

Eigen::Isometry3d KineSolver::DH(std::size_t id, double theta, bool is_deg)
{
    if (is_deg)
    {
        theta *= DEG_2_RAD;
    }
    theta += theta_.at(id);
    double d = d_.at(id);
    double a = a_.at(id);
    double alpha = alpha_.at(id);

    Eigen::Isometry3d tf = Eigen::Isometry3d::Identity();
    tf.linear().matrix() << cos(theta), -sin(theta) * cos(alpha), sin(theta) * sin(alpha),
        sin(theta), cos(theta) * cos(alpha), -cos(theta) * sin(alpha),
        0, sin(alpha), cos(alpha);
    tf.translation() << a * cos(theta), a * sin(theta), d;
    return tf;
}

Eigen::Isometry3d KineSolver::ComputeFK(const Eigen::VectorXd &jnts, bool is_deg)
{
    Eigen::Isometry3d tf = Eigen::Isometry3d::Identity();
    if (jnts.size() != dof_)
    {
        std::cout << "[Warning] Input jnts size should be equal with dof" << std::endl;
        return tf;
    }
    for (std::size_t i = 0; i < dof_; i++)
    {
        tf = tf * DH(i, jnts[i], is_deg);
    }
    return tf;
}

void KineSolver::ComputeAllIK(const Eigen::Isometry3d &tf, std::vector<Eigen::VectorXd> *res_vec)
{
    count = 0;
    // calculate Pw
    Eigen::Vector3d p_e = tf.translation();
    Eigen::Vector3d a_e = tf.linear().col(2);
    double d6 = d_.at(5);
    p_w = p_e - d6 * a_e;

    // std::cout << "a_e:" << a_e.transpose() << std::endl;
    // std::cout << "pw:" << p_w.transpose() << std::endl;

    // calculate theta1/theta2/theta3
    double a2 = a_.at(1);
    double a1 = a_.at(0);
    double a3 = sqrt(std::pow(d_.at(3), 2) + std::pow(a_.at(2), 2));
    double belta = std::atan2(a_.at(2), d_.at(3));

    // solve theta1
    double theta1_1;
    double theta1_2;
    if (p_w[1] >= 0)
    {
        theta1_1 = atan2(p_w[1], p_w[0]) - theta_.at(0);
        theta1_2 = atan2(-p_w[1], -p_w[0]) - theta_.at(0);
    }
    else
    {
        theta1_2 = atan2(p_w[1], p_w[0]) - theta_.at(0);
        theta1_1 = atan2(-p_w[1], -p_w[0]) - theta_.at(0);
    }

    // solve theta3
    double s3 = (p_w[0] * p_w[0] + p_w[1] * p_w[1] + p_w[2] * p_w[2] - std::pow(a2, 2) - std::pow(a3, 2)) / (2 * a2 * a3);
    double c3 = sqrt(1 - std::pow(s3, 2));
    double theta3_1 = atan2(s3, c3) - theta_.at(2) - belta;
    double theta3_2 = atan2(s3, -c3) - theta_.at(2) - belta;

    // theta3_1 += M_PI / 2;
    // theta3_2 += M_PI / 2;
    LimitJoints(theta3_1);
    LimitJoints(theta3_2);

    std::vector<double> vec1 = SolveTrigonometricEquation(-a3 * cos(theta3_1 + belta), a2 + a3 * sin(theta3_1 + belta), p_w[2]);
    double theta2_1 = vec1[0] - theta_.at(1);
    double theta2_2 = vec1[1] - theta_.at(1);

    std::vector<double> vec2 = SolveTrigonometricEquation(-a3 * cos(theta3_2 + belta), a2 + a3 * sin(theta3_2 + belta), p_w[2]);
    double theta2_3 = vec2[0] - theta_.at(1);
    double theta2_4 = vec2[1] - theta_.at(1);

    LimitJoints(theta2_1);
    LimitJoints(theta2_2);
    LimitJoints(theta2_3);
    LimitJoints(theta2_4);

    std::cout << "theta1_1:" << theta1_1 << std::endl;
    std::cout << "theta1_2:" << theta1_2 << std::endl;

    std::cout << "theta2_1:" << theta2_1 << std::endl;
    std::cout << "theta2_2:" << theta2_2 << std::endl;
    std::cout << "theta2_3:" << theta2_3 << std::endl;
    std::cout << "theta2_4:" << theta2_4 << std::endl;

    std::cout << "theta3_1:" << theta3_1 << std::endl;
    std::cout << "theta3_2:" << theta3_2 << std::endl;

    // Solve sphere wrist

    SolveSphereWrist(theta1_1, theta2_1, theta3_1, tf, res_vec);
    SolveSphereWrist(theta1_1, theta2_2, theta3_1, tf, res_vec);

    SolveSphereWrist(theta1_2, theta2_2, theta3_1, tf, res_vec);
    SolveSphereWrist(theta1_2, theta2_1, theta3_1, tf, res_vec);

    SolveSphereWrist(theta1_1, theta2_3, theta3_2, tf, res_vec);
    SolveSphereWrist(theta1_1, theta2_4, theta3_2, tf, res_vec);

    SolveSphereWrist(theta1_2, theta2_4, theta3_2, tf, res_vec);
    SolveSphereWrist(theta1_2, theta2_3, theta3_2, tf, res_vec);

    // std::cout << "ik size:" << res_vec->size() << std::endl;
}

void KineSolver::SolveSphereWrist(double theta1, double theta2, double theta3, const Eigen::Isometry3d &tf, std::vector<Eigen::VectorXd> *res_vec)
{
    count++;
    Eigen::Matrix3d R_0_3, R;
    R_0_3 = (DH(0, theta1) * DH(1, theta2) * DH(2, theta3)).linear();
    if (((DH(0, theta1) * DH(1, theta2) * DH(2, theta3) * DH(3, 0)).translation() - p_w).norm() > 1e-3)
    {
        return;
    }
    // std::cout << "theta:" << theta1 << "  " << theta2 << "  " << theta3 << std::endl;
    // std::cout << "1-3:" << (DH(0, theta1) * DH(1, theta2) * DH(2, theta3) * DH(3, 0)).translation().transpose() << std::endl;

    R = R_0_3.transpose() * tf.linear();

    double theta4_1 = std::atan2(R(1, 2), (R(0, 2)) - theta_.at(3));
    double theta5_1 = std::atan2((sqrt(R(0, 2) * R(0, 2) + R(1, 2) * R(1, 2))), R(2, 2)) - theta_.at(4);
    double theta6_1 = std::atan2(R(2, 1), -R(2, 0)) - theta_.at(5);

    double theta4_2 = std::atan2(-R(1, 2), (-R(0, 2)) - theta_.at(3));
    double theta5_2 = std::atan2((-sqrt(R(0, 2) * R(0, 2) + R(1, 2) * R(1, 2))), R(2, 2)) - theta_.at(4);
    double theta6_2 = std::atan2(-R(2, 1), R(2, 0)) - theta_.at(5);
    LimitJoints(theta4_1);
    LimitJoints(theta4_2);
    LimitJoints(theta5_1);
    LimitJoints(theta6_1);
    LimitJoints(theta5_2);
    LimitJoints(theta6_2);
    // std::cout << "theta4:" << theta4 << std::endl;
    // std::cout << "theta5_1:" << theta5_1 << std::endl;
    // std::cout << "theta6_1:" << theta6_1 << std::endl;

    // std::cout << "theta5_2:" << theta5_2 << std::endl;
    // std::cout << "theta6_2:" << theta6_2 << std::endl;

    Eigen::VectorXd res1(6), res2(6);
    res1 << theta1, theta2, theta3, theta4_1, theta5_1, theta6_1;
    res2 << theta1, theta2, theta3, theta4_2, theta5_2, theta6_2;

    Eigen::Quaterniond q1(ComputeFK(res1, false).linear());
    Eigen::Quaterniond q2(ComputeFK(res2, false).linear());
    Eigen::Quaterniond q(tf.linear());
    if ((q1.coeffs() - q.coeffs()).norm() < 1e-5)
    {
        res_vec->emplace_back(res1);
        std::cout << "res:" << res1.transpose() << "  " << count << std::endl;
    }
    if ((q2.coeffs() - q.coeffs()).norm() < 1e-5)
    {
        res_vec->emplace_back(res2);
        std::cout << "res:" << res2.transpose() << "  " << count << std::endl;
    }
    // std::cout << "=====================================" << std::endl;

    // std::cout << "R" << R << std::endl;
    // std::cout << "1:" << -sin(theta5_1) * cos(theta6_1) << "  1:" << sin(theta5_1) * sin(theta6_1) << std::endl;
    // std::cout << "2:" << -sin(theta5_2) * cos(theta6_2) << "  2:" << sin(theta5_2) * sin(theta6_2) << std::endl;

    // std::cout << "res1 fk:" << this->ComputeFK(res1, false).matrix() << '\n'
    //            << std::endl;
    // std::cout << "res2 fk:" << this->ComputeFK(res2, false).matrix() << '\n'
    //            << std::endl
    //           << std::endl;

    // std::cout << "res1 fk:" << (DH(3, theta4_1) * DH(4, theta5_1) * DH(5, theta6_1)).linear() << '\n'
    //           << "res1 j:" << res1.transpose() << std::endl;
    // std::cout << "res2 fk:" << (DH(3, theta4_2) * DH(4, theta5_2) * DH(5, theta6_2)).linear() << '\n'
    //           << "res2 j:" << res2.transpose() << std::endl
    //           << std::endl;
}

IKStatus KineSolver::ComputeIK(const Eigen::Isometry3d &tf, const Eigen::VectorXd &init_jnt, Eigen::VectorXd *jnts)
{
    if (init_jnt.size() != dof_)
    {
        return IKStatus::kSolveFail;
    }

    std::vector<Eigen::VectorXd> res_vec;
    double min_dis = std::numeric_limits<double>::max();
    ComputeAllIK(tf, &res_vec);
    for (auto res : res_vec)
    {
        double dis = (init_jnt - res).norm();
        if (dis < min_dis)
        {
            min_dis = dis;
            (*jnts) = res;
        }
    }
    if (IsJointWithinLimit(*jnts) == IKStatus::kExceedJointLimit)
    {
        return IKStatus::kExceedJointLimit;
    }
    return IKStatus::kSuccess;
}

IKStatus KineSolver::IsJointWithinLimit(const Eigen::VectorXd &jnt)
{
    for (auto i = 0; i < jnt_lower_limit_.size(); i++)
    {
        if (jnt(i) > jnt_upper_limit_.at(i) || jnt(i) < jnt_lower_limit_.at(i))
        {
            return IKStatus::kExceedJointLimit;
        }
    }
    return IKStatus::kSuccess;
}

void KineSolver::LimitJoints(double &joint)
{
    if (joint > M_PI)
    {
        joint -= 2 * M_PI;
    }
    else if (joint < -M_PI)
    {
        joint += 2 * M_PI;
    }
}

std::vector<double> KineSolver::SolveTrigonometricEquation(double b, double a, double c)
{
    // solve equation: b*cos(theta)+a*sin(theta) = c
    double r = sqrt(a * a + b * b);
    double theta1 = atan2(c, sqrt(r * r - c * c)) - atan2(b, a);
    double theta2 = atan2(c, -sqrt(r * r - c * c)) - atan2(b, a);

    std::vector<double> res;
    res.emplace_back(theta1);
    res.emplace_back(theta2);

    return res;
}

void KineSolver::SetJntLimit(const std::array<double, 6> &lower, const std::array<double, 6> &upper)
{
    jnt_lower_limit_ = lower;
    jnt_upper_limit_ = upper;
}