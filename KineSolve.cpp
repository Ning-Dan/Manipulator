#include "KineSolver.h"

#include <cmath>
#include <iostream>

KineSolver::KineSolver(const std::array<double, 7> &a, const std::array<double, 7> &alpha, const std::array<double, 7> &d, const std::array<double, 7> &theta)
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

Eigen::Isometry3d KineSolver::mDH(std::size_t id, double theta, bool is_deg)
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
    tf.linear().matrix() << cos(theta), -sin(theta), 0,
        sin(theta) * cos(alpha), cos(theta) * cos(alpha), -sin(alpha),
        sin(theta) * sin(alpha), cos(theta) * sin(alpha), cos(alpha);
    tf.translation() << a, -d * sin(alpha), d * cos(alpha);
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
        tf = tf * mDH(i, jnts[i], is_deg);
    }
    tf = tf * mDH(6, 0, false);
    return tf;
}

Eigen::Isometry3d KineSolver::ComputeBaseEndfTf(const Eigen::VectorXd &jnts, bool is_deg)
{
    Eigen::Isometry3d tf = Eigen::Isometry3d::Identity();
    if (jnts.size() != dof_)
    {
        std::cout << "[Warning] Input jnts size should be equal with dof" << std::endl;
        return tf;
    }
    for (std::size_t i = 0; i < dof_; i++)
    {
        tf = tf * mDH(i, jnts[i], is_deg);
    }
    return tf;
}

void KineSolver::ComputeAllIK(const Eigen::Isometry3d &tf, std::vector<Eigen::VectorXd> *res_vec)
{
    Eigen::Isometry3d tf1(tf);
    tf1 = tf * mDH(6, 0).inverse();

    // calculate Pw
    Eigen::Vector3d p_e = tf1.translation();
    Eigen::Vector3d a_e = tf1.linear().col(2);
    double d6 = d_.at(5);
    p_w = p_e - d6 * a_e;
    // std::cout << "p_w:" << p_w.transpose() << std::endl;

    // calculate theta1/theta2/theta3
    double a2 = a_.at(1);
    belta_ = std::atan(a_.at(3) / d_.at(3));

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
    // std::cout << "theta1_1:" << theta1_1 * 180 / M_PI << std::endl;
    // std::cout << "theta1_2:" << theta1_2 * 180 / M_PI << std::endl;

    double x1 = p_w[0] - a2 * cos(theta1_1);
    double y1 = p_w[1] - a2 * sin(theta1_1);

    double x2 = p_w[0] - a2 * cos(theta1_2);
    double y2 = p_w[1] - a2 * sin(theta1_2);

    std::vector<double> theta1_vec, theta2_vec, theta3_vec;
    theta1_vec.emplace_back(theta1_1);
    theta1_vec.emplace_back(theta1_2);

    SolveTheta3(x1, y1, &theta3_vec);
    SolveTheta3(x2, y2, &theta3_vec);

    for (auto theta3 : theta3_vec)
    {
        SolveTheta2(theta3, &theta2_vec);
    }
    for (auto theta1 : theta1_vec)
    {
        for (auto theta2 : theta2_vec)
        {
            for (auto theta3 : theta3_vec)
            {
                SolveSphereWrist(theta1, theta2, theta3, tf1, res_vec);
            }
        }
    }
    // std::cout << "res vec size:" << res_vec->size() << std::endl;
    // int size_theta3 = theta3_vec.size();
    // if (size_theta3 == 2)
    // {
    //     double theta3_1 = theta3_vec[0];
    //     double theta3_2 = theta3_vec[1];

    //     std::vector<double> theta2_vec;
    //     SolveTheta2(theta3_1, &theta2_vec);
    //     SolveTheta2(theta3_2, &theta2_vec);
    //     double theta2_1 = theta2_vec[0];
    //     double theta2_2 = theta2_vec[1];
    //     double theta2_3 = theta2_vec[2];
    //     double theta2_4 = theta2_vec[3];

    //     std::cout << "theta2_1:" << theta2_1 * 180 / M_PI << std::endl;
    //     std::cout << "theta2_2:" << theta2_2 * 180 / M_PI << std::endl;
    //     std::cout << "theta2_3:" << theta2_3 * 180 / M_PI << std::endl;
    //     std::cout << "theta2_4:" << theta2_4 * 180 / M_PI << std::endl;

    //     SolveSphereWrist(theta1_1, theta2_1, theta3_1, tf1, res_vec);
    //     SolveSphereWrist(theta1_1, theta2_2, theta3_1, tf1, res_vec);

    //     SolveSphereWrist(theta1_2, theta2_2, theta3_1, tf1, res_vec);
    //     SolveSphereWrist(theta1_2, theta2_1, theta3_1, tf1, res_vec);

    //     SolveSphereWrist(theta1_1, theta2_3, theta3_2, tf1, res_vec);
    //     SolveSphereWrist(theta1_1, theta2_4, theta3_2, tf1, res_vec);

    //     SolveSphereWrist(theta1_2, theta2_4, theta3_2, tf1, res_vec);
    //     SolveSphereWrist(theta1_2, theta2_3, theta3_2, tf1, res_vec);
    // }
    // else if (size_theta3 == 4)
    // {
    //     double theta3_1 = theta3_vec[0];
    //     double theta3_2 = theta3_vec[1];
    //     double theta3_3 = theta3_vec[2];
    //     double theta3_4 = theta3_vec[3];

    //     std::vector<double> theta2_vec;
    //     SolveTheta2(theta3_1, &theta2_vec);
    //     SolveTheta2(theta3_2, &theta2_vec);
    //     SolveTheta2(theta3_3, &theta2_vec);
    //     SolveTheta2(theta3_4, &theta2_vec);
    //     double theta2_1 = theta2_vec[0];
    //     double theta2_2 = theta2_vec[1];
    //     double theta2_3 = theta2_vec[2];
    //     double theta2_4 = theta2_vec[3];
    //     double theta2_5 = theta2_vec[4];
    //     double theta2_6 = theta2_vec[5];
    //     double theta2_7 = theta2_vec[6];
    //     double theta2_8 = theta2_vec[7];

    //     SolveSphereWrist(theta1_1, theta2_1, theta3_1, tf1, res_vec);
    //     SolveSphereWrist(theta1_1, theta2_2, theta3_1, tf1, res_vec);

    //     SolveSphereWrist(theta1_2, theta2_2, theta3_1, tf1, res_vec);
    //     SolveSphereWrist(theta1_2, theta2_1, theta3_1, tf1, res_vec);

    //     SolveSphereWrist(theta1_1, theta2_3, theta3_2, tf1, res_vec);
    //     SolveSphereWrist(theta1_1, theta2_4, theta3_2, tf1, res_vec);

    //     SolveSphereWrist(theta1_2, theta2_4, theta3_2, tf1, res_vec);
    //     SolveSphereWrist(theta1_2, theta2_3, theta3_2, tf1, res_vec);

    //     SolveSphereWrist(theta1_1, theta2_5, theta3_3, tf1, res_vec);
    //     SolveSphereWrist(theta1_1, theta2_6, theta3_3, tf1, res_vec);

    //     SolveSphereWrist(theta1_2, theta2_6, theta3_3, tf1, res_vec);
    //     SolveSphereWrist(theta1_2, theta2_5, theta3_3, tf1, res_vec);

    //     SolveSphereWrist(theta1_1, theta2_7, theta3_4, tf1, res_vec);
    //     SolveSphereWrist(theta1_1, theta2_8, theta3_4, tf1, res_vec);

    //     SolveSphereWrist(theta1_2, theta2_8, theta3_4, tf1, res_vec);
    //     SolveSphereWrist(theta1_2, theta2_7, theta3_4, tf1, res_vec);
    // }
}

void KineSolver::SolveTheta3(double x, double y, std::vector<double> *theta3)
{
    // solve theta3
    double a3 = a_.at(2);
    double d4 = sqrt(std::pow(d_.at(3), 2) + std::pow(a_.at(3), 2));
    if (d_.at(3) < 0)
        d4 = -d4;

    double s3 = (x * x + y * y + p_w[2] * p_w[2] - std::pow(a3, 2) - std::pow(d4, 2)) / (2 * a3 * d4);
    double c3 = sqrt(1 - std::pow(s3, 2));
    if (std::abs(s3) > 1.0 || std::abs(c3) > 1.0)
        return;
    double theta3_1 = atan2(s3, c3) - theta_.at(2) - belta_;
    double theta3_2 = atan2(s3, -c3) - theta_.at(2) - belta_;
    // std::cout << "theta3_1:" << theta3_1 * 180 / M_PI << std::endl;
    // std::cout << "theta3_2:" << theta3_2 * 180 / M_PI << std::endl;
    LimitJoints(theta3_1);
    LimitJoints(theta3_2);

    if (theta3->size() != 0)
    {
        if (std::abs(theta3->at(0) - theta3_1) > 1e-3)
            theta3->emplace_back(theta3_1);
        if (std::abs(theta3->at(1) - theta3_2) > 1e-3)
            theta3->emplace_back(theta3_2);
    }
    else
    {
        theta3->emplace_back(theta3_1);
        theta3->emplace_back(theta3_2);
    }
}

void KineSolver::SolveTheta2(double theta3, std::vector<double> *theta2_vec)
{
    double a3 = a_.at(2);
    double d4 = sqrt(std::pow(d_.at(3), 2) + std::pow(a_.at(3), 2));
    if (d_.at(3) < 0)
        d4 = -d4;
    std::vector<double> vec1 = SolveTrigonometricEquation(a3 + d4 * sin(theta3 + belta_), d4 * cos(theta3 + belta_), p_w[2]);
    double theta2_1 = vec1[0];
    double theta2_2 = vec1[1];
    // std::cout << "theta2_1:" << theta2_1 * 180 / M_PI << std::endl;
    // std::cout << "theta2_1:" << theta2_2 * 180 / M_PI << std::endl;
    LimitJoints(theta2_1);
    LimitJoints(theta2_2);

    theta2_vec->emplace_back(theta2_1);
    theta2_vec->emplace_back(theta2_2);
}

void KineSolver::SolveSphereWrist(double theta1, double theta2, double theta3, const Eigen::Isometry3d &tf, std::vector<Eigen::VectorXd> *res_vec)
{
    // std::cout << "theta:" << theta1 * 180 / M_PI << "  " << theta2 * 180 / M_PI << "  " << theta3 * 180 / M_PI << std::endl;

    Eigen::Matrix3d R_0_3, R;
    R_0_3 = (mDH(0, theta1) * mDH(1, theta2) * mDH(2, theta3)).linear();
    if (((mDH(0, theta1) * mDH(1, theta2) * mDH(2, theta3) * mDH(3, 0)).translation() - p_w).norm() > 1e-3)
    {
        return;
    }
    // std::cout << "1-3:" << (DH(0, theta1) * DH(1, theta2) * DH(2, theta3) * DH(3, 0)).translation().transpose() << std::endl;

    R = R_0_3.transpose() * tf.linear();

    double theta4_1 = std::atan2(R(2, 2), (R(0, 2)) - theta_.at(3));
    double theta5_1 = std::atan2((-sqrt(R(0, 2) * R(0, 2) + R(2, 2) * R(2, 2))), R(1, 2)) - theta_.at(4);
    double theta6_1 = std::atan2(R(1, 1), -R(1, 0)) - theta_.at(5);

    double theta4_2 = std::atan2(-R(2, 2), (-R(0, 2)) - theta_.at(3));
    double theta5_2 = std::atan2((sqrt(R(0, 2) * R(0, 2) + R(2, 2) * R(2, 2))), R(1, 2)) - theta_.at(4);
    double theta6_2 = std::atan2(-R(1, 1), R(1, 0)) - theta_.at(5);
    // std::cout << "theta4_1:" << theta4_1 * 180 / M_PI << std::endl;
    // std::cout << "theta5_1:" << theta5_1 * 180 / M_PI << std::endl;
    // std::cout << "theta6_1:" << theta6_1 * 180 / M_PI << std::endl;

    // std::cout << "theta4_2:" << theta4_2 * 180 / M_PI << std::endl;
    // std::cout << "theta5_2:" << theta5_2 * 180 / M_PI << std::endl;
    // std::cout << "theta6_2:" << theta6_2 * 180 / M_PI << std::endl;

    LimitJoints(theta4_1);
    LimitJoints(theta4_2);
    LimitJoints(theta5_1);
    LimitJoints(theta6_1);
    LimitJoints(theta5_2);
    LimitJoints(theta6_2);

    Eigen::VectorXd res1(6), res2(6);
    res1 << theta1, theta2, theta3, theta4_1, theta5_1, theta6_1;
    res2 << theta1, theta2, theta3, theta4_2, theta5_2, theta6_2;

    Eigen::Quaterniond q1(ComputeBaseEndfTf(res1, false).linear());
    Eigen::Quaterniond q2(ComputeBaseEndfTf(res2, false).linear());
    Eigen::Quaterniond q(tf.linear());
    if ((q1.coeffs() - q.coeffs()).norm() < 1e-5)
    {
        res_vec->emplace_back(res1);
        // std::cout << "res:" << res1.transpose() * 180 / M_PI << std::endl;
    }
    if ((q2.coeffs() - q.coeffs()).norm() < 1e-5)
    {
        res_vec->emplace_back(res2);
        // std::cout << "res:" << res2.transpose() * 180 / M_PI << std::endl;
    }
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
