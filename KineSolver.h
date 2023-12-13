#include <array>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <vector>

#define M_PI 3.14159265358979323846   /* pi */
#define M_PI_2 1.57079632679489661923 /* pi/2 */
#define M_PI_4 0.78539816339744830962 /* pi/4 */
#define DEG_2_RAD M_PI / 180
#define RAD_2_DEG 180 / M_PI

enum class IKStatus : size_t
{
  kSuccess = 0,
  kSolveFail = 1,
  kExceedJointLimit = 2,
};

class KineSolver
{
public:
  KineSolver(const std::array<double, 7> &a, const std::array<double, 7> &alpha, const std::array<double, 7> &d, const std::array<double, 7> &theta);
  ~KineSolver();

  Eigen::Isometry3d ComputeFK(const Eigen::VectorXd &jnts, bool is_deg = false);
  Eigen::Isometry3d ComputeBaseEndfTf(const Eigen::VectorXd &jnts, bool is_deg = false);

  IKStatus ComputeIK(const Eigen::Isometry3d &tf, const Eigen::VectorXd &init_jnt, Eigen::VectorXd *jnts);
  void ComputeAllIK(const Eigen::Isometry3d &tf, std::vector<Eigen::VectorXd> *res_vec);

  void SetJntLimit(const std::array<double, 6> &lower, const std::array<double, 6> &upper);

private:
  Eigen::Isometry3d DH(std::size_t id, double theta, bool is_deg = false);
  Eigen::Isometry3d mDH(std::size_t id, double theta, bool is_deg = false);
  void SolveSphereWrist(double theta1, double theta2, double theta3, const Eigen::Isometry3d &tf, std::vector<Eigen::VectorXd> *res_vec);
  IKStatus IsJointWithinLimit(const Eigen::VectorXd &jnt);
  void LimitJoints(double &joint);
  std::vector<double> SolveTrigonometricEquation(double a, double b, double c);
  void SolveTheta3(double x, double y,std::vector<double> *theta3);
  void SolveTheta2(double theta3,std::vector<double> *theta2);

  std::array<double, 6> jnt_upper_limit_;
  std::array<double, 6> jnt_lower_limit_;

  std::size_t dof_ = 6;
  std::array<double, 7> a_;
  std::array<double, 7> alpha_;
  std::array<double, 7> d_;
  std::array<double, 7> theta_;
  Eigen::Vector3d p_w{0, 0, 0};
  double belta_{0.0};
};