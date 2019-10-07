#ifndef GENERIC_QUATERNION_IMU_MODEL_SYSTEMMODEL_HPP_
#define GENERIC_QUATERNION_IMU_MODEL_SYSTEMMODEL_HPP_

#include <kalman/LinearizedSystemModel.hpp>
#include "MathsUtils.hpp"

namespace GenericQuaternionImuModel
{

/**
 * @brief System state vector-type for a generic quaternion based system model
 *   of an inertial measurement unit (IMU).
 *
 * This is a system state for an inertial measurment unit (IMU) characterized by 
 * its position, velocity, orientation, and biases in its acceleration and 
 * angular velocity measurements.
 *
 * @param T Numeric scalar type
 */
template<typename T>
class State
{
  KALMAN_VECTOR(State, T, 16)

public:
  DEFINE_CONST_GET(T, px, PX);
  DEFINE_CONST_GET(T, py, PY);
  DEFINE_CONST_GET(T, pz, PZ);
  DEFINE_CONST_GET(T, vx, VX);
  DEFINE_CONST_GET(T, vy, VY);
  DEFINE_CONST_GET(T, vz, VZ);
  DEFINE_CONST_GET(T, qw, QW);
  DEFINE_CONST_GET(T, qx, QX);
  DEFINE_CONST_GET(T, qy, QY);
  DEFINE_CONST_GET(T, qz, QZ);
  DEFINE_CONST_GET(T, bax, BAX);
  DEFINE_CONST_GET(T, bay, BAY);
  DEFINE_CONST_GET(T, baz, BAZ);
  DEFINE_CONST_GET(T, bwx, BWX);
  DEFINE_CONST_GET(T, bwy, BWY);
  DEFINE_CONST_GET(T, bwz, BWZ);

  DEFINE_REF_GET(T, px, PX);
  DEFINE_REF_GET(T, py, PY);
  DEFINE_REF_GET(T, pz, PZ);
  DEFINE_REF_GET(T, vx, VX);
  DEFINE_REF_GET(T, vy, VY);
  DEFINE_REF_GET(T, vz, VZ);
  DEFINE_REF_GET(T, qw, QW);
  DEFINE_REF_GET(T, qx, QX);
  DEFINE_REF_GET(T, qy, QY);
  DEFINE_REF_GET(T, qz, QZ);
  DEFINE_REF_GET(T, bax, BAX);
  DEFINE_REF_GET(T, bay, BAY);
  DEFINE_REF_GET(T, baz, BAZ);
  DEFINE_REF_GET(T, bwx, BWX);
  DEFINE_REF_GET(T, bwy, BWY);
  DEFINE_REF_GET(T, bwz, BWZ);

  const typename ConstFixedSegmentReturnType<3>::Type getPosition() const { return this->internal_.template segment<3>(PX); }
  typename FixedSegmentReturnType<3>::Type position() { return this->internal_.template segment<3>(PX); }
  const typename ConstFixedSegmentReturnType<3>::Type getVelocity() const { return this->internal_.template segment<3>(VX); }
  typename FixedSegmentReturnType<3>::Type velocity() { return this->internal_.template segment<3>(VX); }
  const typename ConstFixedSegmentReturnType<3>::Type getBa() const { return this->internal_.template segment<3>(BAX); }
  typename FixedSegmentReturnType<3>::Type ba() { return this->internal_.template segment<3>(BAX); }
  const typename ConstFixedSegmentReturnType<3>::Type getBw() const { return this->internal_.template segment<3>(BWX); }
  typename FixedSegmentReturnType<3>::Type bw() { return this->internal_.template segment<3>(BWX); }

  const Quaternion<T> getOrientation() const {
    Quaternion<T> q;
    q.w() = qw(); q.vec() = this->internal_.template segment<3>(QX);
    return q;
  }

  void setOrientation(const Quaternion<T>& q) {
    qw() = q.w(); this->internal_.template segment(QX, 3) = q.vec();
  }

  void normalizeQ() {
    this->internal_.template segment<4>(QW).normalize();
  }

  enum StateIndices {
    PX, PY, PZ,
    VX, VY, VZ,
    QW, QX, QY, QZ,
    BAX, BAY, BAZ,
    BWX, BWY, BWZ
  };
};

template<typename T>
class NoiseState
{
  KALMAN_VECTOR(NoiseState, T, 12)

public:
  DEFINE_CONST_GET(T, ax_n, AX_N);
  DEFINE_CONST_GET(T, ay_n, AY_N);
  DEFINE_CONST_GET(T, az_n, AZ_N);
  DEFINE_CONST_GET(T, wx_n, WX_N);
  DEFINE_CONST_GET(T, wy_n, WY_N);
  DEFINE_CONST_GET(T, wz_n, WZ_N);
  DEFINE_CONST_GET(T, ax_d, AX_D);
  DEFINE_CONST_GET(T, ay_d, AY_D);
  DEFINE_CONST_GET(T, az_d, AZ_D);
  DEFINE_CONST_GET(T, wx_d, WX_D);
  DEFINE_CONST_GET(T, wy_d, WY_D);
  DEFINE_CONST_GET(T, wz_d, WZ_D);
  
  DEFINE_REF_GET(T, ax_n, AX_N);
  DEFINE_REF_GET(T, ay_n, AY_N);
  DEFINE_REF_GET(T, az_n, AZ_N);
  DEFINE_REF_GET(T, wx_n, WX_N);
  DEFINE_REF_GET(T, wy_n, WY_N);
  DEFINE_REF_GET(T, wz_n, WZ_N);
  DEFINE_REF_GET(T, ax_d, AX_D);
  DEFINE_REF_GET(T, ay_d, AY_D);
  DEFINE_REF_GET(T, az_d, AZ_D);
  DEFINE_REF_GET(T, wx_d, WX_D);
  DEFINE_REF_GET(T, wy_d, WY_D);
  DEFINE_REF_GET(T, wz_d, WZ_D);
  
  typename FixedSegmentReturnType<3>::Type accNoise() 
    { return this->internal_.template segment<3>(AX_N); }
  typename FixedSegmentReturnType<3>::Type gyroNoise() 
    { return this->internal_.template segment<3>(WX_N); }
  typename FixedSegmentReturnType<3>::Type accDrift() 
    { return this->internal_.template segment<3>(AX_D); }
  typename FixedSegmentReturnType<3>::Type gyroDrift() 
    { return this->internal_.template segment<3>(WX_D); }

  Eigen::Matrix<T, 3, 3> accNoiseCov() 
    { return accNoise().array().square().matrix().asDiagonal(); }
  Eigen::Matrix<T, 3, 3> gyroNoiseCov() 
    { return gyroNoise().array().square().matrix().asDiagonal(); }
  Eigen::Matrix<T, 3, 3> accDriftCov() 
    { return accDrift().array().square().matrix().asDiagonal(); }
  Eigen::Matrix<T, 3, 3> gyroDriftCov() 
    { return gyroDrift().array().square().matrix().asDiagonal(); }

  enum NoiseIndices {
    AX_N, AY_N, AZ_N,
    WX_N, WY_N, WZ_N,
    AX_D, AY_D, AZ_D,
    WX_D, WY_D, WZ_D
  };
};

/**
 * @brief System control-input vector-type for an inertial measurement unit (IMU).
 *
 * This is a system control-input of an inertial measurment unit (IMU) given in terms
 * of its acceleration and angular velocity measurements.
 *
 * @param T Numeric scalar type
 */
template<typename T>
class Control
{
  KALMAN_VECTOR(Control, T, 6)
  
public:
  DEFINE_CONST_GET(T, ax, AX);
  DEFINE_CONST_GET(T, ay, AY);
  DEFINE_CONST_GET(T, az, AZ);
  DEFINE_CONST_GET(T, wx, WX);
  DEFINE_CONST_GET(T, wy, WY);
  DEFINE_CONST_GET(T, wz, WZ);

  DEFINE_REF_GET(T, ax, AX);
  DEFINE_REF_GET(T, ay, AY);
  DEFINE_REF_GET(T, az, AZ);
  DEFINE_REF_GET(T, wx, WX);
  DEFINE_REF_GET(T, wy, WY);
  DEFINE_REF_GET(T, wz, WZ);

  typename FixedSegmentReturnType<3>::Type a() { return this->internal_.template segment<3>(AX); }
  typename FixedSegmentReturnType<3>::Type w() { return this->internal_.template segment<3>(WX); }
  const typename ConstFixedSegmentReturnType<3>::Type getA() const { return this->internal_.template segment<3>(AX); }
  const typename ConstFixedSegmentReturnType<3>::Type getW() const { return this->internal_.template segment<3>(WX); }

  enum InputIndices {
    AX, AY, AZ,
    WX, WY, WZ
  };
};

/**
 * @brief System model for update of an inertial measurement unit (IMU).
 *
 * This is the system model defining how the state of the imu evolves over time.
 *
 * @param T Numeric scalar type
 */
template<typename T, template<class> class CovarianceBase = Kalman::StandardBase>
class SystemModel : public Kalman::LinearizedSystemModel<State<T>, Control<T>, CovarianceBase>
{
public:
  //! State type shortcut definition
  typedef GenericQuaternionImuModel::State<T> S;
  
  //! Control type shortcut definition
  typedef GenericQuaternionImuModel::Control<T> C;
  
  //! Noise state type shortcut definition
  typedef GenericQuaternionImuModel::NoiseState<T> NS;

  //! Covariance type shortcut definition
  typedef Kalman::Covariance<State<T>> Covariance;

  void setupProcessCov(const S& x, const NS& noise_state) {
    this->noise_state_ = noise_state;
    process_cov_base_.setZero();
    std::cout << "noise_state:" << this->noise_state_.get().transpose() << std::endl;

    // variance of "change in position" depends on acceleration measurement by dt times
    process_cov_base_.template block<3, 3>(S::StateIndices::PX, S::StateIndices::PX) = noise_state_.accNoiseCov();

    // variance of "change in velocty" depends on acceleration measurement by 1.0 times
    // Note that velocity and position depends on rotated acceleration however
    // rotation matrices resulting from state orientation end up as R.T * R in the filter
    // resulting in an identity so we do not need to process them here
    process_cov_base_.template block<3, 3>(S::StateIndices::VX, S::StateIndices::VX) = noise_state_.accNoiseCov();

    cache_.diff_pq_q = 
      Eigen::Matrix<T, 4, 3>(MathsUtils::diff_pq_q<T>(x.getOrientation()).template block<4, 3>(0, 1));

    // variance of "change in orientation" depends on rotated gyro measurements
    process_cov_base_.template block<4, 4>(S::StateIndices::QW, S::StateIndices::QW) = 
      cache_.diff_pq_q * 0.5 * noise_state_.gyroNoiseCov() * cache_.diff_pq_q.transpose();

    // variance of "change in acceleration bias" depends on acceleration drift-dot by 1.0 times
    process_cov_base_.template block<3, 3>(S::StateIndices::BAX, S::StateIndices::BAX) = noise_state_.accDriftCov();

    // variance of "change in gyro bias" depends on gyro drift-dot by 1.0 times
    process_cov_base_.template block<3, 3>(S::StateIndices::BWX, S::StateIndices::BWX) = noise_state_.gyroDriftCov();
    this->setCovariance(process_cov_base_);
  }

  /**
   * @brief Definition of (non-linear) state transition function
   *
   * This function defines how the system state is propagated through time,
   * i.e. it defines in which state \f$\hat{x}_{k+1}\f$ is system is expected to 
   * be in time-step \f$k+1\f$ given the current state \f$x_k\f$ in step \f$k\f$ and
   * the system control input \f$u\f$.
   *
   * @param [in] x The system state in current time-step
   * @param [in] u The control vector input
   * @returns The (predicted) system state in the next time-step
   */
  
  S f(const S& x, const C& u) const
  {
      //! Predicted state vector after transition
      S x_dot;
      x_dot.get().setZero();

      // update position
      x_dot.position() = x.getVelocity(); // position differentiation = velocity

      // update velocity
      Quaternion<T> q = x.getOrientation();
      std::cout << "qw:" << q.w() << std::endl;
      std::cout << "qv:" << q.vec() << std::endl;
      std::cout << "q1w:" << cache_.acc_in_body_q.w() << std::endl;
      std::cout << "q1v:" << cache_.acc_in_body_q.vec() << std::endl;
      Quaternion<T> acc_in_inertial =  q * cache_.acc_in_body_q * q.inverse(); // Rotate acceleration from body frame to inertial frame
      x_dot.velocity() = acc_in_inertial.vec() - gravity_; // Remove gravitational acceleration and set as velocity differentiation
      std::cout << "acc_in_inertial.vec():" << acc_in_inertial.vec().transpose() << std::endl;
      std::cout << "x_dot.velocity():" << x_dot.velocity().transpose() << std::endl;
      // update orientation
      Quaternion<T> q_dot = q * cache_.gyro_q; // Change in orientation = previous orientation multiplied by half of angular velocity
      q_dot.w() /= 2.0; q_dot.vec() /= 2.0;
      x_dot.setOrientation(q_dot); // Set differential of velocity

      std::cout << "x_dot:" << x_dot.get().transpose() << std::endl;

      // biases don't update...
      S x_new(x_dot.get() * this->dt_); // update state with difference of time
      x_new.normalizeQ();
      return x_new;
  }

  void setDiffTime(const double& dt) { dt_ = dt; }
  
protected:
  void updateJacobians( const S& x, const C& u )
  {
    // F = df/dx (Jacobian of state transition w.r.t. the state)
    this->F.setIdentity();
    this->W.setIdentity();

    // Remove bias from input acceleration measurement
    cache_.acc_in_body_q = Quaternion<T>(0, 0, 0, 0);
    cache_.acc_in_body_q.vec() = u.getA() - x.getBa();

    // Remove gyro bias from input gyro measurement
    cache_.gyro_q = Quaternion<T>(0, 0, 0, 0);
    cache_.gyro_q.vec() = u.getW() - x.getBw();

    Quaternion<T> q = x.getOrientation();
    // diff position term by velocity state variables
    this->F.template block<3, 3>(0, 3) = Eigen::Matrix<T, 3, 3>::Identity(); 
    // diff velocity term by orientation state variables
    this->F.template block<3, 4>(3, 6) = MathsUtils::diff_qvqstar_q<T>(q, cache_.acc_in_body_q.vec()); 
    // diff velocity term by acceleration bias state variables
    this->F.template block<3, 3>(3, 10) = -MathsUtils::diff_qvqstar_v<T>(q); 
    // diff orientation term by orientation state variables
    this->F.template block<4, 4>(6, 6) = 0.5*MathsUtils::diff_pq_p<T>(cache_.gyro_q); 
    // diff orientation term by orientation state variables
    cache_.diff_pq_q = 
      Eigen::Matrix<T, 4, 3>(MathsUtils::diff_pq_q<T>(q).template block<4, 3>(0, 1));
    this->F.template block<4, 3>(6, 13) = -0.5*MathsUtils::diff_pq_q<T>(q).template block<4, 3>(0, 1);

    updateProcessCov(x);
    this->F += this->F * this->dt_;
    this->W *= this->dt_;
  }

  void updateProcessCov(const S& x) {
    Covariance process_cov = process_cov_base_;
    process_cov.template block<3, 3>(S::StateIndices::PX, S::StateIndices::PX) *= this->dt_;

    process_cov.template block<4, 4>(S::StateIndices::QW, S::StateIndices::QW) = 
      cache_.diff_pq_q * 0.5 * this->noise_state_.gyroNoiseCov() * cache_.diff_pq_q.transpose();
    this->setCovariance(process_cov);
  }

private:
  struct Cache {
    Quaternion<T> acc_in_body_q = {Quaternion<T>(0, 0, 0, 0)};
    Quaternion<T> gyro_q = {Quaternion<T>(0, 0, 0, 0)};
    Eigen::Matrix<T, 4, 3> diff_pq_q;
  } cache_;

  bool system_initialized_ = {false};
  double dt_;
  const Vector<T, 3> gravity_ = Vector<T, 3>(0, 0, 9.80655);
  Covariance process_cov_base_;
  NS noise_state_;
};

} // namespace GenericQuaternionImuModel
/*if(!system_initialized_)
{
  system_initialized_ = true; 
  auto phy = atan2(u.ay(), u.az());
  auto theta = atan2(-u.ax(), u.az());
  Vector<T, 3> rpy(phy, theta, 0);
  Quaternion<T> q = euler2quaternion(rpy);
  x.qw() = q.w();
  x.segment<3>(QX) = q.vec();
  return;
} part of main.cpp*/
#endif