#ifndef MATHS_UTILS_HPP_
#define MATHS_UTILS_HPP_

#include <Eigen/Dense>

template <typename T, size_t N>
using Vector = Eigen::Matrix<T, N, 1>;
template <typename T>
using Quaternion = Eigen::Quaternion<T>;

namespace MathsUtils
{
  template <typename T>
  Eigen::Matrix<T, 3, 3> skew_symmetric(Vector<T, 3> v)
  {
    Eigen::Matrix<T, 3, 3> m;
    m <<    0, -v[2],  v[1],
         v[2],     0, -v[0],
        -v[1],  v[0],     0;
    return m;
  }

  // Differentiation of qvqstar in http://web.cs.iastate.edu/~cs577/handouts/quaternion.pdf
  // qvqstar = (qw^2 - norm(qv))v + 2(qv.v)qv + 2qw(qv x v)
  // d(qvqstar) / dqw = 2 qw v + 2 (qv x v) // col - 0
  // d(qvqstar) / dqv = 2 (-v qv' + v.v I + qvv' - qw [v]_{skew}) // col - 1,2,3
  template <typename T>
  Eigen::Matrix<T, 3, 4> diff_qvqstar_q(Quaternion<T> q, Vector<T, 3> v)
  {
    auto& qw = q.w();
    Vector<T, 3> qv = q.vec();
    Eigen::Matrix<T, 3, 4> D(3, 4);
    D.col(0) = 2 * (qw * v + skew_symmetric(qv) * v);
    D.template block<3, 3>(0, 1) = 
      2 * (-v * qv.transpose() + 
           v.dot(qv)*Eigen::Matrix<T, 3, 3>::Identity() + 
           qv*v.transpose() - 
           qw*skew_symmetric(v));
    return D;
  }

  /**
   * Returns 
   *  [
   *    [qw*qw + qx*qx - qz*qz - qy*qy, -qz*qw + qy*qx - qw*qz + qx*qy,	qy*qw + qz*qx + qx*qz + qw*qy
   *     qx*qy + qw*qz + qz*qw + qy*qx,	 qy*qy - qz*qz + qw*qw - qx*qx,	qz*qy + qy*qz - qx*qw - qw*qx
   *     qx*qz - qw*qy + qz*qx - qy*qw,	 qy*qz + qz*qy + qw*qx + qx*qw,	qz*qz - qy*qy - qx*qx + qw*qw]
   *  ]
   */ 
  template <typename T>
  Eigen::Matrix<T, 3, 3> diff_qvqstar_v(Quaternion<T> q)
  {
    auto& qw = q.w();
    Vector<T, 3> qv = q.vec();
    Eigen::Matrix<T, 3, 3> D;
    D = (qw*qw - qv.dot(qv))*Eigen::Matrix<T, 3, 3>::Identity() + 2*qv*qv.transpose() + 2*qw*skew_symmetric(qv);
    return D; 
  }

  template <typename T>
  Eigen::Matrix<T, 4, 4> diff_pq_p(Quaternion<T> q)
  {
    auto& qw = q.w();
    Vector<T, 3> qv = q.vec();
    Eigen::Matrix<T, 4, 4> D;
    D(0, 0) = qw;
    D.template block<1, 3>(0, 1) = -qv.transpose();
    D.template block<3, 1>(1, 0) = qv;
    D.template block<3, 3>(1, 1) = Eigen::Matrix<T, 3, 3>::Identity()*qw - skew_symmetric(qv);
    return D;
  }

  template <typename T>
  //diff_(p*q) /diff_q
  Eigen::Matrix<T, 4, 4> diff_pq_q(Quaternion<T> p)
  {
    auto& pw = p.w();
    Vector<T, 3> pv = p.vec();
    Eigen::Matrix<T, 4, 4> D;
    D(0, 0) = pw;
    D.template block<1, 3>(0, 1) = -pv.transpose();
    D.template block<3, 1>(1, 0) = pv;
    D.template block<3, 3>(1, 1) = Eigen::Matrix<T, 3, 3>::Identity()*pw + skew_symmetric(pv);
    return D;
  }
}
#endif