// this MUST be first, otherwise there might be problems on windows
// see: https://stackoverflow.com/questions/6563810/m-pi-works-with-math-h-but-not-with-cmath-in-visual-studio/6563891#6563891
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <random>
#include <chrono>

#include <kalman/ExtendedKalmanFilter.hpp>
#include <kalman/Types.hpp>
#include "SystemModel.hpp"

typedef float T;

// Some type shortcuts
typedef GenericQuaternionImuModel::State<T> State;
typedef GenericQuaternionImuModel::NoiseState<T> NoiseState;
typedef GenericQuaternionImuModel::Control<T> Control;
typedef GenericQuaternionImuModel::SystemModel<T> SystemModel;

//typedef GenericQuaternionImuModel::PositionMeasurement<T> PositionMeasurement;
//typedef GenericQuaternionImuModel::OrientationMeasurement<T> OrientationMeasurement;
//typedef GenericQuaternionImuModel::PositionMeasurementModel<T> PositionModel;
//typedef GenericQuaternionImuModel::OrientationMeasurementModel<T> OrientationModel;

int main(int argc, char** argv)
{
    // Simulated (true) system state
    State x;
    x.get().setZero();
    x.qw() = 1.0;
    
    // Control input
    Control u;
    u.ax() = 0.0;
    u.ay() = 0.0;
    u.az() = 9.81;
    u.wx() = 0.0;
    u.wy() = 0.0;
    u.wz() = 0.0;    

    NoiseState noise_state;
    noise_state.accNoise() = Vector<T, 3>(0.35, 0.35, 0.35);
    noise_state.gyroNoise() = Vector<T, 3>(0.01, 0.01, 0.01);
    noise_state.accDrift() = Vector<T, 3>(1e-3, 1e-3, 1e-3);
    noise_state.gyroDrift() = Vector<T, 3>(1e-3, 1e-3, 1e-3);

    // System
    SystemModel sys;
    std::cout << "x:" << x.get().transpose() << std::endl;
    sys.setupProcessCov(x, noise_state);
    std::cout << "syscov:\n" << sys.getCovariance() << std::endl;
    // Measurement models
    // Set position landmarks at (-10, -10) and (30, 75)
    //PositionModel pm(-10, -10, 30, 75);
    //OrientationModel om;
    
    // Random number generation (for noise simulation)
    std::default_random_engine generator;
    generator.seed( std::chrono::system_clock::now().time_since_epoch().count() );
    std::normal_distribution<T> noise(0, 1);
    
    // Some filters for estimation
    // Pure predictor without measurement updates
    Kalman::ExtendedKalmanFilter<State> predictor;
    // Extended Kalman Filter
    //Kalman::ExtendedKalmanFilter<State> ekf;
    // Unscented Kalman Filter
    //Kalman::UnscentedKalmanFilter<State> ukf(1);
    
    sys.setDiffTime(0.01);

    // Init filters with true system state
    predictor.init(x);

    u.a() = Vector<T, 3>(0.0, 0.0, -9.80655);
    u.w() = Vector<T, 3>(0.0, 0.0, 0.0);

    u.ax() += noise_state.ax_n()*noise(generator);
    u.ay() += noise_state.ay_n()*noise(generator);
    u.az() += noise_state.az_n()*noise(generator);
    u.wx() += noise_state.wx_n()*noise(generator);
    u.wy() += noise_state.wy_n()*noise(generator);
    u.wz() += noise_state.wz_n()*noise(generator);

    predictor.predict(sys, u);
    //ekf.init(x);
    //ukf.init(x);
    
    /*// Standard-Deviation of noise added to all state vector components during state transition
    T systemNoise = 0.1;
    // Standard-Deviation of noise added to all measurement vector components in orientation measurements
    //T orientationNoise = 0.025;
    // Standard-Deviation of noise added to all measurement vector components in distance measurements
    //T distanceNoise = 0.25;
    
    // Simulate for 100 steps
    const size_t N = 100;
    for(size_t i = 1; i <= N; i++)
    {
        // Generate some control input
        u.v() = 1. + std::sin( T(2) * T(M_PI) / T(N) );
        u.dtheta() = std::sin( T(2) * T(M_PI) / T(N) ) * (1 - 2*(i > 50));
        
        // Simulate system
        x = sys.f(x, u);
        
        // Add noise: Our robot move is affected by noise (due to actuator failures)
        x.x() += systemNoise*noise(generator);
        x.y() += systemNoise*noise(generator);
        x.theta() += systemNoise*noise(generator);
        
        // Predict state for current time-step using the filters
        auto x_pred = predictor.predict(sys, u);
        auto x_ekf = ekf.predict(sys, u);
        auto x_ukf = ukf.predict(sys, u);
        
        // Orientation measurement
        {
            // We can measure the orientation every 5th step
            OrientationMeasurement orientation = om.h(x);
            
            // Measurement is affected by noise as well
            orientation.theta() += orientationNoise * noise(generator);
            
            // Update EKF
            x_ekf = ekf.update(om, orientation);
            
            // Update UKF
            x_ukf = ukf.update(om, orientation);
        }
        
        // Position measurement
        {
            // We can measure the position every 10th step
            PositionMeasurement position = pm.h(x);
            
            // Measurement is affected by noise as well
            position.d1() += distanceNoise * noise(generator);
            position.d2() += distanceNoise * noise(generator);
            
            // Update EKF
            x_ekf = ekf.update(pm, position);
            
            // Update UKF
            x_ukf = ukf.update(pm, position);
        }
        
        // Print to stdout as csv format
        std::cout   << x.x() << "," << x.y() << "," << x.theta() << ","
                    << x_pred.x() << "," << x_pred.y() << "," << x_pred.theta()  << ","
                    << x_ekf.x() << "," << x_ekf.y() << "," << x_ekf.theta()  << ","
                    << x_ukf.x() << "," << x_ukf.y() << "," << x_ukf.theta()
                    << std::endl;
    }
    */
    return 0;
}
