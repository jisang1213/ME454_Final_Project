#pragma once

/// A class to store and process link information

class dspLink
{
private:
    /// mass of link
    double Mass = 0;
    /// rotational inertia of link
    double RotInertia = 0;
    /// position of x and y in world coordinate
    double px = 0, py = 0;
    /// velocity of x and y in world coordinate
    double vx = 0, vy = 0;
    /// acceleration of x and y in world coordinate
    double ax = 0, ay = 0;
    /// angle, angular velocity, angular acceleration in world coordinate
    double theta = 0, omega = 0, alpha = 0;

public:
    dspLink();
    /// get mass of link
    double GetMass();
    /// get rotational inertia of link
    double GetInertia();
    /// get angle of link
    double GetTheta();
    /// get angular velocity of link
    double GetOmega();
    /// get angular acceleration of link
    double GetAlpha();
    /// set mass and rotational inertia of link
    bool Set(double mass, double rotInertia);
    /// set position and angle of link
    bool SetQ(double px, double py, double theta);
    /// set velocity and angular velocity of link
    bool SetQDot(double vx, double vy, double omega);
    /// set position and angular acceleration of link
    bool SetQDDot(double ax, double ay, double alpha);
    /// get position and angle of link
    bool GetQ(double &px, double &py, double &theta);
    /// get velocity and angular velocity of link
    bool GetQDot(double &vx, double &vy, double &omega);
    /// get acceleration and angular acceleration of link
    bool GetQDDot(double &ax, double &ay, double &alpha);
};
