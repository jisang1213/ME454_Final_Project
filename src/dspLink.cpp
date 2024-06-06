#include "dspLink.h"

dspLink::dspLink() {
}

double dspLink::GetMass() {
    return Mass;
}

double dspLink::GetInertia() {
    return RotInertia;
}

double dspLink::GetTheta() {
    return theta;
}

double dspLink::GetOmega() {
    return omega;
}

double dspLink::GetAlpha() {
    return alpha;
}

bool dspLink::Set(double mass, double rotInertia) {
    Mass = mass;
    RotInertia = rotInertia;
    return true;
}

bool dspLink::SetQ(double px0, double py0, double theta0) {
    px = px0;
    py = py0;
    theta = theta0;
    return true;
}

bool dspLink::SetQDot(double vx0, double vy0, double omega0) {
    vx = vx0;
    vy = vy0;
    omega = omega0;
    return true;
}

bool dspLink::SetQDDot(double ax, double ay, double alpha) {
    return true;
}

bool dspLink::GetQ(double& px0, double& py0, double& theta0) {
    px0 = px;
    py0 = py;
    theta0 = theta;
    return true;
}

bool dspLink::GetQDot(double& vx0, double& vy0, double& omega0) {
    vx0 = vx;
    vy0 = vy;
    omega0 = omega;
    return true;
}

bool dspLink::GetQDDot(double& ax, double& ay, double& alpha) {
    return true;
}