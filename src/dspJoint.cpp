#include "dspJoint.h"

dspJoint::dspJoint() {

}

bool dspJoint::Set(int type0, int body1_id, int body2_id, dspLink *pointer1, dspLink *pointer2 ) {
    type = type0;
    link1_id = body1_id;
    link2_id = body2_id;
    link1_pointer = pointer1;
    link2_pointer = pointer2;
    return true;
}

bool dspJoint::SetPointQ(double axisAngle) {
    double link1_theta = link1_pointer->GetTheta();
    double link2_theta = link2_pointer->GetTheta();
    double globalAxisX = 1.0;
    double globalAxisY = 0.0;
    link1_sx2 = link1_sx + globalAxisX *cos(axisAngle - link1_theta) - globalAxisY * sin(axisAngle - link1_theta);
    link1_sy2 = link1_sy + globalAxisX *sin(axisAngle - link1_theta) + globalAxisY * cos(axisAngle - link1_theta);
    link2_sx2 = link2_sx + globalAxisX *cos(axisAngle - link2_theta) - globalAxisY * sin(axisAngle - link2_theta);
    link2_sy2 = link2_sy + globalAxisX *sin(axisAngle - link2_theta) + globalAxisY * cos(axisAngle - link2_theta);
    return true;
}

bool dspJoint::SetPointP(double body1_px, double body1_py, double body2_px, double body2_py) {
    link1_sx = body1_px;
    link1_sy = body1_py;
    link2_sx = body2_px;
    link2_sy = body2_py;
    return true;
}

void dspJoint::GetP1_Link1_LCS(Eigen::Vector2d& vec) {
    vec << link1_sx, link1_sy;
}

void dspJoint::GetP2_Link2_LCS(Eigen::Vector2d& vec) {
    vec << link2_sx, link2_sy;
}

void dspJoint::GetQ1_Link1_LCS(Eigen::Vector2d& vec) {
    vec << link1_sx2, link1_sy2;
}

void dspJoint::GetQ2_Link2_LCS(Eigen::Vector2d& vec) {
    vec << link2_sx2, link2_sy2;
}

int dspJoint::GetType() {
    return type;
}

int dspJoint::GetLink1ID() {
    return link1_id;
}

int dspJoint::GetLink2ID() {
    return link2_id;
}