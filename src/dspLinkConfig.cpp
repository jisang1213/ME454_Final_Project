#include "dspLinkConfig.h"

dspLink* dspLinkConfig::GetLink(int id) {
    return &link_vec[id];
}

dspJoint* dspLinkConfig::GetJoint(int id) {
    return &joint_vec[id];
}

int dspLinkConfig::GetNumJoint() {
    return joint_num;
}

int dspLinkConfig::GetNumLink() {
    return link_num;
}

bool dspLinkConfig::AddLink(int id, double mass, double rotinertia) {
    link_num++;
    dspLink link_now = dspLink();
    link_now.Set(mass, rotinertia);
    link_vec.push_back(link_now);
    return true;
}

bool dspLinkConfig::AddJointRevolute(int id, int body1_id, double body1x, double body1y, int body2_id, double body2x, double body2y) {
    joint_num++;
    dspJoint joint_now = dspJoint();
    joint_now.Set(0, body1_id, body2_id, GetLink(body1_id), GetLink(body2_id));
    joint_now.SetPointP(body1x, body1y, body2x, body2y);
    joint_vec.push_back(joint_now);
    return true;
}

bool dspLinkConfig::AddJointPrismatic(int id, int body1_id, double body1x, double body1y, int body2_id, double body2x, double body2y, double axisAngle) {
    joint_num++;
    dspJoint joint_now = dspJoint();
    joint_now.Set(1, body1_id, body2_id, GetLink(body1_id), GetLink(body2_id));
    joint_now.SetPointP(body1x, body1y, body2x, body2y);
    joint_now.SetPointQ(axisAngle);
    joint_vec.push_back(joint_now);
    return true;
}
