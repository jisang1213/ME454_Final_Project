#pragma once
#include "dspLink.h"
#include "dspJoint.h"
#include <vector>

/// A class to store and process configuration of links and joints

class dspLinkConfig
{
private:
    /// make container of dspLink in order
    std::vector<dspLink> link_vec;
    /// make container of dspJoint in order
    std::vector<dspJoint> joint_vec;
    /// number of link in system
    int link_num = 0;
    /// number of joint in system
    int joint_num = 0;

public:
    /// get pointer of dspLink, which has id
    dspLink* GetLink(int id);
    /// get pointer of dspJoint, which has id
    dspJoint* GetJoint(int id);
    /// get number of links in system
    int GetNumLink();
    /// get number of joints in system
    int GetNumJoint();
    /// get information of new link and update into link_vec
    bool AddLink(int id, double mass, double rotInertia);
    /// get information of new joint and update into joint_vec
    /// input should be positions of joint in each link's local coordinate (revolute)
    bool AddJointRevolute(int id, int body1_id, double body1x, double body1y, int body2_id, double body2x, double body2y);
    /// input should be positions of joint in each link's local coordinate and angle of translational axis (prismatic)
    bool AddJointPrismatic(int id, int body1_id, double body1x, double body1y, int body2_id, double body2x, double body2y, double axisAngle);
};
