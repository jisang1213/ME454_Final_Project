#pragma once
#include "dspLink.h"
#include <iostream>
#include <Eigen/Dense>

/// A class to store and process joint information
class dspJoint
{
private:
    /// joint type [0 : revolute, 1: prismatic]
    int type = 0;
    /// position of P1 in local coordinate link1
    double link1_sx = 0.0;
    double link1_sy = 0.0;
    /// position of P2 in local coordinate link2
    double link2_sx = 0.0;
    double link2_sy = 0.0;
    /// position of Q1 in local coordinate link1
    double link1_sx2 = 0.0;
    double link1_sy2 = 0.0;
    /// position of Q2 in local coordinate link2
    double link2_sx2 = 0.0;
    double link2_sy2 = 0.0;
    /// angle of translational axis in global coordinate
    double axisAngle = 0.0;
    /// id of connected links
    int link1_id = 0;
    int link2_id = 0;
    /// pointer of connected links
    dspLink* link1_pointer = nullptr;
    dspLink* link2_pointer = nullptr;

public:
    dspJoint();
    /// set joint type, id of parent_link, id of child_link and pointers of each link
    bool Set(int type, int link1_id, int link2_id, dspLink* pointer1, dspLink* pointer2);
    /// set P points in each link's local coordinate
    bool SetPointP(double link1_sx, double link1_sy, double link2_sx, double link2_sy);
    /// set Q points in each link's local coordinate
    /// input : angle of translational axis in global coordinate
    bool SetPointQ(double axisAngle);
    /// get position of P1 in local coordinate of link1
    void GetP1_Link1_LCS(Eigen::Vector2d&);
    /// get position of P2 in local coordinate of link2
    void GetP2_Link2_LCS(Eigen::Vector2d&);
    /// get position of Q1 in local coordinate of link1
    void GetQ1_Link1_LCS(Eigen::Vector2d&);
    /// get position of Q2 in local coordinate of link2
    void GetQ2_Link2_LCS(Eigen::Vector2d&);
    /// get type of joint
    int GetType();
    /// get id of link1
    int GetLink1ID();
    /// get id of link2
    int GetLink2ID();
};
