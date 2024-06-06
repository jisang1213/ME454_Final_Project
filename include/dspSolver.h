#pragma once
#include <Eigen/Dense>
#include "dspLinkConfig.h"

/// A class to calculate forward dynamics simulation for the current configuration

class dspSolver
{
public:
    dspSolver(dspLinkConfig& linkConfig);

private:
    /// pointer of dspLinkConfig structure
    dspLinkConfig& LinkConfig;
    /// mass matrix, inverse mass matrix
    Eigen::MatrixXd M, M_inv;
    /// jacobian matrix
    Eigen::MatrixXd J;
    /// first derivative of jacobian matrix
    Eigen::MatrixXd J_dot;
    /// matrix definition if you need
    Eigen::MatrixXd MatLeft;
    /// vector definition if you need
    Eigen::VectorXd VecRight;

    /// link information vector q, which is a n-dimensional vector (n = 3 * number of links)
    /// you must set its order like below
    ///    [x position of link, y position of link, angle of link, ... ]
    Eigen::VectorXd q;
    /// first derivative of vector q
    Eigen::VectorXd qdot;
    /// second derivative of vector q
    Eigen::VectorXd qddot;
    /// external force/torque vector
    Eigen::VectorXd F_ext;
    /// lambda vector
    Eigen::VectorXd lambda;
        
private:
    /// make and update mass matrix : M
    bool Make_M(void);
    /// make and update Jacobian matrix which has first-order partial derivatives of vector q's components : J
    bool Make_J(void);
    /// make and update time derivative of Jacobian matrix  : J_dot
    bool Make_J_dot(void);
    /// make and update a n-dimensional vector for external force/torque : F_ext
    /// it contains x direction of force, y direction of force, and torque
    bool Make_F_ext(void);

public:
    /// calculate error of your constraint definition
    bool CalculateConstraintError(Eigen::VectorXd&);
    /// calculate second derivative of vector q : qddot
    bool CalcLinAlg(void);
    /// load and update information what you need in class variables
    /// current vector q and its derivative should be loaded : q, q_dot
    bool UpdateCurrentInfo(void);
    /// update vector q and its derivative : q, q_dot
    /// also save them in class variables for next step
    bool UpdateNextInfo(double timestep);

    /// get q information of current state
    bool GetQ(Eigen::VectorXd&);
    /// get first derivative of q information of current state
    bool GetQDot(Eigen::VectorXd&);
    /// update q information into each dspLink class
    bool SetQ(void);
    /// update first derivative of q information into each dspLink class
    bool SetQDot(void);
};
