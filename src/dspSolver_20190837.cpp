#include "dspSolver.h"
#include <iostream>

// Implement this function
dspSolver::dspSolver(dspLinkConfig& linkConfig) : LinkConfig(linkConfig) {
    /// Solver variable initialization
    /// The link with index 0 is the ground link
    /// There are only two joint types : revolute and prismatic
    int link_num = LinkConfig.GetNumLink();
    int joint_num = LinkConfig.GetNumJoint();

    /// You have to initialize following variables
    int size = 3*(link_num-1);
    q.resize(size);
    qdot.resize(size);
    qddot.resize(size);
    M.resize(size, size);
    M_inv.resize(size, size);
    J.resize(2*joint_num, size);
    J_dot.resize(2*joint_num, size);
    F_ext.resize(size);
    lambda.resize(2*joint_num);

    q.setZero();
    qdot.setZero();
    qddot.setZero();
    M.setZero();
    M_inv.setZero();
    J.setZero();
    J_dot.setZero();
    F_ext.setZero();
    lambda.setZero();
}

// Implement this function
bool dspSolver::CalculateConstraintError(Eigen::VectorXd& error_c) {
    /// calculate error what you define as constraint
    for(int i=0; i<LinkConfig.GetNumJoint(); i++){
        //two constraints per joint
        dspJoint* joint = LinkConfig.GetJoint(i);

        //for each joint, we need each link's pointer, their angles/positions, and the position of P (and Q if prismatic)
        int link1_id = joint->GetLink1ID(), link2_id = joint->GetLink2ID();
        dspLink* link1 = LinkConfig.GetLink(link1_id), *link2 = LinkConfig.GetLink(link2_id);

        double x1,y1,theta1, x2,y2,theta2;
        link1->GetQ(x1,y1,theta1);
        link2->GetQ(x2,y2,theta2);

        //get P1 and P2
        Eigen::Vector2d P1, P2;
        joint->GetP1_Link1_LCS(P1);
        joint->GetP2_Link2_LCS(P2);

        //get postion of P1 and P2 in GCS
        double P1x = x1 + P1.x()*cos(theta1) - P1.y()*sin(theta1);
        double P1y = y1 + P1.x()*sin(theta1) + P1.y()*cos(theta1);
        double P2x = x2 + P2.x()*cos(theta2) - P2.y()*sin(theta2);
        double P2y = y2 + P2.x()*sin(theta2) + P2.y()*cos(theta2);

        if(joint->GetType() == 0){
            //revolute joint

            //first constraint (x)
            double C1 = P1x - P2x; // = 0

            //second constraint (y)
            double C2 = P1y - P2y; // = 0

            error_c.segment(2*i, 2) << C1, C2;
        }

        else if(joint->GetType() == 1){
            //prismatic joint - get Q too
            Eigen::Vector2d Q1, Q2;
            joint->GetQ1_Link1_LCS(Q1);
            joint->GetQ2_Link2_LCS(Q2);

            //get postion of Q1 and Q2 in GCS
            double Q1x = x1 + Q1.x()*cos(theta1) - Q1.y()*sin(theta1);
            double Q1y = y1 + Q1.x()*sin(theta1) + Q1.y()*cos(theta1);
            double Q2x = x2 + Q2.x()*cos(theta2) - Q2.y()*sin(theta2);
            double Q2y = y2 + Q2.x()*sin(theta2) + Q2.y()*cos(theta2);

            //first constraint (translational contraint)
            double C1 = (P1y - Q1y)*(P2x - P1x) + (Q1x - P1x)*(P2y - P1y); // = 0

            //second constraint (rotational constraint)
            double C2 = (P1y - Q1y)*(Q2x - P2x) + (Q1x - P1x)*(Q2y - P2y); // = 0

            error_c.segment(2*i, 2) << C1, C2;
        }
    }
    return true;
}

// Implement this function
bool dspSolver::Make_M(void) { 
    /// make mass matrix : M
    for(int i=1; i<LinkConfig.GetNumLink(); i++){
        dspLink* link = LinkConfig.GetLink(i);
        double mass = link->GetMass();
        double inertia = link->GetInertia();
        int idx = i-1;
        M(3*idx,3*idx) = mass;
        M(3*idx+1,3*idx+1) = mass;
        M(3*idx+2,3*idx+2) = inertia;
    }
    return true;
}

// Implement this function
bool dspSolver::Make_J(void) {
    /// make Jacobian matrix which has first-order partial derivatives of vector q's components : J

//J has rows = #constraints and cols = #elements in q

    for(int i=0; i<LinkConfig.GetNumJoint(); i++){
        //two constraints per joint
        dspJoint* joint = LinkConfig.GetJoint(i);

        //for each joint, we need each link's pointer, their angles, and the position of P (and Q if prismatic)
        int link1_id = joint->GetLink1ID(), link2_id = joint->GetLink2ID();
        dspLink* link1 = LinkConfig.GetLink(link1_id), *link2 = LinkConfig.GetLink(link2_id);

        double x1,y1,theta1, x2,y2,theta2;
        link1->GetQ(x1,y1,theta1);
        link2->GetQ(x2,y2,theta2);

        //get P1 and P2
        Eigen::Vector2d P1, P2;
        joint->GetP1_Link1_LCS(P1);
        joint->GetP2_Link2_LCS(P2);
        double P1x = P1.x(), P1y = P1.y(), P2x = P2.x(), P2y = P2.y();

        double partial_x1, partial_y1, partial_theta1, partial_x2, partial_y2, partial_theta2;
        int row = 2*i;

        //if one of the two revolute joints is link0, then ignore the partial derivative terms of terms of that link.
        //The Jacobian should only contain the derivatives for link index 1 and above.
        //Each joint has 2 constraints (2 rows in the Jacobian)
        //The corresponding column index in the jacobian is 3*(link_idx-1) + 0,1,2

        if(joint->GetType() == 0){
            //revolute joint

            //C1(x) first row
            //link1
            partial_x1 = 1;
            partial_y1 = 0;
            partial_theta1 = - P1y*cos(theta1) - P1x*sin(theta1);
            //link2
            partial_x2 = -1;
            partial_y2 = 0;
            partial_theta2 = P2y*cos(theta2) + P2x*sin(theta2);

            if(link1_id !=0){
                J(row, 3*(link1_id-1)) = partial_x1;
                J(row, 3*(link1_id-1)+1) = partial_y1;
                J(row, 3*(link1_id-1)+2) = partial_theta1;
            }
            if(link2_id !=0){
                J(row, 3*(link2_id-1)) = partial_x2;
                J(row, 3*(link2_id-1)+1) = partial_y2;
                J(row, 3*(link2_id-1)+2) = partial_theta2;
            }

            //C2(y) second row
            row++;
            //link1
            partial_x1 = 0;
            partial_y1 = 1;
            partial_theta1 = P1x*cos(theta1) - P1y*sin(theta1);
            //link2
            partial_x2 = 0;
            partial_y2 = -1;
            partial_theta2 = P2y*sin(theta2) - P2x*cos(theta2);

            if(link1_id !=0){
                J(row, 3*(link1_id-1)) = partial_x1;
                J(row, 3*(link1_id-1)+1) = partial_y1;
                J(row, 3*(link1_id-1)+2) = partial_theta1;
            }
            if(link2_id !=0){
                J(row, 3*(link2_id-1)) = partial_x2;
                J(row, 3*(link2_id-1)+1) = partial_y2;
                J(row, 3*(link2_id-1)+2) = partial_theta2;
            }
        }

        else if(joint->GetType() == 1){
            //prismatic joint - get Q too
            Eigen::Vector2d Q1, Q2;
            joint->GetQ1_Link1_LCS(Q1);
            joint->GetQ2_Link2_LCS(Q2);
            double Q1x= Q1.x(), Q1y = Q1.y(), Q2x = Q2.x(), Q2y = Q2.y();

            //C1 first row
            //link1
            partial_x1 = Q1y*cos(theta1) - P1y*cos(theta1) - P1x*sin(theta1) + Q1x*sin(theta1);
            partial_y1 = P1x*cos(theta1) - Q1x*cos(theta1) - P1y*sin(theta1) + Q1y*sin(theta1);
            partial_theta1 = (P1y*cos(theta1) + P1x*sin(theta1))*(P1y*cos(theta1) - Q1y*cos(theta1) + P1x*sin(theta1) - Q1x*sin(theta1)) + (P1x*cos(theta1) - P1y*sin(theta1))*(P1x*cos(theta1) - Q1x*cos(theta1) - P1y*sin(theta1) + Q1y*sin(theta1)) - (P1x*cos(theta1) - Q1x*cos(theta1) - P1y*sin(theta1) + Q1y*sin(theta1))*(x1 - x2 + P1x*cos(theta1) - P2x*cos(theta2) - P1y*sin(theta1) + P2y*sin(theta2)) - (P1y*cos(theta1) - Q1y*cos(theta1) + P1x*sin(theta1) - Q1x*sin(theta1))*(y1 - y2 + P1y*cos(theta1) - P2y*cos(theta2) + P1x*sin(theta1) - P2x*sin(theta2));
            //link2
            partial_x2 = P1y*cos(theta1) - Q1y*cos(theta1) + P1x*sin(theta1) - Q1x*sin(theta1);
            partial_y2 = Q1x*cos(theta1) - P1x*cos(theta1) + P1y*sin(theta1) - Q1y*sin(theta1);
            partial_theta2 = - (P2y*cos(theta2) + P2x*sin(theta2))*(P1y*cos(theta1) - Q1y*cos(theta1) + P1x*sin(theta1) - Q1x*sin(theta1)) - (P2x*cos(theta2) - P2y*sin(theta2))*(P1x*cos(theta1) - Q1x*cos(theta1) - P1y*sin(theta1) + Q1y*sin(theta1));

            if(link1_id !=0){
                J(row, 3*(link1_id-1)) = partial_x1;
                J(row, 3*(link1_id-1)+1) = partial_y1;
                J(row, 3*(link1_id-1)+2) = partial_theta1;
            }
            if(link2_id !=0){
                J(row, 3*(link2_id-1)) = partial_x2;
                J(row, 3*(link2_id-1)+1) = partial_y2;
                J(row, 3*(link2_id-1)+2) = partial_theta2;
            }

            //C2 second row
            row++;
            //link1
            partial_x1 = 0;
            partial_y1 = 0;
            partial_theta1 = - (P1x*cos(theta1) - Q1x*cos(theta1) - P1y*sin(theta1) + Q1y*sin(theta1))*(P2x*cos(theta2) - Q2x*cos(theta2) - P2y*sin(theta2) + Q2y*sin(theta2)) - (P1y*cos(theta1) - Q1y*cos(theta1) + P1x*sin(theta1) - Q1x*sin(theta1))*(P2y*cos(theta2) - Q2y*cos(theta2) + P2x*sin(theta2) - Q2x*sin(theta2));
            //link2
            partial_x2 = 0;
            partial_y2 = 0;
            partial_theta2 = (P1x*cos(theta1) - Q1x*cos(theta1) - P1y*sin(theta1) + Q1y*sin(theta1))*(P2x*cos(theta2) - Q2x*cos(theta2) - P2y*sin(theta2) + Q2y*sin(theta2)) + (P1y*cos(theta1) - Q1y*cos(theta1) + P1x*sin(theta1) - Q1x*sin(theta1))*(P2y*cos(theta2) - Q2y*cos(theta2) + P2x*sin(theta2) - Q2x*sin(theta2));
            //zero for the rest

            if(link1_id !=0){
                J(row, 3*(link1_id-1)) = partial_x1;
                J(row, 3*(link1_id-1)+1) = partial_y1;
                J(row, 3*(link1_id-1)+2) = partial_theta1;
            }
            if(link2_id !=0){
                J(row, 3*(link2_id-1)) = partial_x2;
                J(row, 3*(link2_id-1)+1) = partial_y2;
                J(row, 3*(link2_id-1)+2) = partial_theta2;
            }
        }
    }
    return true;
}

// Implement this function
bool dspSolver::Make_J_dot(void) {
    /// make time derivative of Jacobian matrix  : J_dot

//J has rows = #constraints and cols = #elements in q

    for(int i=0; i<LinkConfig.GetNumJoint(); i++){
        //two constraints per joint
        dspJoint* joint = LinkConfig.GetJoint(i);

        //for each joint, we need each link's pointer, their angles, and the position of P (and Q if prismatic)
        int link1_id = joint->GetLink1ID(), link2_id = joint->GetLink2ID();
        dspLink *link1 = LinkConfig.GetLink(link1_id), *link2 = LinkConfig.GetLink(link2_id);

        double x1,y1,theta1, x2,y2,theta2;
        link1->GetQ(x1,y1,theta1);
        link2->GetQ(x2,y2,theta2);

        double x1dot,y1dot,theta1dot, x2dot,y2dot,theta2dot;
        link1->GetQDot(x1,y1,theta1);
        link2->GetQDot(x2,y2,theta2);

        //get P1 and P2
        Eigen::Vector2d P1, P2;
        joint->GetP1_Link1_LCS(P1);
        joint->GetP2_Link2_LCS(P2);
        double P1x = P1.x(), P1y = P1.y(), P2x = P2.x(), P2y = P2.y();

        double partial_x1_dt, partial_y1_dt, partial_theta1_dt, partial_x2_dt, partial_y2_dt, partial_theta2_dt;
        int row = 2*i;

        if(joint->GetType() == 0){
            //revolute joint

            //C1 first row
            //link 1
            partial_x1_dt = 0;
            partial_y1_dt = 0;
            partial_theta1_dt = (P1y*sin(theta1) - P1x*cos(theta1))*theta1dot;
            //link 2
            partial_x2_dt = 0;
            partial_y2_dt = 0;
            partial_theta2_dt = -(P2y*sin(theta2) - P2x*cos(theta2))*theta2dot;

            if(link1_id !=0){
                J_dot(row, 3*(link1_id-1)) = partial_x1_dt;
                J_dot(row, 3*(link1_id-1)+1) = partial_y1_dt;
                J_dot(row, 3*(link1_id-1)+2) = partial_theta1_dt;
            }
            if(link2_id !=0){
                J_dot(row, 3*(link2_id-1)) = partial_x2_dt;
                J_dot(row, 3*(link2_id-1)+1) = partial_y2_dt;
                J_dot(row, 3*(link2_id-1)+2) = partial_theta2_dt;
            }

            //C2 second row
            row++;
            //link 1
            partial_x1_dt = 0;
            partial_y1_dt = 0;
            partial_theta1_dt = -(P1x*sin(theta1) + P1y*cos(theta1))*theta1dot;
            //link 2
            partial_x2_dt = 0;
            partial_y2_dt = 0;
            partial_theta2_dt = (P2x*sin(theta2) + P2y*cos(theta2))*theta2dot;

            if(link1_id !=0){
                J_dot(row, 3*(link1_id-1)) = partial_x1_dt;
                J_dot(row, 3*(link1_id-1)+1) = partial_y1_dt;
                J_dot(row, 3*(link1_id-1)+2) = partial_theta1_dt;
            }
            if(link2_id !=0){
                J_dot(row, 3*(link2_id-1)) = partial_x2_dt;
                J_dot(row, 3*(link2_id-1)+1) = partial_y2_dt;
                J_dot(row, 3*(link2_id-1)+2) = partial_theta2_dt;
            }
        }

        else if(joint->GetType() == 1){
            //prismatic joint - get Q too
            Eigen::Vector2d Q1, Q2;
            joint->GetQ1_Link1_LCS(Q1);
            joint->GetQ2_Link2_LCS(Q2);
            double Q1x= Q1.x(), Q1y = Q1.y(), Q2x = Q2.x(), Q2y = Q2.y();

            //C1 first row
            //link 1
            partial_x1_dt = theta1dot*(P1y*sin(theta1) - Q1y*sin(theta1) - P1x*cos(theta1) + Q1x*cos(theta1));
            partial_y1_dt = -theta1dot*(P1x*sin(theta1) - Q1x*sin(theta1) + P1y*cos(theta1) - Q1y*cos(theta1));
            partial_theta1_dt = (P1y*sin(theta1) - Q1y*sin(theta1) - P1x*cos(theta1) + Q1x*cos(theta1))*(x1dot - x2dot - P1y*cos(theta1)*theta1dot + P2y*cos(theta2)*theta2dot - P1x*sin(theta1)*theta1dot + P2x*sin(theta2)*theta2dot) - (P1x*sin(theta1) - Q1x*sin(theta1) + P1y*cos(theta1) - Q1y*cos(theta1))*(y1dot - y2dot + P1x*cos(theta1)*theta1dot - P2x*cos(theta2)*theta2dot - P1y*sin(theta1)*theta1dot + P2y*sin(theta2)*theta2dot) + (P1x*sin(theta1) - Q1x*sin(theta1) + P1y*cos(theta1) - Q1y*cos(theta1))*(x1 - x2 - P1y*sin(theta1) + P2y*sin(theta2) + P1x*cos(theta1) - P2x*cos(theta2))*theta1dot + (P1y*sin(theta1) - Q1y*sin(theta1) - P1x*cos(theta1) + Q1x*cos(theta1))*(y1 - y2 + P1x*sin(theta1) - P2x*sin(theta2) + P1y*cos(theta1) - P2y*cos(theta2))*theta1dot;
            //link 2
            partial_x2_dt = -theta1dot*(P1y*sin(theta1) - Q1y*sin(theta1) - P1x*cos(theta1) + Q1x*cos(theta1));
            partial_y2_dt = theta1dot*(P1x*sin(theta1) - Q1x*sin(theta1) + P1y*cos(theta1) - Q1y*cos(theta1));
            partial_theta2_dt = (P2x*sin(theta2) + P2y*cos(theta2))*(P1y*sin(theta1) - Q1y*sin(theta1) - P1x*cos(theta1) + Q1x*cos(theta1))*theta1dot - (P2x*sin(theta2) + P2y*cos(theta2))*(P1y*sin(theta1) - Q1y*sin(theta1) - P1x*cos(theta1) + Q1x*cos(theta1))*theta2dot - (P2y*sin(theta2) - P2x*cos(theta2))*(P1x*sin(theta1) - Q1x*sin(theta1) + P1y*cos(theta1) - Q1y*cos(theta1))*theta1dot + (P2y*sin(theta2) - P2x*cos(theta2))*(P1x*sin(theta1) - Q1x*sin(theta1) + P1y*cos(theta1) - Q1y*cos(theta1))*theta2dot;

            if(link1_id !=0){
                J_dot(row, 3*(link1_id-1)) = partial_x1_dt;
                J_dot(row, 3*(link1_id-1)+1) = partial_y1_dt;
                J_dot(row, 3*(link1_id-1)+2) = partial_theta1_dt;
            }
            if(link2_id !=0){
                J_dot(row, 3*(link2_id-1)) = partial_x2_dt;
                J_dot(row, 3*(link2_id-1)+1) = partial_y2_dt;
                J_dot(row, 3*(link2_id-1)+2) = partial_theta2_dt;
            }

            //C2:
            row++;
            //link 1
            partial_x1_dt = 0;
            partial_y1_dt = 0;
            partial_theta1_dt = theta1dot*(P1y*sin(theta1) - Q1y*sin(theta1) - P1x*cos(theta1) + Q1x*cos(theta1))*(P2x*sin(theta2) - Q2x*sin(theta2) + P2y*cos(theta2) - Q2y*cos(theta2)) - theta1dot*(P1x*sin(theta1) - Q1x*sin(theta1) + P1y*cos(theta1) - Q1y*cos(theta1))*(P2y*sin(theta2) - Q2y*sin(theta2) - P2x*cos(theta2) + Q2x*cos(theta2)) + theta2dot*(P1x*sin(theta1) - Q1x*sin(theta1) + P1y*cos(theta1) - Q1y*cos(theta1))*(P2y*sin(theta2) - Q2y*sin(theta2) - P2x*cos(theta2) + Q2x*cos(theta2)) - theta2dot*(P1y*sin(theta1) - Q1y*sin(theta1) - P1x*cos(theta1) + Q1x*cos(theta1))*(P2x*sin(theta2) - Q2x*sin(theta2) + P2y*cos(theta2) - Q2y*cos(theta2));

            //link 2
            partial_x2_dt = 0;
            partial_y2_dt = 0;
            partial_theta2_dt = (P1x*sin(theta1) - Q1x*sin(theta1) + P1y*cos(theta1) - Q1y*cos(theta1))*(P2y*sin(theta2) - Q2y*sin(theta2) - P2x*cos(theta2) + Q2x*cos(theta2))*theta1dot - (P1y*sin(theta1) - Q1y*sin(theta1) - P1x*cos(theta1) + Q1x*cos(theta1))*(P2x*sin(theta2) - Q2x*sin(theta2) + P2y*cos(theta2) - Q2y*cos(theta2))*theta1dot - (P1x*sin(theta1) - Q1x*sin(theta1) + P1y*cos(theta1) - Q1y*cos(theta1))*(P2y*sin(theta2) - Q2y*sin(theta2) - P2x*cos(theta2) + Q2x*cos(theta2))*theta2dot + (P1y*sin(theta1) - Q1y*sin(theta1) - P1x*cos(theta1) + Q1x*cos(theta1))*(P2x*sin(theta2) - Q2x*sin(theta2) + P2y*cos(theta2) - Q2y*cos(theta2))*theta2dot;
            
            if(link1_id !=0){
                J_dot(row, 3*(link1_id-1)) = partial_x1_dt;
                J_dot(row, 3*(link1_id-1)+1) = partial_y1_dt;
                J_dot(row, 3*(link1_id-1)+2) = partial_theta1_dt;
            }
            if(link2_id !=0){
                J_dot(row, 3*(link2_id-1)) = partial_x2_dt;
                J_dot(row, 3*(link2_id-1)+1) = partial_y2_dt;
                J_dot(row, 3*(link2_id-1)+2) = partial_theta2_dt;
            }
        }
    }
    return true;
}

// Implement this function
bool dspSolver::Make_F_ext(void) {
    /// make a n-dimensional vector for external force/torque : F_ext
    /// it contains x direction of force, y direction of force, and torque

    for(int i=1; i<LinkConfig.GetNumLink(); i++){
        dspLink* link = LinkConfig.GetLink(i);
        //set gravity force
        double mass = link->GetMass();
        F_ext[3*(i-1)+1] = -9.81*mass;
    }

    return true;
}

// Implement this function
bool dspSolver::CalcLinAlg(void) {
    /// calculate second derivative of vector q : qddot

    M_inv = M.inverse();
    lambda = (J*M_inv*J.transpose()).inverse()* (-J_dot*qdot - J*M_inv*F_ext);
    qddot = M_inv*(F_ext + J.transpose()*lambda);

    return true;
}

// Implement this function
bool dspSolver::UpdateCurrentInfo(void) {
    /// load and update information what you need in class variables
    /// current vector q and its derivative should be loaded : q, q_dot

    for(int i=1; i<LinkConfig.GetNumLink(); i++){
        dspLink* link = LinkConfig.GetLink(i);
        double x,y,theta;

        //load q from links
        link->GetQ(x,y,theta);
        q.segment(3*(i-1),3) << x,y,theta;

        //load qdot from links
        link->GetQDot(x,y,theta);
        qdot.segment(3*(i-1),3) << x,y,theta;
    }

    Make_M();
    Make_J();
    Make_J_dot();
    Make_F_ext();
    return true;
}

// Implement this function
bool dspSolver::UpdateNextInfo(double timestep) {
    /// update vector q and its derivative : q, q_dot
    /// also save them in class variables for next step
    qdot += qddot*timestep;
    q += qdot*timestep;

    SetQ();
    SetQDot();

    return true;
}

// Implement this function
bool dspSolver::SetQDot(void) {
    /// update new q_dot information in each link structure
    for(int i=1; i<LinkConfig.GetNumLink(); i++){
        dspLink* link = LinkConfig.GetLink(i);
        double velx, vely, w;
        velx = qdot(3*(i-1));
        vely = qdot(3*(i-1)+1);
        w = qdot(3*(i-1)+2);
        link->SetQDot(velx, vely, w);
    }

    return true;
}

// Implement this function
bool dspSolver::SetQ(void) {
    /// update new q information in each link structure
    for(int i=1; i<LinkConfig.GetNumLink(); i++){
        dspLink* link = LinkConfig.GetLink(i);
        double x,y,theta;
        x = q(3*(i-1));
        y = q(3*(i-1)+1);
        theta = q(3*(i-1)+2);
        link->SetQ(x,y,theta);
    }

    return true;
}

// Do not change this function
bool dspSolver::GetQ(Eigen::VectorXd& q_record) {
    /// get current q information
    q_record = q;
    return true;
}

// Do not change this function
bool dspSolver::GetQDot(Eigen::VectorXd& qdot_record) {
    /// get current q_dot information 
    qdot_record = qdot;
    return true;
}
