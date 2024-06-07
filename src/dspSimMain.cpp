#include <iostream>
#include <fstream>
#include "dspLinkConfig.h"
#include "dspSolver.h"

int main(void)
{
    std::ofstream q_results;
    std::ofstream qdot_results;
    std::ofstream error_results;
    dspLinkConfig LinkConfig;
    /// index of the ground link is fixed to 0

    // /// Example of linkage configuration for a four bar linkage
    // LinkConfig.AddLink(0, 1000000.0, 1000000.0);
    // LinkConfig.AddLink(1, 1.0, 1.0 / 12.0 * 1.0 * 1.0 * 1.0);
    // LinkConfig.AddLink(2, 1.0, 1.0 / 12.0 * 1.0 * 4.0 * 4.0);
    // LinkConfig.AddLink(3, 1.0, 1.0 / 12.0 * 1.0 * 2.5 * 2.5);
    // double theta1 = 1.578;
    // double theta2 = 0.3533;
    // double theta3 = 1.2649 + EIGEN_PI;
    // LinkConfig.GetLink(0)->SetQ(0, 0, 0);
    // LinkConfig.GetLink(1)->SetQ(0.5 * cos(theta1), 0.5 * sin(theta1), theta1);
    // LinkConfig.GetLink(2)->SetQ(1.0 * cos(theta1) + 2.0 * cos(theta2), 1.0 * sin(theta1) + 2.0 * sin(theta2), theta2);
    // LinkConfig.GetLink(3)->SetQ(1.0 * cos(theta1) + 4.0 * cos(theta2) + 1.25 * cos(theta3), 1.0 * sin(theta1) + 4.0 * sin(theta2) + 1.25 * sin(theta3), theta3);
    // LinkConfig.GetLink(1)->SetQDot(0, 0, 0);
    // LinkConfig.GetLink(2)->SetQDot(0, 0, 0);
    // LinkConfig.GetLink(3)->SetQDot(0, 0, 0);
    // LinkConfig.AddJointRevolute(0, 0, 0.0, 0.0, 1, -0.5, 0.0);
    // LinkConfig.AddJointRevolute(1, 1, 0.5, 0.0, 2, -2.0, 0.0);
    // LinkConfig.AddJointRevolute(2, 2, 2.0, 0.0, 3, -1.25, 0.0);
    // LinkConfig.AddJointRevolute(3, 3, 1.25, 0.0, 0, 3.0, 0.0);

    // Example of linkage configuration for a crank-shaft system
    LinkConfig.AddLink(0, 1000000.0, 1000000.0);
    LinkConfig.AddLink(1, 1.0, 1.0 );
    LinkConfig.AddLink(2, 1.0, 1.0 );
    LinkConfig.AddLink(3, 2.0, 1.30125); // composite inertia 1.30125 = link_box : 0.270, link_ground : 1, distance :0.25
    double theta1 = EIGEN_PI / 2 - 0.5;
    double theta2 = - (EIGEN_PI / 2 - (EIGEN_PI - 1.25 - 0.5));
    double theta3 = 0;
    LinkConfig.GetLink(0)->SetQ(0, 0, 0);
    LinkConfig.GetLink(1)->SetQ(0.5 * cos(theta1), 0.1 + 0.5 * sin(theta1), theta1);
    LinkConfig.GetLink(2)->SetQ(1.0 * cos(theta1) + 2.0 * cos(theta2), 0.1 + 1.0 * sin(theta1) + 2.0 * sin(theta2), theta2);
    LinkConfig.GetLink(3)->SetQ(1.0 * cos(theta1) + 4.0 * cos(theta2) + 0.0 * cos(theta3), 0.1 + 1.0 * sin(theta1) + 4.0 * sin(theta2) + 0.0 * sin(theta3), theta3);
    LinkConfig.GetLink(1)->SetQDot(0, 0, 0);
    LinkConfig.GetLink(2)->SetQDot(0, 0, 0);
    LinkConfig.GetLink(3)->SetQDot(0, 0, 0);
    LinkConfig.AddJointRevolute(0, 0, 0.0, 0.1, 1, -0.5, 0.0);
    LinkConfig.AddJointRevolute(1, 1, 0.5, 0.0, 2, -2.0, 0.0);
    LinkConfig.AddJointRevolute(2, 2, 2.0, 0.0, 3, 0.0, 0.0);
    LinkConfig.AddJointPrismatic(3, 3, 0.0, -(0.1 + 1.0 * sin(theta1) + 4.0 * sin(theta2) + 0.0 * sin(theta3)),
                                0, 1.0 * cos(theta1) + 4.0 * cos(theta2) + 0.0 * cos(theta3), 0.0,
                                0.0); //angle


    dspSolver Solver(LinkConfig);
    Eigen::VectorXd Constraint_error;
    Constraint_error.setZero(2 * LinkConfig.GetNumJoint());
    Eigen::VectorXd q_record(3*(LinkConfig.GetNumLink()-1)), qdot_record(3*(LinkConfig.GetNumLink()-1)); // except ground


    /// test time = 20 sec, record period = 0.01 sec, time step = 0.001 sec
    /// if you don't initialize variables in dspSolver class, recording will not be done.
    int n_sim = 20000;
    bool record_ = true;
    double t = 0, h = 0.001;

    if(record_){
        q_results.open("q_results_cpp.csv");
        qdot_results.open("qdot_results_cpp.csv");
        error_results.open("c_results_cpp.csv");
    }

    std::cout << "Simulation Start" << std::endl;
    for (int i = 0; i < n_sim; i++) {
        Solver.UpdateCurrentInfo();  // q_pre, qdot_pre
        Solver.CalcLinAlg();         // qddot_next
        Solver.UpdateNextInfo(h);     // q_next, qdot_next
        if (record_ and i % 10 == 0)
        {
            Solver.GetQ(q_record);
            Solver.GetQDot(qdot_record);
            Solver.CalculateConstraintError(Constraint_error);

            std::cout << t << " sec" << std::endl;
            q_results << t << ",";
            qdot_results << t << ",";
            error_results << t << ",";
            for (int rec_i = 0; rec_i < 3 * (LinkConfig.GetNumLink() - 1); rec_i++) {
                q_results << q_record(rec_i) << ",";
            }
            q_results << "\n";
            for (int rec_i = 0; rec_i < 3 * (LinkConfig.GetNumLink() - 1); rec_i++) {
                qdot_results << qdot_record(rec_i) << ",";
            }
            qdot_results << "\n";
            error_results << Constraint_error(0);
            for (int j = 1; j < 2 * LinkConfig.GetNumJoint(); j++) {
                error_results << "," << Constraint_error(j);
            }error_results << "\n";
        }
        t = t + h;
    }
    if(record_){
        q_results.close();
        qdot_results.close();
        error_results.close();
    }
    std::cout << "Simulation End"<<std::endl;

    return 0;
}