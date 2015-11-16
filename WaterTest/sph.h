#include <iostream>
#include "GLUT/glut.h"
#include <cmath>
#include <Eigen/core>
#include <Eigen/Geometry>
#include <random>

class sph {
public:
    sph();
  //  int pppp;
    Eigen::MatrixX2d x;
    Eigen::MatrixX2d u;
    Eigen::MatrixXi neighbours;
    Eigen::VectorXd rho;
    Eigen::VectorXd p;
    Eigen::MatrixX2d f;

    void render();
    void update();
    
    void step();
    void searchNeighbour(int i);
    double w(double x_ij);
    double w_grad(double x_ij);
    double w_grad2(double x_ij);
    void computeDensity(int i);
    void computePressure(int i);
    void computeForce(int i);
    
    
};
