
//TODO: implement 2dim derivatives and forces (-> currently not corract)





#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <random>


using  namespace Eigen;//Eigen::Matrix2Xd, Eigen::MatrixXd,Eigen::MatrixXi,Eigen::VectorXd;


const double m=1.0;
const int N=10;
const double h=0.2;
const double mu=1.;
const double cutOff=.2;
const double k=1.;
const double rho_0=0.1;
const double g=9.81;
const int N_steps=10;

void searchNeighbour(int i,MatrixX2d& x,MatrixXi& neighbours)
{
    Vector2d temp;
	//TODO: adapt to cell list!
    for (int j=0; j< N ; ++j) { //loop over all particles
        if (j!=i)
        {
            temp=Vector2d(x(i,0)-x(j,0),x(i,1)-x(j,1));
            if(temp.norm()<cutOff){// if ((x(i)-x(j)).norm())<cutOff) {
                neighbours(i,j)=1;
                neighbours(j,i)=1;
            }else{  neighbours(i,j)=0;
                neighbours(j,i)=0;}
        }
    }
	
}



double w(double x_ij) //kernel poly6
{
    if((x_ij<=h)&&(x_ij>=0))
    {
        return 315/(64*M_PI*std::pow(h,9))*std::pow((h*h-x_ij*x_ij),3);
    }
    return 0.;
}

double w_grad(double x_ij) //derivative of poly6 kernel
{
  
    if((x_ij<=h)&&(x_ij>=0))
    {
        return -315*6/(64*M_PI*std::pow(h,9))*std::pow((h*h-x_ij*x_ij),2)*x_ij;
    }
    return 0.;
}
                       
double w_grad2(double x_ij)
{
    if((x_ij<=h)&&(x_ij>=0))
    {
        return 315/(64*M_PI*std::pow(h,9))*(-std::pow((h*h-x_ij*x_ij),2)*6*x_ij+(h*h-x_ij*x_ij)*24*x_ij*x_ij);
    }
    return 0.;
}

void computeDensity(VectorXd& rho,int i,MatrixX2d& x,MatrixXi& neighbours)
{
    Vector2d temp;
    for (int j; j<N; ++j) {
        temp = Vector2d(x(i,0)-x(j,0),x(i,1)-x(j,1));
        rho(i)+=m* w(temp.norm()) *neighbours(i,j); //check with ,neighobours' that other particle is within i.a. radius
        
    }
}

void computePressure(VectorXd& p,int i,Vector2d& rho)
{
    double rho_0=0.1;
    p(i)=k*(rho(i)-rho_0);
}

void computeForce(int i,MatrixX2d& f,VectorXd& p,MatrixX2d& u, MatrixX2d& x,VectorXd& rho,MatrixXi& neighbours)
{
    
    VectorXd f_grav(p.size());
    
    MatrixX2d f_vis=MatrixX2d::Zero(p.size(),2);
    MatrixX2d f_p(p.size(),2);
    Vector2d temp(2);
    
 //   std::cout << "f_vis "<< f_vis.size() << "\n";
   //  std::cout << "u "<< u.size() << "\n";
     //std::cout << "rho "<< rho.size() << "\n";
    for (int j=0; j<N; ++j)
    {
        if (neighbours(i,j)==1)
        {
            temp(0)=x(i,0)-x(j,0);
            temp(1)=x(i,1)-x(j,1);
            f_p(i,1)+= m/rho(j)*(p(i)+p(j))/2.*w_grad(std::abs(x(i,1)-x(j,1)));
            f_p(i,0)+= m/rho(j)*(p(i)+p(j))/2.*w_grad(std::abs(x(i,0)-x(j,0)));

           f_vis(i,1)+=mu*(m/rho(j))*u(j,1)-u(i,1)*w_grad2(std::abs(x(i,1)-x(j,1)));
            f_vis(i,0)+=mu*(m/rho(j))*u(j,0)-u(i,0)*w_grad2(std::abs(x(i,0)-x(j,0)));

        }
        
    }
    for (int j=0;j<N; ++j) {
        f(i,1)= -f_p(i,1) +f_vis(i,1) ;//+f_grav(i);
        f(i,0)= -f_p(i,0) +f_vis(i,0); //+f_grav(i);

    }
}
                       
                       
void Step(MatrixX2d& u,MatrixX2d& x,double dt,MatrixXi& neighbours,VectorXd& p,VectorXd& rho,MatrixX2d& f)
{
    for (int i=0; i<N; ++i) {
        searchNeighbour(i,x,neighbours);
    }
    for (int i=0; i<N; ++i) {
        computeDensity(rho,i,x,neighbours);
       p(i)=k*(rho(i)-rho_0);// computePressure(i);
        
    }
  //  std::cout <<p.size();
    for (int i=0; i<N; ++i)
    {
        computeForce( i,f,p,u,x,rho,neighbours);
    }
    for (int i=0; i<N; ++i) {
        u(i,1) += dt*f(i,1)/rho(i);
        u(i,0) += dt*f(i,0)/rho(i);
        
        if((x(i,1)+ dt*u(i,1))<0) {     //collision detection with boundaries
          x(i,1)=0.;
        }else
        {
            x(i,1) += dt*u(i,1);
        }
        
        if((x(i,0)+ dt*u(i,0))<0) {
            x(i,0)=0.;
        }else
        {
            x(i,0) += dt*u(i,0);
        }
    }
    
    
}


int main()
{
    std::mt19937 mt_rand(42);
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    int x_size;
    MatrixX2d u=MatrixX2d::Zero(N,2);//= MatrixX2d::Random(N,2);//std::vector<double> u(2*N);
    MatrixX2d x=MatrixX2d::Zero(N,2);
    
    VectorXd rho(N);
    VectorXd p(N);
    
    for (int i=0; i<N; ++i) {
        x(i,0)=dis(mt_rand);
        x(i,1)=dis(mt_rand);
        u(i,0)=dis(mt_rand);
        u(i,1)=dis(mt_rand);
        p(i)=rho(i)=1.;
    }
    MatrixXi neighbours=MatrixXi::Zero(N,N);//std::vector<double> neighbours(N*N,0);

    MatrixX2d f(N,2);
    
    for (int i=0; i<N; ++i) {
        std::cout <<"( "<< x(i,0)<<", " <<x(i,1)<<" )"<<std::endl;
    }
    for (int s=0; s<N_steps; ++s)
    {
        Step(u,x,0.1,neighbours,p,rho,f);
    }
    std::cout <<"\n\n ";
    for (int i=0; i<N; ++i) {
            std::cout <<"( "<< x(i,0)<<", " <<neighbours(i,1)<<" )"<<std::endl;
    }
    /*for (int i=0; i<N; ++i) {
        std::cout <<"( "<<neighbours(i,1)<<", " <<neighbours(i,2)<<", " <<neighbours(i,3)<<", " <<neighbours(i,4)<<", " <<neighbours(i,5)<<", " <<neighbours(i,6)<<", " <<neighbours(i,7)<<" )"<<std::endl;
    }*/
    return 0;
}