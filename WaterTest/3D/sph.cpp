
//TODO: implement 2dim derivatives and forces (-> currently not corract)

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <random>
#include "sph.h"


using  namespace Eigen;//Eigen::Matrix2Xd, Eigen::MatrixXd,Eigen::MatrixXi,Eigen::VectorXd;



const int N=200;
const double d=1.;
const double h=2*d;
const double rho_0=1000;
const double m=d*d*d*rho_0;
const double mu=10;
const double cutOff=h;
const int boundary=11;
const double k=1100;

const double g=9.81;
const int N_steps=10;
const double dt=0.01; //1/h



sph::sph()//MatrixX2d& x, MatrixX2d& u, VectorXd& rho,VectorXd& p,MatrixXi& neighbours)
{
    std::mt19937 mt_rand(42);
    std::uniform_real_distribution<double> dis(0.0, 1.0);
   
    u=MatrixX3d::Zero(N,3);//= MatrixX2d::Random(N,2);//std::vector<double> u(2*N);
     x=MatrixX3d::Zero(N,3);
    f=MatrixX3d::Zero(N,3);

    rho=VectorXd::Zero(N);
    p=VectorXd::Zero(N);
    double q=0,qq=0;
    for (int i=0; i<N; ++i) {
        x(i,0)=d*(i%10);
        if(i%10==0){++q;}
        if(i%50==0)
        {
            ++qq;
            q=0;
        }
        x(i,1)=d*(q);
        x(i,2)=d*((qq));
        u(i,0)=10*dis(mt_rand);
        u(i,1)=5*dis(mt_rand);
        u(i,2)=0.;//dis(mt_rand);
        p(i)=0.;
        rho(i)=rho_0;
        //std::cout <<x(i,0);
    }
     neighbours=MatrixXi::Zero(N,N);
   // pppp=1000;
}

void sph::searchNeighbour(int i)//,MatrixX2d& x,MatrixXi& neighbours)
{
    Vector3d temp=Vector3d::Zero(3);
	//TODO: adapt to cell list!
    for (int j=0; j< N ; ++j) { //loop over all particles
        if (j!=i)
        {
         //   std::cout << x.size() <<"\n";
           // temp=Vector3d(x(i,0)-x(j,0),x(i,1)-x(j,1),x(i,2)-x(j,2));
            temp(0) =x(i,0)-x(j,0);
            temp(1)=x(i,1)-x(j,1);
            temp(2)=x(i,2)-x(j,2);
            if(temp.norm()<cutOff){// if ((x(i)-x(j)).norm())<cutOff) {
                neighbours(i,j)=1;
                neighbours(j,i)=1;
            }else{ neighbours(i,j)=0;
                neighbours(j,i)=0;}
        }
    }
	
}



double sph::w(double x_ij) //kernel
{
    if((x_ij<=h)&&(x_ij>=0))
    {
        return 315/(64*M_PI*std::pow(h,9))*std::pow((h*h-x_ij*x_ij),3); //poly6
       // return 15./(M_PI*std::pow(h,6))*std::pow(std::abs(h-(x_ij)),3); //spiky
    }
    return 0.;
}

double sph::w_grad(double x_ij) //derivative of poly6 kernel
{
  
    if((x_ij<=h)&&(x_ij>=0))
    {
        //return -315*6/(64*M_PI*std::pow(h,9))*std::pow((h*h-x_ij*x_ij),2)*x_ij;
       // return -45./(M_PI*std::pow(h,6))*1./x_ij*std::pow((h-x_ij),2);
        return -45./(M_PI*std::pow(h,6))*std::pow((h-x_ij),2);
    }
    return 0.;
}
                       
double sph::w_grad2(double x_ij)
{
    if((x_ij<=h)&&(x_ij>=0))
    {
       // return 315/(64*M_PI*std::pow(h,9))*(-std::pow((h*h-x_ij*x_ij),2)*6*x_ij+(h*h-x_ij*x_ij)*24*x_ij*x_ij); //poly6
       // return -90./(M_PI*std::pow(h,6))*1./(std::abs(x_ij))*(h-std::abs(x_ij))*(h-2*std::abs(x_ij)); //spiky
        return 45./(M_PI*std::pow(h,6))*(h-x_ij);
    }
    return 0.;
}

void sph::computeDensity(int i)//VectorXd& rho,int i,MatrixX2d& x,MatrixXi& neighbours)
{
    
    Vector3d temp=Vector3d::Zero(3);
    double rho_temp=0.;
    rho(i)=0.;
    for (int j=0; j<N; ++j)
    {
        if ((neighbours(i,j)==1)&& (j!=i))
        {
            temp(0) =x(i,0)-x(j,0);
            temp(1)=x(i,1)-x(j,1);
            temp(2)=x(i,2)-x(j,2);
            rho(i)+=m*w(temp.norm()); //check with ,neighobours' that other particle is within i.a. radius
        }

    }
   // rho(i)=rho_temp;
}

void sph::computePressure(int i)//VectorXd& p,int i,Vector2d& rho)
{
    if((rho(i)-rho_0)>0.){
        p(i)=k*(rho(i)-rho_0);
    }else{p(i)=0.;}
}

void sph::computeForce(int i)//,MatrixX2d& f,VectorXd& p,MatrixX2d& u, MatrixX2d& x,VectorXd& rho,MatrixXi& neighbours)
{
    
    VectorXd f_grav(p.size());
    
    MatrixX3d f_vis=MatrixX3d::Zero(p.size(),3);
    MatrixX3d f_p=MatrixX3d::Zero(p.size(),3);
    Vector3d temp= Vector3d::Zero(3);
    
    for (int j=0; j<N; ++j)
    {
        if ((neighbours(i,j)==1)&& (j!=i))
        {
            //temp = Vector3d(x(i,0)-x(j,0),x(i,1)-x(j,1),x(i,2)-x(j,2));
            temp(0) =x(i,0)-x(j,0);
            temp(1)=x(i,1)-x(j,1);
            temp(2)=x(i,2)-x(j,2);
            
            f_p(i,1)+= m*(p(i)+p(j))/(rho(j)*2.)*(x(i,1)-x(j,1))*w_grad(temp.norm());//std::abs(x(i,1)-x(j,1)));
            f_p(i,0)+= m*(p(i)+p(j))/(rho(j)*2.)*(x(i,0)-x(j,0))*w_grad(temp.norm());
            
            f_p(i,2)+=m*(p(i)+p(j))/(rho(j)*2.)*(x(i,2)-x(j,2))*w_grad(temp.norm());
            
            
            
            f_vis(i,1)+=mu*(m/rho(j))*(u(j,1)-u(i,1))*w_grad2(temp.norm());
            f_vis(i,0)+=mu*(m/rho(j))*(u(j,0)-u(i,0))*w_grad2(temp.norm());
            f_vis(i,2)+=mu*(m/rho(j))*(u(j,2)-u(i,2))*w_grad2(temp.norm());
            
        }
        
    }
    
    //  for (int j=0;j<N; ++j) {
    f(i,1)= -f_p(i,1) +f_vis(i,1) ;//+f_grav(i);
    f(i,0)= -f_p(i,0) +f_vis(i,0); //+f_grav(i);
    f(i,2)= -f_p(i,2) +f_vis(i,2)-g*rho(i);//*rho(i);
    //}

}
                       
                       
void sph::step()//MatrixX2d& u,MatrixX2d& x,double dt,MatrixXi& neighbours,VectorXd& p,VectorXd& rho,MatrixX2d& f)
{
    for (int i=0; i<N; ++i) {
        searchNeighbour(i);//,x,neighbours);
    }
    for (int i=0; i<N; ++i) {
        computeDensity(i);//rho,i,x,neighbours);
        computePressure(i);
       // double tempp=k*(rho(i)-rho_0);
       /* if(tempp>0)
        {*/
      //  p(i)=tempp;// computePressure(i);
      /*  }else {p(i)=0.;}
        */
    }
  //  std::cout <<p.size();
    for (int i=0; i<N; ++i)
    {
        computeForce(i);//,f,p,u,x,rho,neighbours);
    }
    for (int i=0; i<N; ++i) {
        u(i,1) += dt*f(i,1)/rho(i);
        u(i,0) += dt*f(i,0)/rho(i);
        u(i,2) += dt*f(i,2)/rho(i);
        
        if((x(i,1)+ dt*u(i,1))<0.) {     //collision detection with boundaries
          x(i,1)=0.;
            u(i,1)=-u(i,1);
        }else if((x(i,1)+dt*u(i,1))>=boundary)
        {
            x(i,1)=boundary;
            u(i,1)=-u(i,1);
        }else
        {
            x(i,1) += dt*u(i,1);
        }
        
        if((x(i,0)+ dt*u(i,0))<=0) {
            x(i,0)=0.;
            u(i,0)=-u(i,0);;
        }else if((x(i,0)+dt*u(i,0))>=boundary)
        {
            x(i,0)=boundary;
            u(i,0)=-u(i,0);
        }else
        {
            x(i,0) += dt*u(i,0);
        }
        
        if((x(i,2)+ dt*u(i,2))<=0.) {
            x(i,2)=0.;
            u(i,2)=-u(i,2);
        }else if((x(i,2)+dt*u(i,2))>=boundary)
        {
            x(i,2)=boundary;
            u(i,2)=-u(i,2);
        }else
        {
            x(i,2) += dt*u(i,2);
        }
        //x(i,0) += dt*u(i,0);
        //x(i,1) += dt*u(i,1);
    //x(i,2) += dt*u(i,2);

      //  std::cout <<"x:  "<<x(i,0)<<", "<<x(i,1)<<", "<<x(i,2)<<" )"<<"        f:  "<<f(i,0)<<", "<<f(i,1)<<", "<<f(i,2)<<"      p:  "<<p(i)<<"\n";
       // std::cout<<p(i)<<"   "<<rho(i)<<"\n";
      //  std::cout <<"( "<<neighbours(i,1)<<", " <<neighbours(i,2)<<", " <<neighbours(i,3)<<", " <<neighbours(i,4)<<", " <<neighbours(i,5)<<", " <<neighbours(i,6)<<", " <<neighbours(i,7)<<" )"<<std::endl;
    }
    
    
    
}
void sph::render()//MatrixX2d& x)
{
   // glColor3f(1,10,0.5);
    for(int i=0; i<N; i++)
    {
        glPushMatrix();
        {
            glTranslated(x(i,0), x(i,1), x(i,2));
            glutSolidSphere(0.3, 50, 50);
        }
        glPopMatrix();
    }
}


/*int main()
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
  //  std::cout <<"\n\n ";
    for (int i=0; i<N; ++i) {
            std::cout <<"( "<< x(i,0)<<", " <<x(i,1)<<" )"<<std::endl;
    }
    for (int i=0; i<N; ++i) {
        std::cout <<"( "<<neighbours(i,1)<<", " <<neighbours(i,2)<<", " <<neighbours(i,3)<<", " <<neighbours(i,4)<<", " <<neighbours(i,5)<<", " <<neighbours(i,6)<<", " <<neighbours(i,7)<<" )"<<std::endl;
    }
    return 0;
}*/
