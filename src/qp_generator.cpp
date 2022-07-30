#include "qp_generator.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;    
using namespace Eigen;

Vector3d startVel;
Vector3d startAcc;
#define inf 1>>30

int m=3; // number of segments in the trajectory
vector<double> TrajectoryGenerator::getCost() { return qp_cost; }

TrajectoryGenerator::TrajectoryGenerator(){}

TrajectoryGenerator::~TrajectoryGenerator(){}

Eigen::MatrixXd TrajectoryGenerator::PolyQPGeneration(
            const Eigen::MatrixXd &Path,
            const Eigen::Vector3d &Vel,
            const Eigen::Vector3d &Acc,
            const Eigen::VectorXd &Time,
            const int &type) // 1 for using minumum snap for intial; 2 for using straight line for initial
{     
      /*   Get initial trajectory which is a straight line( 0 zero end velocity and acceleration ) minimum snap trajectory or truly minumum snap trajectory  */

      startVel = Vel; //[0; 0; 0]
      startAcc = Acc; //[0; 0; 0]
      m = Time.size(); //轨迹段数
      MatrixXd PolyCoeff(m, 3 * 6);
      VectorXd Px(6 * m), Py(6 * m), Pz(6 * m);

      int num_f, num_p; // number of fixed and free variables
      int num_d;        // number of all segments' derivatives
      const static auto Factorial = [](int x){
          int fac = 1;

          for(int i = x; i > 0; i--)
              fac = fac * i;
            
          return fac;
      };

      /*   Produce Mapping Matrix A to the entire trajectory. 将系数c0,...,c5映射为端点处的各阶导数，A*η = Dx  */
      MatrixXd Ab;
      MatrixXd A = MatrixXd::Zero(m * 6, m * 6);

      for(int k = 0; k < m; k++){
          Ab = Eigen::MatrixXd::Zero(6, 6);
          for(int i = 0; i < 3; i++){
              Ab(2 * i, i) = Factorial(i); //c0 + c1*t + ... + c5*t^5,第i阶导数在t=0时ci对应的系数
              for(int j = i; j < 6; j++)
                  Ab( 2 * i + 1, j ) = Factorial(j) / Factorial( j - i ) * pow( Time(k), j - i );//第i阶导数在t=Time(k)时cj对应的系数
          }

          A.block(k * 6, k * 6, 6, 6) = Ab;    
      }

      _A = A;

      /*   Produce the dereivatives in X, Y and Z axis directly.  */
      //等式约束，分别是端点处的位置，速度，加速度
      VectorXd Dx = VectorXd::Zero(m * 6);
      VectorXd Dy = VectorXd::Zero(m * 6);
      VectorXd Dz = VectorXd::Zero(m * 6);

      for(int k = 1; k < (m + 1); k ++ ){
          Dx((k-1)*6) = Path(k - 1, 0); Dx((k-1)*6 + 1) = Path(k, 0); //0-5：p0, p1, v0, v1, a0, a1
          Dy((k-1)*6) = Path(k - 1, 1); Dy((k-1)*6 + 1) = Path(k, 1); //同时中间路径点的速度和加速度为0
          Dz((k-1)*6) = Path(k - 1, 2); Dz((k-1)*6 + 1) = Path(k, 2); 
          
          if( k == 1){
              Dx((k-1)*6 + 2) = Vel(0);
              Dy((k-1)*6 + 2) = Vel(1); 
              Dz((k-1)*6 + 2) = Vel(2);

              Dx((k-1)*6 + 4) = Acc(0);
              Dy((k-1)*6 + 4) = Acc(1); 
              Dz((k-1)*6 + 4) = Acc(2);
          }
      }

      /*   Produce the Minimum Jerk cost function, the Hessian Matrix   */
      MatrixXd H = MatrixXd::Zero( m * 6, m * 6 );
      
      for(int k = 0; k < m; k ++){
          for(int i = 3; i < 6; i ++){
              for(int j = 3; j < 6; j ++){
                  H( k*6 + i, k*6 + j ) = i * (i - 1) * (i - 2) * j * (j - 1) * (j - 2) / (i + j - 5) * pow( Time(k), (i + j - 5) );
              }
          }
      }

      _Q = H; // Only minumum Jerk term is considered here inf the Hessian matrix

      /*   Produce Selection Matrix C'   */
      MatrixXd Ct; // The transpose of selection matrix C
      MatrixXd C;  // The selection matrix C

      if(type == 1){ // generating minumum Jerk curve
          num_f = 2 * m + 4; //3 + 3 + (m - 1) * 2 = 2m + 4
          num_p = 2 * m - 2; //(m - 1) * 2 = 2m - 2
          num_d = 6 * m; //每段轨迹的首尾端点处共6个状态
          Ct = MatrixXd::Zero(num_d, num_f + num_p); 
          //处理的是第1段轨迹
          Ct( 0, 0 ) = 1; Ct( 2, 1 ) = 1;         Ct( 4, 2 ) = 1; // stack the start point
          Ct( 1, 3 ) = 1; Ct( 3, 2 * m + 4 ) = 1; Ct( 5, 2 * m + 5 ) = 1; 
        
          //处理的是第m段轨迹
          Ct(6 * (m - 1) + 0, 2 * m + 0) = 1; 
          Ct(6 * (m - 1) + 1, 2 * m + 1) = 1; // Stack the end point
          Ct(6 * (m - 1) + 2, 4 * m + 0) = 1;
          Ct(6 * (m - 1) + 3, 2 * m + 2) = 1; // Stack the end point
          Ct(6 * (m - 1) + 4, 4 * m + 1) = 1;
          Ct(6 * (m - 1) + 5, 2 * m + 3) = 1; // Stack the end point

          for(int j = 2; j < m; j ++ ){ //j表示第j段轨迹
              Ct( 6 * (j - 1) + 0, 2 + 2 * (j - 1)         + 0 ) = 1;
              Ct( 6 * (j - 1) + 1, 2 + 2 * (j - 1)         + 1 ) = 1;
              Ct( 6 * (j - 1) + 2, 2 * m + 4 + 2 * (j - 2) + 0 ) = 1;
              Ct( 6 * (j - 1) + 3, 2 * m + 4 + 2 * (j - 1) + 0 ) = 1;
              Ct( 6 * (j - 1) + 4, 2 * m + 4 + 2 * (j - 2) + 1 ) = 1;
              Ct( 6 * (j - 1) + 5, 2 * m + 4 + 2 * (j - 1) + 1 ) = 1;
          }

          C = Ct.transpose();

          VectorXd Dx1 = C * Dx; //Dx1对应的是[df dp]
          VectorXd Dy1 = C * Dy;
          VectorXd Dz1 = C * Dz;

          MatrixXd Q; // The Hessian of the total cost function
          Q = H;      // Now only minumum jerk is used in the cost
          MatrixXd R   = C * A.transpose().inverse() * Q * A.inverse() * Ct;
          VectorXd Dxf(2 * m + 4), Dyf(2 * m + 4), Dzf(2 * m + 4);
          
          Dxf = Dx1.segment( 0, 2 * m + 4 );
          Dyf = Dy1.segment( 0, 2 * m + 4 );
          Dzf = Dz1.segment( 0, 2 * m + 4 );

          MatrixXd Rff(2 * m + 4, 2 * m + 4);
          MatrixXd Rfp(2 * m + 4, 2 * m - 2);
          MatrixXd Rpf(2 * m - 2, 2 * m + 4);
          MatrixXd Rpp(2 * m - 2, 2 * m - 2);

          Rff = R.block(0, 0, 2 * m + 4, 2 * m + 4);
          Rfp = R.block(0, 2 * m + 4, 2 * m + 4, 2 * m - 2);
          Rpf = R.block(2 * m + 4, 0,         2 * m - 2, 2 * m + 4);
          Rpp = R.block(2 * m + 4, 2 * m + 4, 2 * m - 2, 2 * m - 2);

          VectorXd Dxp(2 * m - 2), Dyp(2 * m - 2), Dzp(2 * m - 2);
          Dxp = - (Rpp.inverse() * Rfp.transpose()) * Dxf;
          Dyp = - (Rpp.inverse() * Rfp.transpose()) * Dyf;
          Dzp = - (Rpp.inverse() * Rfp.transpose()) * Dzf;

          Dx1.segment(2 * m + 4, 2 * m - 2) = Dxp;
          Dy1.segment(2 * m + 4, 2 * m - 2) = Dyp;
          Dz1.segment(2 * m + 4, 2 * m - 2) = Dzp;

          Px = (A.inverse() * Ct) * Dx1;
          Py = (A.inverse() * Ct) * Dy1;
          Pz = (A.inverse() * Ct) * Dz1;

          _Dx = Ct * Dx1; //_DX = [p_1s; v_1s; a_1s; p_1e; v_1e; a_1e; ...; p_ms; v_ms; a_ms; p_me; v_me; a_me]
          _Dy = Ct * Dy1;
          _Dz = Ct * Dz1;
          
          _Pxi = Px; _Pyi = Py; _Pzi = Pz;
      }

      else{ // generating piecewise straight line 
          num_f = 3 * m + 3; // All fixed 3 + 3 + (m - 1) * 3
          num_p = 0;         // No free derivatives
          num_d = 6 * m;
          
          //下面这块没用
          Ct = MatrixXd::Zero(num_d, num_d); 

          for(int j = 1; j < (m + 1); j ++ ){
              Ct( 6 * (j - 1) + 0, 6 * (j - 1) + 0 ) = 1;
              Ct( 6 * (j - 1) + 1, 6 * (j - 1) + 3 ) = 1;
              Ct( 6 * (j - 1) + 2, 6 * (j - 1) + 1 ) = 1;
              Ct( 6 * (j - 1) + 3, 6 * (j - 1) + 4 ) = 1;
              Ct( 6 * (j - 1) + 4, 6 * (j - 1) + 2 ) = 1;
              Ct( 6 * (j - 1) + 5, 6 * (j - 1) + 5 ) = 1;
          }
          C = Ct.transpose();
          //上面这块没用

          Px = A.inverse() * Dx;
          Py = A.inverse() * Dy;
          Pz = A.inverse() * Dz;

          _Dx = Dx;
          _Dy = Dy;
          _Dz = Dz;

          _Pxi = Px; _Pyi = Py; _Pzi = Pz;
      }
      
      for(int i = 0; i < m; i ++){
          PolyCoeff.block(i, 0,  1, 6) = Px.segment( i * 6, 6 ).transpose();
          PolyCoeff.block(i, 6,  1, 6) = Py.segment( i * 6, 6 ).transpose();
          PolyCoeff.block(i, 12, 1, 6) = Pz.segment( i * 6, 6 ).transpose();
      }

      ROS_WARN("[Generator] Unconstrained QP solved");
      return PolyCoeff; //PolyCoeff每行对应的是：[Px(0-5) Py(0-5) Pz(0-5)]
}  

void TrajectoryGenerator::StackOptiDep(){
      
      double num_f = 6;          // 3 + 3 : only start and target position has fixed derivatives   
      double num_p = 3 * m - 3;  // All other derivatives are free   
      double num_d = 6 * m;
            
      MatrixXd Ct;    
      Ct = MatrixXd::Zero(num_d, num_f + num_p); 
      Ct( 0, 0 ) = 1; Ct( 2, 1 ) = 1; Ct( 4, 2 ) = 1;  // Stack the start point
      Ct( 1, 6 ) = 1; Ct( 3, 7 ) = 1; Ct( 5, 8 ) = 1; 

      Ct(6 * (m - 1) + 0, 3 * m + 0 ) = 1; 
      Ct(6 * (m - 1) + 2, 3 * m + 1 ) = 1;
      Ct(6 * (m - 1) + 4, 3 * m + 2 ) = 1;

      Ct(6 * (m - 1) + 1, 3) = 1; // Stack the end point
      Ct(6 * (m - 1) + 3, 4) = 1;
      Ct(6 * (m - 1) + 5, 5) = 1;

      for(int j = 2; j < m; j ++ ){
          Ct( 6 * (j - 1) + 0, 6 + 3 * (j - 2) + 0 ) = 1;
          Ct( 6 * (j - 1) + 1, 6 + 3 * (j - 1) + 0 ) = 1;
          Ct( 6 * (j - 1) + 2, 6 + 3 * (j - 2) + 1 ) = 1;
          Ct( 6 * (j - 1) + 3, 6 + 3 * (j - 1) + 1 ) = 1;
          Ct( 6 * (j - 1) + 4, 6 + 3 * (j - 2) + 2 ) = 1;
          Ct( 6 * (j - 1) + 5, 6 + 3 * (j - 1) + 2 ) = 1;
      }

      _C = Ct.transpose();
      _L = _A.inverse() * Ct; //η = _L * [df dp]'

      MatrixXd B = _A.inverse();
      _R = _C * B.transpose() * _Q * (_A.inverse()) * Ct;

      _Rff.resize(6, 6);
      _Rfp.resize(6, 3 * m - 3);
      _Rpf.resize(3 * m - 3, 6);
      _Rpp.resize(3 * m - 3, 3 * m - 3);

      _Rff =  _R.block(0, 0, 6, 6);
      _Rfp =  _R.block(0, 6, 6, 3 * m - 3);
      _Rpf =  _R.block(6, 0, 3 * m - 3, 6);
      _Rpp =  _R.block(6, 6, 3 * m - 3, 3 * m - 3);
      ROS_WARN("[Generator] Stack finish");
}

std::pair< Eigen::MatrixXd, Eigen::MatrixXd > TrajectoryGenerator::getInitialD()
{  
      // Return Format: row(0) -> X | row(1) -> Y | row(2) -> Z

      VectorXd Dxf = VectorXd::Zero(6); //[ps; vs; as; pe; ve; ae]
      VectorXd Dyf = VectorXd::Zero(6);
      VectorXd Dzf = VectorXd::Zero(6);

      VectorXd Dxp = VectorXd::Zero(3 * m - 3);
      VectorXd Dyp = VectorXd::Zero(3 * m - 3);
      VectorXd Dzp = VectorXd::Zero(3 * m - 3);

      Dxf(0) = _Dx(0); Dxf(3) = _Dx( _Dx.size() - 5 ); 
      Dyf(0) = _Dy(0); Dyf(3) = _Dy( _Dy.size() - 5 ); 
      Dzf(0) = _Dz(0); Dzf(3) = _Dz( _Dz.size() - 5 ); 

      Dxf(1) = startVel(0);
      Dyf(1) = startVel(1);
      Dzf(1) = startVel(2);
      
      Dxf(2) = startAcc(0);
      Dyf(2) = startAcc(1);
      Dzf(2) = startAcc(2); //终点的速度和加速度为0

      for (int k = 1; k < m; k ++ ){
          for (int i = 0; i < 3; i ++ ){
                Dxp( (k - 1) * 3 + i) = _Dx((k - 1) * 6 + 2 * i + 1);
                Dyp( (k - 1) * 3 + i) = _Dy((k - 1) * 6 + 2 * i + 1);
                Dzp( (k - 1) * 3 + i) = _Dz((k - 1) * 6 + 2 * i + 1);
          }
      }

      MatrixXd Dp(3, Dxp.size()), Df(3, Dxf.size());
      Dp.row(0) = Dxp;
      Dp.row(1) = Dyp;
      Dp.row(2) = Dzp;

      Df.row(0) = Dxf;
      Df.row(1) = Dyf;
      Df.row(2) = Dzf;

      return make_pair(Dp, Df);
}

Eigen::MatrixXd TrajectoryGenerator::getA(){  return _A;  }
Eigen::MatrixXd TrajectoryGenerator::getQ(){  return _Q;  }
Eigen::MatrixXd TrajectoryGenerator::getC(){  return _C;  }
Eigen::MatrixXd TrajectoryGenerator::getL(){  return _L;  }
Eigen::MatrixXd TrajectoryGenerator::getR(){  return _R;  }

Eigen::MatrixXd TrajectoryGenerator::getRff(){  return _Rff;  }
Eigen::MatrixXd TrajectoryGenerator::getRpp(){  return _Rpp;  }
Eigen::MatrixXd TrajectoryGenerator::getRpf(){  return _Rpf;  }
Eigen::MatrixXd TrajectoryGenerator::getRfp(){  return _Rfp;  }