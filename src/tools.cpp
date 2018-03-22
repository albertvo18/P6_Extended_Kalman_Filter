#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;


// Constructor
Tools::Tools() {}

// Destructor
Tools::~Tools() {}

// From Eigen Library, helper method to CalculateRMSE
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    if(estimations.size() == 0){
      cout << "VALUE ERROR:  CalculateRMSE():  Estimation Vector is 0" << endl;
      return rmse;
    }

    if(ground_truth.size() == 0){
      cout << "VALUE ERROR:  CalculateRMSE():  Ground Truth Vector is 0" << endl;
      return rmse;
    }

    unsigned int estimation_size = estimations.size();
    if(estimation_size != ground_truth.size()){
      cout << "VALUE ERROR:  CalculateRMSE():  Estimation Size not equal to Ground Truth Vector " << endl;
      return rmse;
    }

    for(unsigned int i=0; i < estimations.size(); ++i){
      VectorXd diff = estimations[i] - ground_truth[i];
      diff = diff.array()*diff.array();
      rmse += diff;
    }

    rmse = rmse / estimation_size;
    rmse = rmse.array().sqrt();
    return rmse;
}

// From Eigen Library, helper method to CalculateJacobian Matrix
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);

  if ( x_state.size() != 4 ) {
    cout << "VALUE ERROR: CalculateJacobian(): state vector is not equal to 4." << endl;
    return Hj;
  }
	//parameters
	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

	//compute terms
	double c1 = px*px+py*py;
	double c2 = sqrt(c1);
	double c3 = (c1*c2);

	//check division by zero
	if(fabs(c1) < 0.0001){
		cout << "VALUE ERROR: CalculateJacobian(): Division by Zero" << endl;
		return Hj;
	}

	//compute Jacobian 
	Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;
}
