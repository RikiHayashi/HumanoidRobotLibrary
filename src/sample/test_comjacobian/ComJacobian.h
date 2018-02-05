#ifndef _COMJACOBIAN_H_
#define _COMJACOBIAN_H_

#include<iostream>
#include<vector>
#include<Eigen/Core>

#include "../../kinematics/Jacobian.h"
#include "../../kinematics/Kinematics.h"
#include "../../kinematics/Link.h"
#include <Eigen/SVD>

using namespace std;

class CoMJacobian
{
public:
	struct Link *ulink;
	Kinematics *kine;
	CoMJacobian(Link *ulink)
	{
		this->ulink = ulink;
		kine = new Kinematics(ulink);
	}
public:
	Matrix<double,3,1> calcCoMPosition(int to);
	Matrix<double,3,1> calcfixedlinkPosition(int to, int end);

	MatrixXd calcfixedJacobian(int to, vector<int> idx);
	MatrixXd calcCoMJacobian(int to, vector<int> idx);

	Matrix<double,3,3> skewsymmetric(Matrix<double,3,1> a);

	MatrixXd totalCoMJacobian(int targetlink);
	void fixedJacobian(vector<int> idx);
	Matrix<double,3,6> resultCJ(int targetlink, Matrix<double,3,1> err);
	void calcInverse(int to, Link target,Matrix<double,3,1> eff);
public:
	int bcf = WAIST; //Body Center Frame
	float total = 0;
	MatrixXd J_fix;
	//MatrixXd J_ref = MatrixXd::Zero(3,6);

	MatrixXd J_fix_v;
	MatrixXd J_fix_w;

	Matrix<double,3,1> p_g;
	Matrix<double,3,1> p_f;
	Matrix<double,3,3> R_0;

	MatrixXd comJ;
	Matrix<double,3,3> K;
};

#endif
