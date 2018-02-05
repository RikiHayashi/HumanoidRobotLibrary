#include "ComJacobian.h"
#include "../../util/func.h"

#include <iostream>
using namespace std;

double xv[3];
int getf;

int main(int argc, char* argv[]){
	Link ulink[JOINT_NUM];
	Kinematics kine(ulink);
	SetJointInfo(ulink);
	CoMJacobian CJ(ulink);

	for(int i=0; i<JOINT_NUM;i++)
		ulink[i].q = 0.0;
	
	ulink[RLEG_JOINT2].q = deg2rad(15);
	ulink[RLEG_JOINT3].q = deg2rad(-30);
	ulink[RLEG_JOINT4].q = deg2rad(15);
	ulink[LLEG_JOINT2].q = deg2rad(15);
	ulink[LLEG_JOINT3].q = deg2rad(-30);
	ulink[LLEG_JOINT4].q = deg2rad(15);

	kine.calcForwardKinematics(WAIST);
	Link RLEG_LINK = ulink[RLEG_JOINT5];
	Link LLEG_LINK = ulink[LLEG_JOINT5];

	FILE *fp;

	fp = fopen("out.csv","r");

	for(int count=0;;count++){
		getf = fgetc(fp);
		if(getf == EOF){
			fclose(fp);
			break;
		}

		fscanf(fp, "%lf,%lf", &xv[1], &xv[2]);

		Matrix<double,3,1> eff;
		eff << 0.0, 0.0, xv[2];
		
		//calculation Forward Kinematics
		CJ.kine->calcForwardKinematics(WAIST);

		//calculation Inverce Kinematics
		Link RLink;
		Link LLink;

		//calculation ComJacobian
		CJ.calcInverse(RLEG_JOINT5, RLink, eff);
		CJ.calcInverse(LLEG_JOINT5, LLink, eff);

		CJ.kine->calcForwardKinematics(WAIST);
	}
}
