#include "ComJacobian.h"

	template <typename t_matrix>
t_matrix PseudoInverse(const t_matrix& m, const double &tolerance=1.e-6)
{
	using namespace Eigen;
	typedef JacobiSVD<t_matrix> TSVD;
	unsigned int svd_opt(ComputeThinU | ComputeThinV);
	if(m.RowsAtCompileTime!=Dynamic || m.ColsAtCompileTime!=Dynamic)
		svd_opt= ComputeFullU | ComputeFullV;
	TSVD svd(m, svd_opt);
	const typename TSVD::SingularValuesType &sigma(svd.singularValues());
	typename TSVD::SingularValuesType sigma_inv(sigma.size());
	for(long i=0; i<sigma.size(); ++i)
	{
		if(sigma(i) > tolerance)
			sigma_inv(i)= 1.0/sigma(i);
		else
			sigma_inv(i)= 0.0;
	}
	return svd.matrixV()*sigma_inv.asDiagonal()*svd.matrixU().transpose();
}

MatrixXd CoMJacobian::totalCoMJacobian(int targetlink)
{
	MatrixXd J;
	MatrixXd J_ref = MatrixXd::Zero(3,6);
	//total = 0.0;
	total = calcTotalMass(ulink,WAIST);
	vector<int> idy = kine->FindRoute(targetlink);
	const int ksize = idy.size();
	for(int to=0; to<ksize; to++){
		vector<int> idx = kine->FindRoute(idy[to]);

		const int jsize = idx.size();

		J.resize(3, jsize);

		J = calcCoMJacobian(bcf, idx);

		J_ref = (ulink[to].m*J) + J_ref;

		//total += ulink[to].m;
	}
	return J_ref/total;
}

MatrixXd CoMJacobian::calcCoMJacobian(int from, vector<int> idx)
{
	size_t jsize = idx.size();
	Matrix<double,3,1> target = ulink[idx.back()].c;
	//Matrix<double,3,1> target = ulink[idx.back()].c - ulink[from].c;
	Matrix<double,3,6> J = MatrixXd::Zero(3,6);

	for(size_t i=0; i<jsize; i++)
	{
		int j = idx[i];
		Matrix<double,3,1> a = ulink[j].R * ulink[j].a;
		//Matrix<double,3,1> b = a.cross(target - (ulink[j].c - ulink[from].c));
		Matrix<double,3,1> b = a.cross(target - (ulink[j].c));

		J(0,i) = b(0); J(1,i) = b(1); J(2,i) = b(2);
	}
	return J;
}

void CoMJacobian::fixedJacobian(vector<int> idx)
{
	const int jsize = idx.size();

	J_fix.resize(6, jsize);
	J_fix = calcfixedJacobian(bcf, idx);

	J_fix_v.resize(3, jsize);
	J_fix_w.resize(3, jsize);

	for(size_t i=0; i<jsize; i++)
	{
		J_fix_v(0,i) = J_fix(0,i); J_fix_v(1,i) = J_fix(1,i); J_fix_v(2,i) = J_fix(2,i);
		J_fix_w(0,i) = J_fix(3,i); J_fix_w(1,i) = J_fix(4,i); J_fix_w(2,i) = J_fix(5,i);
	}
}

Matrix<double,3,1> CoMJacobian::calcCoMPosition(int from)
{
	return -ulink[from].p + ulink[WAIST].c;
}

Matrix<double,3,1> CoMJacobian::calcfixedlinkPosition(int from, int end)
{
	return -ulink[from].p + ulink[end].p;
}

Matrix<double,3,3> CoMJacobian::skewsymmetric(Matrix<double,3,1> a)
{
	Matrix<double,3,3> outer;
	outer <<	0, -a(2), a(1),
				a(2), 0, -a(0),
				-a(1), a(0), 0;

	return outer;
}

MatrixXd CoMJacobian::calcfixedJacobian(int from, vector<int> idx)
{
	size_t jsize = idx.size();
	Matrix<double,3,1> target = ulink[idx.back()].p - ulink[from].p;
	Matrix<double,6,12> J = MatrixXd::Zero(6,12);

	for(size_t i=0; i<jsize; i++)
	{
		int j = idx[i];
		Matrix<double,3,1> a = ulink[j].R * ulink[j].a;
		Matrix<double,3,1> b = a.cross(target - (ulink[j].p - ulink[from].p));
		J(0,i) = b(0); J(1,i) = b(1); J(2,i) = b(2);
		J(3,i) = a(0); J(4,i) = a(1); J(5,i) = a(2);
	}
	return J;
}

Matrix<double,3,6> CoMJacobian::resultCJ(int targetlink, Matrix<double,3,1> eff)
{
	Matrix<double,3,6> J;
	J = totalCoMJacobian(targetlink);
	vector<int> idx = kine->FindRoute(targetlink);
	fixedJacobian(idx);

	//p_g = calcCoMPosition(bcf);
	p_g = calcCoM(ulink);
	//cout << p_g(1)*10 *-100;
	p_f = calcfixedlinkPosition(bcf,targetlink);
	R_0 <<	1,0,0,
			0,1,0,
			0,0,1;

	K = skewsymmetric(p_g - p_f);

	return R_0 * (J - J_fix_v + (K * J_fix_w));
}

void CoMJacobian::calcInverse(int targetlink, Link target, Matrix<double,3,1> eff)
{
	size_t jsize;
	vector<int> idx;
	Matrix<double,6,1> err;
	MatrixXd J(3,6);

	idx = kine->FindRoute(targetlink);
	Matrix<double,6,3> ref_J = Matrix<double,6,3>::Zero();
	Matrix<double,6,1> dp = Matrix<double,6,1>::Zero();

	J = resultCJ(targetlink, eff);
	cout << "ComJacobian\n" << J << endl;
	ref_J = PseudoInverse(J);

	dp = ref_J * eff;

	for(size_t nn=0;nn<idx.size();nn++){
		int j = idx[nn];
		kine->ulink[j].q += dp(nn)*0.01;
	}
}

