/*#include<cstdio>
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<ceres/ceres.h>
#include<gflags/gflags.h>
#include<glog/logging.h>
using namespace std;
using namespace ceres;

class DistanceFromCircleCost {
public:
	DistanceFromCircleCost(double xx,double yy,double zz,double nx,double ny,double nz):xx_(xx),yy_(yy),zz_(zz),nx_(nx),ny_(ny),nz_(nz){}
	template <typename T>
	bool operator()(const T* const x,const T* const y,const T* const z, 
		const T* const nx, const T* const ny, const T* const nz,
		const T* const m,T* residual)const
	{
		T deltax = T(xx_) - x[0];
		T deltay = T(yy_) - y[0];
		T deltaz = T(zz_) - z[0];
		T r = m[0] * m[0];
		residual[0] = r*r - deltax*deltax - deltay*deltay - deltaz*deltaz;
		residual[1] = deltax*nx[0] + deltay*ny[0] + deltaz*nz[0];
		residual[2] = nx[0] * nx[0] + ny[0] * ny[0] + nz[0] * nz[0] - T(nx_*nx_ + ny_*ny_ + nz_*nz_);
		return true;
	}
private:
	double xx_, yy_,zz_,nx_,ny_,nz_;
};

int main(int argc,char** argv)
{
	google::InitGoogleLogging(argv[0]);

	double x, y, z, nx, ny, nz, r;

	string filename="C:\\Users\\gaoyixue\\Desktop\\1.txt";
	ifstream fin(filename.c_str());
	string line;
	getline(fin, line);
	char *pend;

	x = strtod(line.c_str(), &pend);
	y = strtod(pend, &pend);
	z = strtod(pend, &pend);
	r = strtod(pend, &pend);
	nx = strtod(pend, &pend);
	ny = strtod(pend, &pend);
	nz = strtod(pend, nullptr);

	fprintf(stderr, "got x,y,z,r,nx,ny,nz %lg %lg %lg %lg %lg %lg %lg\n", x, y, z, r, nx, ny, nz);

	double initial_x = x;
	double initial_y = y;
	double initial_z = z;
	double initial_r = r;
	double initial_nx = nx;
	double initial_ny = ny;
	double initial_nz = nz;

	double m = ceres::sqrt(r);
	Problem problem;
	double xx, yy, zz;
	while (getline(fin,line))
	{
		xx = strtod(line.c_str(), &pend);
		yy = strtod(pend, &pend);
		zz = strtod(pend, nullptr);
		//std::cout << "got (" << xx << "," << yy << ")\n";
		CostFunction *costfunc = new AutoDiffCostFunction<DistanceFromCircleCost, 3, 1, 1, 1, 1, 1, 1, 1>(new DistanceFromCircleCost(xx, yy, zz,nx,ny,nz));
		problem.AddResidualBlock(costfunc, new CauchyLoss(0.5), &x, &y, &z, &nx, &ny, &nz, &m);
	}

	//std::cout << "got " << numpoints << " points.\n";
	

	Solver::Options options;
	options.linear_solver_type = DENSE_QR;
	options.max_num_iterations = 500;
	Solver::Summary summary;

	Solve(options, &problem, &summary);

	r = m*m;
	std::cout << summary.BriefReport() << "\n";
	std::cout << "x:" << initial_x << "->" << x << "\n";
	std::cout << "y:" << initial_y << "->" << y << "\n";
	std::cout << "z:" << initial_z << "->" << z << "\n";
	std::cout << "r:" << initial_r << "->" << r << "\n";
	std::cout << "nx:" << initial_nx << "->" << nx << "\n";
	std::cout << "ny:" << initial_ny << "->" << ny << "\n";
	std::cout << "nz:" << initial_nz << "->" << nz << "\n";
	return 0;
}*/
