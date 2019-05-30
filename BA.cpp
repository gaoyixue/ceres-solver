#include<cmath>
#include<cstdio>
#include<iostream>
#include<ceres/ceres.h>
#include<ceres/rotation.h>
using namespace std;
using namespace ceres;

class BALproblem
{
public:
	~BALproblem() {
		delete[] observations;
		delete[] parameters;
		delete[] camera_index;
		delete[] point_index;
	}

	int get_num_camera(){return num_cameras_;}
	int get_num_point() { return num_points_; }
	int get_num_observations() { return num_observations_; }
	int get_num_parameters() { return num_parameters_; }
	double get_observations(int i) { return observations[i]; }


	double *get_camera(int i)
	{
		return parameters + 9 * camera_index[i];
	}

	double *get_point(int i)
	{
		return parameters + 9 * num_cameras_ + 3 * point_index[i];
	}

	bool loadFile(char *filename)
	{
		FILE* fd=fopen(filename, "r");
		if (fd==NULL)
		{
			printf("open file is error..\n");
			return false;
		}
		fscanf(fd, "%d", &num_cameras_);
		fscanf(fd, "%d", &num_points_);
		fscanf(fd, "%d", &num_observations_);
		num_parameters_ = 9 * num_cameras_ + 3 * num_points_;
		point_index = new int[num_observations_];
		camera_index = new int[num_observations_];
		parameters = new double[num_parameters_];
		observations = new double[2 * num_observations_];

		for (int i = 0; i < num_observations_; ++i)
		{
			fscanf(fd, "%d", camera_index + i);//写成camera_index[i]会引发异常。可以写成camera_index+i或者&camera_index[i]
			fscanf(fd, "%d", point_index + i);
			fscanf(fd, "%lf", observations + 2 * i);
			fscanf(fd, "%lf", observations + 2 * i + 1);
			//cout << camera_index[i] << "\n";
			//cout << point_index[i] << "\n";
			//cout << observations[2*i] << "\n";
			//cout << observations[2*i+1] << "\n";
		}

		for (int i = 0; i < num_parameters_; ++i)
		{
			fscanf(fd, "%lf", parameters+i);
		}
		return true;
	}

private:
	int num_cameras_;
	int num_points_;
	int num_observations_;
	int num_parameters_;

	double *observations;
	double *parameters;
	int *camera_index;
	int *point_index;
};

struct costFunc
{
	costFunc(double x,double y):observed_x(x), observed_y(y){}

	template<typename T>
	bool operator()(const T* const camera,const T* const point,T* residuals)const
	{
		//point指的是三维点，需要经过坐标变换变到图像坐标系
		//camera共有9个参数，表示R的Rodrigues矢量，一个T，f,k1,k2
		//p=R*point+T,需要将Rodrigues矢量变为旋转矩阵R
		//xx=-p.x/p.z,yy=-p.y/p.z
		//xx'=f*r*xx,yy'=f*r*yy.其中r=1+k1*r2+k2*(r2*r2),r2=xx*xx+yy*yy
		T p[3];
		ceres::AngleAxisRotatePoint(camera, point, p);
		p[0] = p[0] + camera[3];
		p[1] = p[1] + camera[4];
		p[2] = p[2] + camera[5];//已将三维点point转换到相机坐标系下
		T xp = -p[0] / p[2];
		T yp = -p[1] / p[2];
		T r2 = xp*xp + yp*yp;

		T r = 1.0 + r2  * camera[7] + camera[8] * r2*r2;

		// Compute final projected point position.

		T predicted_x = camera[6] * r * xp;

		T predicted_y = camera[6] * r * yp;

		// The error is the difference between the predicted and observed position.

		residuals[0] = predicted_x - observed_x;

		residuals[1] = predicted_y - observed_y;
		return true;
	}
	double observed_x, observed_y;//传入的观察图像坐标
};

int main(int argc, char**argv)
{
	google::InitGoogleLogging(argv[0]);

	BALproblem balproblem;
	if (balproblem.loadFile("E:\\浏览器下载文件\\problem-93-61203-pre.txt"))
	{
		cout << "load file success\n";
	}

	Problem problem;
	for (int i = 0; i < balproblem.get_num_observations(); ++i)
	{
		CostFunction* costfunc = new AutoDiffCostFunction<costFunc, 2, 9, 3>
			(new costFunc(balproblem.get_observations(2 * i), balproblem.get_observations(2 * i + 1)));
		problem.AddResidualBlock(costfunc, NULL, balproblem.get_camera(i), balproblem.get_point(i));
	}

	Solver::Options options;
	options.linear_solver_type = DENSE_SCHUR;
	options.minimizer_progress_to_stdout = true;
	options.num_threads = 4;
	Solver::Summary summary;
	Solve(options, &problem, &summary);
	std::cout << summary.BriefReport() << "\n";
	cout << "优化后的相机：" << *balproblem.get_camera(0) << " " << *(balproblem.get_camera(0) + 1) << " " << *(balproblem.get_camera(0) + 2);
	return 0;
}