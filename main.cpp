#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Eigenvalues" // used to decompose matricies

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::EigenSolver;
using namespace std;

//Created by Ben Penwell and Adam Landis
//Pattern Recognition, Project 3
//April 9, 2019

class Image{
public:
	void ExtractEigenValues(MatrixXd m, int dimension1, int dimension2);
private:
	vector<float> eigenVectors;
};

void printMatrix(MatrixXd m, int dimension1, int dimension2);
const int DIM = 3;

int main()
{
	MatrixXd g(DIM, DIM);

	g << 1, 0, 0, 2, 1, 1, 0, 0, 1;
	printMatrix(g, DIM, DIM);

	Image testImage;
	testImage.ExtractEigenValues(g, DIM, DIM);
}

void printMatrix(MatrixXd m, int dimension1, int dimension2)
{
	for (int i = 0; i < dimension1; ++i)
	{
		for (int j = 0; j < dimension2; ++j)
		{
			cout << m(i,j) << " ";
		}
		cout << endl;
	}
}

void Image::ExtractEigenValues(MatrixXd m, int dimension1, int dimension2)
{
	MatrixXd mm_t = m*m.transpose();
	
	//cout << "m:" << endl;
	//printMatrix(m, m.rows(), m.cols());
	
	//cout << "m.transpose:" << endl;
	//printMatrix(m.transpose(), m.rows(), m.cols());
	
	cout << "m*m.transpose:" << endl;
	printMatrix(mm_t, mm_t.rows(), mm_t.cols());

	EigenSolver<MatrixXd> EigenSolver;
	EigenSolver.compute(mm_t,false); //Initializes eigensolver with something to de-compose
	printMatrix(EigenSolver.eigenvectors(), EigenSolver.eigenvectors().rows(), EigenSolver.eigenvectors().cols());

}