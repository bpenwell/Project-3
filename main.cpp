#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Eigenvalues" // used to decompose matricies

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::EigenSolver;
using namespace std;

//Created by Ben Penwell and Adam Landis
//Pattern Recognition, Project 3
//April 9, 2019


class Image{
public:
	void ExtractEigenVectors(MatrixXd m, int dimension1, int dimension2);
	void SortEigenVectors();
	void ComputeColumnValue(VectorXd c);
private:
	MatrixXcd eigenVectors;
	MatrixXcd eigenValues;
	vector<VectorXd> sortedEigenVectors;

};

void printMatrix(MatrixXd m, int dimension1, int dimension2);
const int DIM = 3;

int main()
{
	MatrixXd g(DIM, DIM);

	g << 1, 0, 0, 2, 1, 1, 0, 0, 1;
	//printMatrix(g, DIM, DIM);

	Image testImage;
	testImage.ExtractEigenVectors(g, DIM, DIM);
	//testImage.SortEigenVectors();
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

void Image::ComputeColumnValue(VectorXd c)
{
	for (int i = 0; i < c.size(); ++i)
	{
		/* code */
	}
}


void Image::SortEigenVectors()
{
	cout << "Rows: " << eigenVectors.rows() << endl;
	for (int i = 0; i < eigenVectors.cols(); ++i)
	{
		complex<double> doubleBarOperation=0.0; //sqrt(x1^2 + .... + xn^2)
		for (int j = 0; j < eigenVectors.rows(); ++j)
		{
			doubleBarOperation += (eigenVectors(j,i)*eigenVectors(j,i));
			cout << doubleBarOperation << " ";
		}
		complex<double> value = sqrt(doubleBarOperation);
		cout << endl << "Col value: " << i << endl;
		cout << value << endl;
		if(sortedEigenVectors.size() == 0)
		{
			VectorXd column(eigenVectors.rows());
		}
		else
		{

		}
		for (int i = 0; i < sortedEigenVectors.size(); ++i)
		{
			/* code */
		}
	}
}

void Image::ExtractEigenVectors(MatrixXd m, int dimension1, int dimension2)
{
	MatrixXd m_tm = m.transpose()*m;

	cout << "m.transpose*m:" << endl;
	printMatrix(m_tm, m_tm.rows(), m_tm.cols());

	EigenSolver<MatrixXd> EigenSolver;
	EigenSolver.compute(m_tm,true); //Initializes eigensolver with something to de-compose
	eigenVectors = EigenSolver.eigenvectors();
	cout << endl << "EigenVector Matrix: " << endl << eigenVectors << endl;

	eigenValues = EigenSolver.eigenvalues();
	cout << endl << "EigenValues Matrix: " << endl << eigenValues << endl;
}
