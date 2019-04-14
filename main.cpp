#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <dirent.h>
#include <sstream>
#include "eigen3/Eigen/Dense"
#include "image.h"
#include "eigen3/Eigen/Eigenvalues" // used to decompose matricies

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::EigenSolver;
using namespace std;

//Created by Ben Penwell and Adam Landis
//Pattern Recognition, Project 3
//April 9, 2019

const int ID_LENGTH = 5;
const int DATE_LENGTH = 6;
const int SET_LENGTH = 2;

class Image{
public:
	void ExtractEigenVectors(MatrixXd m, int dimension1, int dimension2);

	void SetIdentifier(int input) {identifier = input; };
	void SetDateTaken(int input) {dateTaken = input; };
	void SetDataset(string input) {dataset = input; };
	void SetWearingGlasses(bool input) {wearingGlasses = input; };
	void SetFileName(string input) {fileName = input; };
	int GetIdentifier() {return identifier; };
	int GetDateTaken() {return dateTaken; };
	string GetFileName() {return fileName; };
	string GetDataset() {return dataset; };
	bool GetWearingGlasses() {return wearingGlasses; };


private:

	string fileName;
	int identifier; //nnnnn
	int dateTaken; //yymmdd
	string dataset;
	bool wearingGlasses;

	MatrixXcd eigenVectors;
	MatrixXcd eigenValues;
};

void readImageHeader(char[], int&, int&, int&, bool&);
void readImage(char[], ImageType&);
void writeImage(char[], ImageType&);
void PrintMatrix(MatrixXd m, int dimension1, int dimension2);
vector<Image> PopulateImages(string directory);
const int DIM = 3;

int main()
{
	MatrixXd g(DIM, DIM);
	string trainingDataset = "Faces_FA_FB/fa_H";

	g << 1, 0, 0, 2, 1, 1, 0, 0, 1;
	//printMatrix(g, DIM, DIM);

	vector<Image> ImageVector;
	//testImage.ExtractEigenVectors(g, DIM, DIM);

	if(1)
	{
		ImageVector = PopulateImages(trainingDataset);
	}
}

void PrintMatrix(MatrixXd m, int dimension1, int dimension2)
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

vector<Image> PopulateImages(string directory)
{
	vector<Image> returnVector;
	DIR *dir;
	struct dirent *ent;
		if ((dir = opendir (directory.c_str())) != NULL) 
		{
		  /* print all the files and directories within directory */
			while ((ent = readdir (dir)) != NULL) 
			{
				string temp, fileName = ent->d_name;
				//cout << "fileName" << fileName << endl;
				int imageID, imageDate, i;
				string imageDataset;
				bool imageGlasses;
				for (i=0; i < ID_LENGTH; ++i)
				{
				    temp += fileName[i];
				}
				imageID = atoi(temp.c_str());
				i++;

				temp = "";
				for (; i < ID_LENGTH+DATE_LENGTH+1; ++i)
				{
					temp += fileName[i];
				}
				imageDate = atoi(temp.c_str());
				i++;

				temp = "";
				for (; i < ID_LENGTH+DATE_LENGTH+SET_LENGTH+2; ++i)
				{
					temp += fileName[i];
				}
				imageDataset = temp;

				//cout << "fileName[i] " << fileName[i] << endl; 
				if(fileName[i] == '.')
				{
					imageGlasses = false;
				}
				else
				{
					imageGlasses = true;
				}

				Image newImage;
				newImage.SetIdentifier(imageID);
				newImage.SetDateTaken(imageDate);
				newImage.SetFileName(ent->d_name);
				newImage.SetDataset(imageDataset);
				newImage.SetWearingGlasses(imageGlasses);
				//cout << newImage.GetIdentifier() << " " << newImage.GetDateTaken() << " " << newImage.GetFileName() << " " << newImage.GetDataset() << " " << newImage.GetWearingGlasses() << endl;
				
				string currentFile = directory + '/' + ent->d_name;
				cout << "ent->d_name: " << ent->d_name << endl;
				cout << "currentFile: " << currentFile << endl;
				if(fileName != "." && fileName != "..")
				{
					int k, j, M, N, Q;
				 	bool type;
					int val;
					
					readImageHeader((char*) currentFile.c_str(), N, M, Q, type);

					 // allocate memory for the image array
					ImageType image(N, M, Q);
	 				readImage((char*) currentFile.c_str(), image);

	 				MatrixXd imageMatrix(N,M);
					// threshold image
					for(k=0; k<N; k++)
						for(j=0; j<M; j++) {
			            image.getPixelVal(k, j, val);
			            imageMatrix(k,j) = val;
						//cout << val << endl;
					}
					//newImage.ExtractEigenVectors(imageMatrix, N, M);

					//prints pixel values to output file. I am using this to validate the image has been correctly extracted
					ofstream fout;
					fout.open("output.txt", std::ofstream::trunc);
					for (int i = 0; i < N; ++i)
					{
						for (int j = 0; j < M; ++j)
						{
							fout << imageMatrix(i,j) << " ";
						}
						fout << endl;
					}
					fout.close();

					//THIS DOESN'T WORK CORRECTLY
					writeImage((char*)"output.pgm", image);
				}
			}
			closedir (dir);

			return returnVector;
		} 
		else 
		{
		 	/* could not open directory */
			return returnVector;
		}
}

void Image::ExtractEigenVectors(MatrixXd m, int dimension1, int dimension2)
{
	MatrixXd m_tm = m.transpose()*m;

	cout << "m.transpose*m:" << endl;
	PrintMatrix(m_tm, m_tm.rows(), m_tm.cols());

	EigenSolver<MatrixXd> EigenSolver;
	EigenSolver.compute(m_tm,true); //Initializes eigensolver with something to de-compose
	eigenVectors = EigenSolver.eigenvectors();
	cout << endl << "EigenVector Matrix: " << endl << eigenVectors << endl;

	eigenValues = EigenSolver.eigenvalues();
	cout << endl << "EigenValues Matrix: " << endl << eigenValues << endl;
}
