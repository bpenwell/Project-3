#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <dirent.h>
#include <sstream>
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Eigenvalues" // used to decompose matricies
#include "image.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::MatrixXf;
using Eigen::MatrixXcd;
using Eigen::MatrixXcf;
using Eigen::EigenSolver;
using namespace std;

//Created by Ben Penwell and Adam Landis
//Pattern Recognition, Project 3
//April 9, 2019

const int ID_LENGTH = 5;
const int DATE_LENGTH = 6;
const int SET_LENGTH = 2;
const int DIM = 3;
const int IMG_W = 48, IMG_H = 60, IMG_VEC_LEN = IMG_H * IMG_W;
const int NUM_SAMPLES = 1204;

class Image
{
public:
	void ExtractEigenVectors(MatrixXd m, int dimension1, int dimension2);

	// Setters
	void SetIdentifier(int input) { identifier = input; };
	void SetDateTaken(int input) { dateTaken = input; };
	void SetDataset(string input) { dataset = input; };
	void SetWearingGlasses(bool input) { wearingGlasses = input; };
	void SetFileName(string input) { fileName = input; };
	void setFaceVector(VectorXi vector) { faceVector = vector; };

	// Getters
	int GetIdentifier() { return identifier; };
	int GetDateTaken() { return dateTaken; };
	string GetFileName() { return fileName; };
	string GetDataset() { return dataset; };
	bool GetWearingGlasses() { return wearingGlasses; };
	VectorXi getFaceVector() const { return faceVector; };

private:
	string fileName;
	int identifier; //nnnnn
	int dateTaken; //yymmdd
	string dataset;
	bool wearingGlasses;

	VectorXi faceVector;
	MatrixXcd eigenVectors;
	MatrixXcd eigenValues;
};

void readImageHeader(char[], int&, int&, int&, bool&);
void readImage(char[], ImageType&);
void writeImage(char[], ImageType&);
void PrintMatrix(MatrixXd m, int dimension1, int dimension2);
vector<Image> PopulateImages(string directory, int imageWidth, int imageHeight);
vector<Image> obtainTrainingFaces(string directory, int imageWidth, int imageHeight);
VectorXi compAvgFaceVec(const vector<Image> &imageVector);

int main()
{
	string inputString, trainingDataset = "Faces_FA_FB/fa_H", eigenvaluesFile = "eigenValues.txt";
	vector<Image> ImageVector;
	VectorXi avgFaceVector;
	int K=0;
	MatrixXf A(IMG_VEC_LEN, NUM_SAMPLES);
	MatrixXf C(IMG_VEC_LEN, IMG_VEC_LEN);
	MatrixXcf eigenVectors, eigenValues;

	do
	{
		cout << endl
		     << "+=======================================================+\n"
			 << "|Select  0 to obtain training faces (I_1...I_M)         |\n"
			 << "|Select  1 to compute average face vector (Psi)         |\n"
			 << "|Select  2 to compute matrix A ([Phi_i...Phi_M])        |\n"
			 << "|Select  3 to compute covariance matrix C (1/M * AA^T)  |\n"
			 << "|Select  4 to compute the eigenvectors/values of AA^T   |\n"
			 << "|Select  5 to project eigenvalues                       |\n"
		     << "|Select -1 to exit                                      |\n"
		     << "+=======================================================+\n"
		     << endl
		     << "Choice: ";

		cin >> inputString;
		if (inputString == "0")
		{
			ImageVector = obtainTrainingFaces(trainingDataset, IMG_W, IMG_H);
		}
		else if (inputString == "1")
		{
			avgFaceVector = compAvgFaceVec(ImageVector);
			ImageType avgFaceImg(IMG_H, IMG_W, 255);

			for (int i = 0; i < IMG_H; i++)
			{
				for (int j = 0; j < IMG_W; j++)
				{
					int val = avgFaceVector(i * IMG_W + j);
					avgFaceImg.setPixelVal(i, j, val);
					writeImage((char*) "myAvgFaceImg.pgm", avgFaceImg);
				}
			}
		}
		else if (inputString == "2")
		{
			VectorXi phi;
			
			for (int j = 0; j < NUM_SAMPLES; j++)
			{
				phi = ImageVector[j].getFaceVector() - avgFaceVector;

				for (int i = 0; i < IMG_VEC_LEN; i++)
				{
					A(i, j) = phi(i);
				}
			}
		}
		else if (inputString == "3")
		{
			C = A * A.transpose();
			C = C * (1.0f / (float)NUM_SAMPLES);

			// cout << C(0, 0) << " " << C(0, 1) << " " << C(0, 2) << "\n"
			// 	 << C(1, 0) << " " << C(1, 1) << " " << C(1, 2) << "\n"
			// 	 << C(2, 0) << " " << C(2, 1) << " " << C(2, 2) << "\n";
		}
		else if (inputString == "4")
		{
			MatrixXf AT_A(NUM_SAMPLES, NUM_SAMPLES);
			AT_A = A.transpose() * A;

			EigenSolver<MatrixXf> es(AT_A);
			cout << "Computing eigenvalues..." << endl;
			eigenValues = es.eigenvalues();
			cout << "Finished computing eigenvalues!" << endl;

			cout << "Computing eigenvectors..." << endl;
			//eigenVectors = es.eigenvectors();
			cout << "Finished computing eigenvectors!" << endl;

			cout << eigenValues << endl;

			ofstream fout;
			fout.open(eigenvaluesFile.c_str());
			for (int i = 0; i < eigenValues.size(); ++i)
			{
				fout << eigenValues(i,0) << endl;
			}
			fout.close();
		}
		else if (inputString == "5")
		{
			complex<double> threshold;
			complex<double> currentEigenValueNum=0, totalEigenValueNum=0;
			cout << "Select threshold value(0 to 1): ";
			cin >> threshold;

			for (int i = 0; i < eigenValues.size(); ++i)
			{
				totalEigenValueNum += eigenValues(i,0);
			}

			for (int i = 0; i < eigenValues.size(); ++i)
			{
				currentEigenValueNum += eigenValues(i,0);
				if((currentEigenValueNum/totalEigenValueNum).real() >= threshold.real())
				{
					cout << "Found K threshold to save " << threshold.real() << " of info @ K = " << i << endl;
					K = i;
					break;
				}
			}

		}

		cout << endl;
	} while (inputString != "-1");
}

void PrintMatrix(MatrixXd m, int dimension1, int dimension2)
{
	for (int i = 0; i < dimension1; ++i)
	{
		for (int j = 0; j < dimension2; ++j)
		{
			cout << m(i, j) << " ";
		}
		cout << endl;
	}
}

vector<Image> PopulateImages(string directory, int imageWidth, int imageHeight)
{
	vector<Image> returnVector;
	MatrixXd avgFaceMatrix(imageHeight, imageWidth);
	int numFiles = 0;
	int k, j, M, N, Q;
 	bool type;
	int val;
	DIR *dir;
	struct dirent *ent;

	cout << "Populating images..." << endl;

	if ((dir = opendir (directory.c_str())) != NULL) 
	{
		while ((ent = readdir (dir)) != NULL) 
		{
			string temp, fileName = ent->d_name;
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
			for (; i < (ID_LENGTH + DATE_LENGTH + 1); ++i)
			{
				temp += fileName[i];
			}
			imageDate = atoi(temp.c_str());
			i++;

			temp = "";
			for (; i < (ID_LENGTH + DATE_LENGTH + SET_LENGTH + 2); ++i)
			{
				temp += fileName[i];
			}
			imageDataset = temp;

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
			
			string currentFile = directory + '/' + ent->d_name;

			if (fileName != "." && fileName != "..")
			{					
				readImageHeader((char*) currentFile.c_str(), N, M, Q, type);

				// allocate memory for the image array
				ImageType image(N, M, Q);
 				readImage((char*) currentFile.c_str(), image);

 				MatrixXd imageMatrix(N, M);
				
				for (k = 0; k < N; k++)
				{
					for (j = 0; j < M; j++)
					{
			            image.getPixelVal(k, j, val);
			            imageMatrix(k, j) = val;
					}
				}

 				avgFaceMatrix += imageMatrix;
				newImage.ExtractEigenVectors(imageMatrix, N, M);
				returnVector.push_back(newImage);
				//prints pixel values to output file. I am using this to validate the image has been correctly extracted
			}

			numFiles++;
		}
		
		numFiles -= 2; // subtract 2 to correct count due to '.' and '..' being counted

		string avgFaceText = "avgFace.txt", avgFaceImage = "avgFaceImage.pgm";

		int temp;
		ofstream fout;
		ImageType image(imageHeight, imageWidth, Q);
		
		fout.open(avgFaceText.c_str(), std::ofstream::trunc);

		for (int i = 0; i < imageHeight; ++i)
		{
			for (int j = 0; j < imageWidth; ++j)
			{
				avgFaceMatrix(i, j) /= numFiles;
				temp = (int) avgFaceMatrix(i, j);
				image.setPixelVal(i, j, temp);
				fout << temp << " ";
			}

			fout << endl;
		}

		fout.close();
		closedir (dir);

		writeImage((char*)avgFaceImage.c_str(), image);

		cout << "Finished populating images." << endl;

		return returnVector;
	} 
	else 
	{
	 	/* could not open directory */
		return returnVector;
	}
}

vector<Image> obtainTrainingFaces(string directory, int imageWidth, int imageHeight)
{
	vector<Image> returnVector;
	MatrixXd avgFaceMatrix(imageHeight, imageWidth);
	int k, j, M, N, Q;
 	bool type;
	int val;
	DIR *dir;
	struct dirent *ent;

	cout << "Obtaining training faces..." << endl;

	if ((dir = opendir(directory.c_str())) != NULL) 
	{
		while ((ent = readdir(dir)) != NULL) 
		{
			string temp, fileName = ent->d_name;
			int imageID, imageDate, i;
			string imageDataset;
			bool imageGlasses;
			
			for (i = 0; i < ID_LENGTH; ++i)
			{
			    temp += fileName[i];
			}
			imageID = atoi(temp.c_str());
			i++;

			temp = "";
			for (; i < (ID_LENGTH + DATE_LENGTH + 1); ++i)
			{
				temp += fileName[i];
			}
			imageDate = atoi(temp.c_str());
			i++;

			temp = "";
			for (; i < (ID_LENGTH + DATE_LENGTH + SET_LENGTH + 2); ++i)
			{
				temp += fileName[i];
			}
			imageDataset = temp;

			if(fileName[i] == '.')
			{
				imageGlasses = false;
			}
			else
			{
				imageGlasses = true;
			}

			Image currentImage;
			currentImage.SetFileName(ent->d_name);
			currentImage.SetIdentifier(imageID);
			currentImage.SetDateTaken(imageDate);
			currentImage.SetWearingGlasses(imageGlasses);		
			
			string currentFile = directory + '/' + ent->d_name;

			if (fileName != "." && fileName != "..")
			{
				readImageHeader((char*) currentFile.c_str(), N, M, Q, type);

				// allocate memory for the image array
				ImageType tempImage(N, M, Q);

 				readImage((char*) currentFile.c_str(), tempImage);

 				VectorXi imageVector(N * M);
				
				for (k = 0; k < N; k++)
				{
					for (j = 0; j < M; j++)
					{
			            tempImage.getPixelVal(k, j, val);
			            imageVector(k * M + j) = val;
					}
				}

				currentImage.setFaceVector(imageVector);

				returnVector.push_back(currentImage);
			}
		}

		closedir(dir);

		cout << "Finished obtaining training faces." << endl;
	} 
	else
	{
	 	cout << "Error: Could not open directory " << directory << endl;
	}

	return returnVector;
}

VectorXi compAvgFaceVec(const vector<Image> &imageVector)
{
	VectorXi result = VectorXi::Zero(IMG_VEC_LEN);

	for (int i = 0; i < NUM_SAMPLES; i++)
	{
		result += imageVector[i].getFaceVector();
	}

	result /= NUM_SAMPLES;

	return result;
}

void Image::ExtractEigenVectors(MatrixXd m, int dimension1, int dimension2)
{
	MatrixXd m_tm = m.transpose() * m;

	//cout << "m.transpose*m:" << endl;
	//PrintMatrix(m_tm, m_tm.rows(), m_tm.cols());

	EigenSolver<MatrixXd> EigenSolver;
	EigenSolver.compute(m_tm, true); //Initializes eigensolver with something to de-compose
	eigenVectors = EigenSolver.eigenvectors();
	// cout << endl << "EigenVector Matrix: " << endl << eigenVectors << endl;

	eigenValues = EigenSolver.eigenvalues();
	// cout << endl << "EigenValues Matrix: " << endl << eigenValues << endl;
}
