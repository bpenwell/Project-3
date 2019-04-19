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
using Eigen::VectorXcf;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
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
const int IMG_W = 48, IMG_H = 60, IMG_VEC_LEN = IMG_H * IMG_W;
const int NUM_SAMPLES = 1204;

struct EigenValVecPair
{
	double eigenValue;
	VectorXd eigenVector;
};

struct by_eigenValue
{ 
    bool operator()(EigenValVecPair const &a, EigenValVecPair const &b) const 
    { 
        return abs(a.eigenValue) > (b.eigenValue);
    }
};

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
	string inputString;
	string trainingDataset = "Faces_FA_FB/fa_H",
		   eigenvaluesFile = "eigenValues.txt",
		   eigenvectorsFile = "eigenVectors.txt",
		   imageCoefficientsFile = "imageCoefficients.txt";
	vector<Image> ImageVector;
	VectorXi avgFaceVector;
	vector<VectorXi> phi;
	int K=0;
	MatrixXd A(IMG_VEC_LEN, NUM_SAMPLES);
	MatrixXd C(IMG_VEC_LEN, IMG_VEC_LEN);
	MatrixXd eigenVectors;
	VectorXd eigenValues;
	vector<EigenValVecPair> EigenValVecPairs;

	do
	{
		cout << endl
		     << "+=======================================================+\n"
			 << "|Select  0 to obtain training faces (I_1...I_M)         |\n"
			 << "|Select  1 to compute average face vector (Psi)         |\n"
			 << "|Select  2 to compute matrix A ([Phi_i...Phi_M])        |\n"
			 << "|Select  3 to compute the eigenvectors/values of AA^T   |\n"
			 << "|Select  4 to project eigenvalues                       |\n"
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
			for (int j = 0; j < NUM_SAMPLES; j++)
			{
				phi.push_back(ImageVector[j].getFaceVector() - avgFaceVector);

				for (int i = 0; i < IMG_VEC_LEN; i++)
				{
					A(i, j) = phi[j](i);
				}
			}
		}
		else if (inputString == "3")
		{
			MatrixXd AT_A(NUM_SAMPLES, NUM_SAMPLES);
			AT_A = A.transpose() * A;

			EigenSolver<MatrixXd> es(AT_A);
			cout << "Computing eigenvalues..." << endl;
			eigenValues = es.eigenvalues().real();
			cout << "Finished computing eigenvalues!" << endl;

			cout << "Computing eigenvectors..." << endl;
			eigenVectors = es.eigenvectors().real();
			eigenVectors = A * eigenVectors;
			eigenVectors.colwise().normalize();
			cout << "Finished computing eigenvectors!" << endl;

			for (int i = 0; i < eigenValues.rows(); i++)
			{
				EigenValVecPair pair;
				pair.eigenValue = eigenValues(i);
				pair.eigenVector = eigenVectors.col(i);
				EigenValVecPairs.push_back(pair);
			}

			std::sort(EigenValVecPairs.begin(), EigenValVecPairs.end(), by_eigenValue());

			ofstream fout_vals, fout_vecs;

			fout_vals.open(eigenvaluesFile.c_str());
			fout_vecs.open(eigenvectorsFile.c_str());

			for (unsigned i = 0; i < EigenValVecPairs.size(); i++)
			{
				fout_vals << EigenValVecPairs[i].eigenValue;
				fout_vecs << EigenValVecPairs[i].eigenVector.transpose();

				if (i < EigenValVecPairs.size() - 1) // if we're not on last iteration
				{
					fout_vals << endl;
					fout_vecs << endl;
				}
			}

			fout_vals.close();
			fout_vecs.close();

			// cout << sqrt(((eigenVectors.col(0)).dot(eigenVectors.col(0)))) << endl;
		}
		else if (inputString == "4")
		{
			double threshold,currentEigenValueNum=0, totalEigenValueNum=0;
			cout << "Select threshold value(0 to 1): ";
			cin >> threshold;
			vector<double> topEigenValues;
			vector<VectorXd> topEigenVectors;
			ifstream fin_vals, fin_vecs;
			double eigenValue;
			VectorXd eigenVector(IMG_VEC_LEN);

			fin_vals.open(eigenvaluesFile.c_str());
			while(!fin_vals.eof())
			{
				fin_vals >> eigenValue;
				topEigenValues.push_back(eigenValue);
			}
			fin_vals.close();

			fin_vecs.open(eigenvectorsFile.c_str());
			while(!fin_vecs.eof())
			{
				for (unsigned i = 0; i < topEigenValues.size(); i++)
				{
					fin_vecs >> eigenVector(i);
				}
				
				topEigenVectors.push_back(eigenVector);
			}
			fin_vecs.close();
			
			for (unsigned i = 0; i < topEigenValues.size(); ++i)
			{
				totalEigenValueNum += topEigenValues[i];
				cout << topEigenValues[i] << endl;
			}

			for (unsigned i = 0; i < topEigenValues.size(); ++i)
			{
				currentEigenValueNum += topEigenValues[i];
				if((currentEigenValueNum/totalEigenValueNum) >= threshold)
				{
					cout << "Found K threshold to save " << threshold << " of info @ K = " << i << endl;
					K = i;
					break;
				}
			}

			topEigenValues.erase(topEigenValues.begin() + K, topEigenValues.end());
			topEigenVectors.erase(topEigenVectors.begin() + K, topEigenVectors.end());

			/*cout << "Projecting all faces into K dimensions..." << endl;
			VectorXi phi;
			ofstream fout;
			fout.open(imageCoefficientsFile.c_str());
			for (int j = 0; j < NUM_SAMPLES; j++)
			{
				phi = (VectorXi)ImageVector[j].getFaceVector() - avgFaceVector;
				EigenSolver<VectorXi> es(phi);
				cout << "Computing eigenvalues..." << endl;
				VectorXcf imageEigenValues = es.eigenvalues().real();
				cout << "Finished computing eigenvalues!" << endl;
				for (int i = 0; i < K; ++i)
				{
					fout << imageEigenValues(i) << " ";
				}
				fout << endl;
			}*/

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