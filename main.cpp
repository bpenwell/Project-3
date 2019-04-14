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
const int DIM = 3;
const int IMG_W = 48, IMG_H = 60;
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
	void setFaceVector(VectorXd vector) { faceVector = vector; };

	// Getters
	int GetIdentifier() { return identifier; };
	int GetDateTaken() { return dateTaken; };
	string GetFileName() { return fileName; };
	string GetDataset() { return dataset; };
	bool GetWearingGlasses() { return wearingGlasses; };
	VectorXd getFaceVector() { return faceVector; };

private:
	string fileName;
	int identifier; //nnnnn
	int dateTaken; //yymmdd
	string dataset;
	bool wearingGlasses;

	VectorXd faceVector;
	MatrixXcd eigenVectors;
	MatrixXcd eigenValues;
};

void readImageHeader(char[], int&, int&, int&, bool&);
void readImage(char[], ImageType&);
void writeImage(char[], ImageType&);
void PrintMatrix(MatrixXd m, int dimension1, int dimension2);
vector<Image> PopulateImages(string directory, int imageWidth, int imageHeight);
vector<Image> obtainTrainingFaces(string directory, int imageWidth, int imageHeight);

int main()
{
	// MatrixXd g(DIM, DIM);
	string trainingDataset = "Faces_FA_FB/fa_H";

	// g << 1, 0, 0, 2, 1, 1, 0, 0, 1;
	// PrintMatrix(g, DIM, DIM);

	vector<Image> ImageVector;
	//testImage.ExtractEigenVectors(g, DIM, DIM);
	string inputString;
	do
	{
		cout << endl
		     << "+=======================================================+\n"
			 << "|Select  0 to obtain face images                        |\n"
			 << "|Select  1 to generate average face                     |\n"
		     << "|Select -1 to exit                                      |\n"
		     << "+=======================================================+\n"
		     << endl
		     << "Choice: ";

		cin >> inputString;
		if(inputString == "0")
		{
			ImageVector = obtainTrainingFaces(trainingDataset, IMG_W, IMG_H);
		}
		else if(inputString == "1")
		{
			ImageVector = PopulateImages(trainingDataset, 48, 60);
		}

		cout << endl;
	} while(inputString != "-1");


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

 				VectorXd imageVector(N * M);
				
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
