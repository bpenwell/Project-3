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
using Eigen::VectorXf;
using Eigen::VectorXcf;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
<<<<<<< HEAD
=======
using Eigen::MatrixXf;
>>>>>>> master
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
//const int DIM = 3;
const int IMG_W = 48, IMG_H = 60, IMG_VEC_LEN = IMG_H * IMG_W;
const int NUM_SAMPLES = 1204;
const string trainingFaces = "trainingFaces.txt", avgFaceText = "avgFace.txt", avgFaceImage = "avgFaceImage.pgm", eigenCoefficientsFile = "CoefficientsFile.txt";

class Image
{
public:
	void ExtractEigenVectors(MatrixXd m, int dimension1, int dimension2);
	//void ExtractEigenVectors(MatrixXi m, int dimension1, int dimension2);

	// Setters
	void SetIdentifier(int input) { identifier = input; };
	void SetDateTaken(int input) { dateTaken = input; };
	void SetDataset(string input) { dataset = input; };
	void SetWearingGlasses(bool input) { wearingGlasses = input; };
	void SetFileName(string input) { fileName = input; };
	void setFaceVector(VectorXi vector) { faceVector = vector; };

	// Getters
	int GetIdentifier() const { return identifier; };
	int GetDateTaken() const { return dateTaken; };
	string GetFileName() const { return fileName; };
	string GetDataset() const { return dataset; };
	bool GetWearingGlasses() const { return wearingGlasses; };
	VectorXi getFaceVector() const { return faceVector; };
	MatrixXcd GetEigenValues() const { return eigenValues; };

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
int comp_K_Threshold(MatrixXcd eigenvalues, double targetPercentInfoRetained);
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
<<<<<<< HEAD
	// MatrixXd g(DIM, DIM);
	string trainingDataset = "Faces_FA_FB/fa_H";

	vector<Image> ImageVector;
	vector<Image> projectedImageVector;
	string inputString;
=======
	int K=0;
	MatrixXf A(IMG_VEC_LEN, NUM_SAMPLES);
	MatrixXf C(IMG_VEC_LEN, IMG_VEC_LEN);
	MatrixXf eigenVectors;
	VectorXf eigenValues;

>>>>>>> master
	do
	{
		cout << endl
		     << "+=======================================================+\n"
<<<<<<< HEAD
			 << "|Select  0 to obtain training faces                     |\n"
			 << "|Select  1 to compute average face vector               |\n"
			 << "|Select  2 to project the images (PCA)                  |\n"
=======
			 << "|Select  0 to obtain training faces (I_1...I_M)         |\n"
			 << "|Select  1 to compute average face vector (Psi)         |\n"
			 << "|Select  2 to compute matrix A ([Phi_i...Phi_M])        |\n"
			 << "|Select  3 to compute the eigenvectors/values of AA^T   |\n"
			 << "|Select  4 to project eigenvalues                       |\n"
>>>>>>> master
		     << "|Select -1 to exit                                      |\n"
		     << "+=======================================================+\n"
		     << endl
		     << "Choice: ";

		cin >> inputString;
		if (inputString == "0")
		{
			ImageVector = obtainTrainingFaces(trainingDataset, IMG_W, IMG_H);
			
			VectorXi holdFace;
			ofstream fout;
			fout.open(trainingFaces.c_str());
			for (uint i = 0; i < ImageVector.size(); ++i)
			{
				holdFace = ImageVector[i].getFaceVector();
				for (int j = 0; j < holdFace.size(); ++j)
				{
					fout << holdFace[j] << " ";
				}
				fout << endl;
			}
		}
<<<<<<< HEAD
		else if(inputString == "1")
=======
		else if (inputString == "1")
>>>>>>> master
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
<<<<<<< HEAD
		else if(inputString == "2")
		{
			cout << "You must have run selection 0 & 1 in order to run this." << endl;
			cout << "Select threshold value:" << endl;
			double threshold;
			cin >> threshold;
			
			VectorXi differentialImageVector(ImageVector.size());
			//MatrixXi diffMatrix(ImageVector.size(),ImageVector.size());
			vector<VectorXi> differentialImageVectorVector;
			ofstream fout;
			fout.open(eigenCoefficientsFile.c_str());
			for (uint i = 0; i < ImageVector.size(); ++i)
			{
				differentialImageVector = ImageVector[i].getFaceVector() - avgFaceVector;
				differentialImageVectorVector.push_back(differentialImageVector);
				//diffMatrix.col(i) = differentialImageVector;
				fout << differentialImageVector.transpose() << endl;
			}
			fout.close();

			int s = ImageVector.size();
			MatrixXi CovarianceMatrix(s,s);

			//Compute covariance matrix
			//Ensure using A^T*A (MxM matrix)
			//Not A*A^T (N^2xN^2 matrix)

			//NOT RETURNING A SYMMETRIC MATRIX. Needs to be?
			for (uint i = 0; i < s; ++i)
			{
				for (int j = 0; j < s; ++j)
				{
					if(i == 0)
					{
						CovarianceMatrix(i,j) = differentialImageVectorVector[i].row(j)*differentialImageVectorVector[i].transpose().col(j);
					}
					else
					{
						CovarianceMatrix(i,j) += differentialImageVectorVector[i].row(j)*differentialImageVectorVector[i].transpose().col(j);
						cout << "differentialImageVectorVector[i].transpose().col(j): " << differentialImageVectorVector[i].transpose().col(j) << endl;
						cout << "differentialImageVectorVector[i].row(j): " << differentialImageVectorVector[i].row(j) << endl;
					}
				}
			}

			for (int i = 0; i < s; ++i)
			{
				for (int j = 0; j < s; ++j)
				{
					CovarianceMatrix(i,j) /= s;
				}
			}
			//ExtractEigenVectors(CovarianceMatrix, CovarianceMatrix.rows(), CovarianceMatrix.cols());
			cout << CovarianceMatrix << endl;
			cout << "CovarianceMatrix.rows(): " << CovarianceMatrix.rows() << endl;
			cout << "CovarianceMatrix.cols(): " << CovarianceMatrix.cols() << endl;
=======
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
			MatrixXf AT_A(NUM_SAMPLES, NUM_SAMPLES);
			AT_A = A.transpose() * A;


			EigenSolver<MatrixXf> es(AT_A);
			cout << "Computing eigenvalues..." << endl;
			eigenValues = es.eigenvalues().real();
			cout << "Finished computing eigenvalues!" << endl;

			cout << "Computing eigenvectors..." << endl;
			eigenVectors = es.eigenvectors().real();
			eigenVectors = A * eigenVectors;
			eigenVectors.colwise().normalize();
			cout << "Finished computing eigenvectors!" << endl;

			ofstream fout;

			fout.open(eigenvaluesFile.c_str());
			fout << eigenValues;
			fout.close();

			fout.open(eigenvectorsFile.c_str());
			fout << eigenVectors;
			fout.close();

			// cout << sqrt(((eigenVectors.col(0)).dot(eigenVectors.col(0)))) << endl;
		}
		else if (inputString == "4")
		{
			double threshold,currentEigenValueNum=0, totalEigenValueNum=0;
			cout << "Select threshold value(0 to 1): ";
			cin >> threshold;
			vector<double> valuesVector;
			ifstream fin;

			double fileInput;
			fin.open(eigenvaluesFile.c_str());
			while(!fin.eof())
			{
				fin >> fileInput;
				valuesVector.push_back(fileInput);
			}

	  		std::cout.precision(2);
			cout << std::fixed;
			for (int i = 0; i < valuesVector.size(); ++i)
			{
				totalEigenValueNum += valuesVector[i];
				cout << valuesVector[i] << endl;
			}

			for (int i = 0; i < valuesVector.size(); ++i)
			{
				currentEigenValueNum += valuesVector[i];
				if((currentEigenValueNum/totalEigenValueNum) >= threshold)
				{
					cout << "Found K threshold to save " << threshold << " of info @ K = " << i << endl;
					K = i;
					break;
				}
			}

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

>>>>>>> master
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

int comp_K_Threshold(MatrixXcd eigenvalues, double targetPercentInfoRetained)
{
	double fullData=0, partialData=0;
	double temp;
	cout << "eigenvalues (rows,cols): " << eigenvalues.rows() << " " << eigenvalues.cols() << endl;
	for (int i = 0; i < eigenvalues.rows(); ++i)
	{
		for (int j = 0; j < eigenvalues.cols(); ++j)
		{
			fullData += eigenvalues(i,j).real();
			//fullData += temp;
		}
	}

	int k=0;
	for (int i = 0; i < eigenvalues.rows(); ++i)
	{
		for (int j = 0; j < eigenvalues.cols(); ++j)
		{
			k++;
			partialData += eigenvalues(i).real();
			//partialData += temp;
			//cout << (partialData/fullData) << endl;
			if((partialData/fullData) >= targetPercentInfoRetained)
			{
				cout << "Threshold hit at: " << k << endl;
				return k;
			}
		}
	}
	return -1;
}

void Image::ExtractEigenVectors(MatrixXd m, int dimension1, int dimension2)
{
	MatrixXd m_tm = m.transpose() * m;

	//cout << "m.transpose*m:" << endl;
	//PrintMatrix(m_tm, m_tm.rows(), m_tm.cols());

	EigenSolver<MatrixXd> EigenSolver;
	EigenSolver.compute(m_tm, true); //Initializes eigensolver with something to de-compose
	eigenVectors = EigenSolver.eigenvectors();
	//cout << endl << "EigenVector Matrix: " << endl << eigenVectors << endl;

	eigenValues = EigenSolver.eigenvalues();
	//cout << endl << "EigenValues Matrix: " << endl << eigenValues << endl;
}

/*void Image::ExtractEigenVectors(MatrixXi m, int dimension1, int dimension2)
{
	MatrixXi m_tm = m.transpose() * m;

	//cout << "m.transpose*m:" << endl;
	//PrintMatrix(m_tm, m_tm.rows(), m_tm.cols());

	EigenSolver<MatrixXi> EigenSolver;
	EigenSolver.compute(m_tm, true); //Initializes eigensolver with something to de-compose
	eigenVectors = EigenSolver.eigenvectors();
	//cout << endl << "EigenVector Matrix: " << endl << eigenVectors << endl;

	eigenValues = EigenSolver.eigenvalues();
	//cout << endl << "EigenValues Matrix: " << endl << eigenValues << endl;
}*/
