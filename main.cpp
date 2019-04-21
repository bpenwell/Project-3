#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <sys/stat.h>
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
vector<Image> obtainTrainingFaces(string directory, int imageWidth, int imageHeight);
VectorXi compAvgFaceVec(const vector<Image> &imageVector);
void getMinMax(VectorXd,double&,double&);
int main()
{
	string inputString;
	string trainingDataset = "Faces_FA_FB/fa_H",
		   testingDataset = "Faces_FA_FB/fb_H",
		   eigenvaluesFile = "eigenValues.txt",
		   eigenvectorsFile = "eigenVectors.txt",
		   imageCoefficientsFile = "imageCoefficients.txt";
	vector<Image> ImageVector;
	vector<Image> TestImageVector;
	VectorXi avgFaceVector;
	vector<VectorXi> phi;
	unsigned K = 0;
	MatrixXd A(IMG_VEC_LEN, NUM_SAMPLES);
	MatrixXd C(IMG_VEC_LEN, IMG_VEC_LEN);
	MatrixXd eigenVectors_u, eigenVectors_v;
	VectorXd eigenValues;
	vector<EigenValVecPair> EigenValVecPairs;

	do
	{
		cout << endl
		     << "+==============================================================+\n"
			 << "|Select  0 to obtain training images (I_1...I_M)               |\n"
			 << "|Select  1 to compute average face vector (Psi)                |\n"
			 << "|Select  2 to compute matrix A ([Phi_i...Phi_M])               |\n"
			 << "|Select  3 to compute the eigenvectors/values of A^TA          |\n"
			 << "|Select  4 to project eigenvalues (req: 0,1,2)                 |\n"
			 << "|Select  5 to visualize the 10 largest & smallesteigenvectors  |\n"
		     << "|Select  6 to recognize a face                                 |\n"
		     << "|Select  7 to obtain testing images                            |\n"
		     << "|Select -1 to exit                                             |\n"
		     << "+==============================================================+\n"
		     << endl
		     << "Choice: ";

		cin >> inputString;
		if (inputString == "0") //Initialize images
		{
			ImageVector = obtainTrainingFaces(trainingDataset, IMG_W, IMG_H);
			cout << "ImageVector.size() = " << ImageVector.size() << endl;
		}
		else if (inputString == "1") //Generate Psi
		{
			avgFaceVector = compAvgFaceVec(ImageVector);
			ImageType avgFaceImg(IMG_H, IMG_W, 255);

			for (int i = 0; i < IMG_H; i++)
			{
				for (int j = 0; j < IMG_W; j++)
				{
					int val = avgFaceVector(i * IMG_W + j);
					avgFaceImg.setPixelVal(i, j, val);
				}
			}

			writeImage((char*) "myAvgFaceImg.pgm", avgFaceImg);
		}
		else if (inputString == "2") // Generate Phi
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
		else if (inputString == "3") //Generate eigenvalues/vectors
		{
			MatrixXd AT_A(NUM_SAMPLES, NUM_SAMPLES);
			AT_A = A.transpose() * A;

			EigenSolver<MatrixXd> es(AT_A);
			cout << "Computing eigenvalues..." << endl;
			eigenValues = es.eigenvalues().real();
			cout << "Finished computing eigenvalues!" << endl;

			cout << "Computing eigenvectors..." << endl;
			eigenVectors_v = es.eigenvectors().real();
			eigenVectors_u = A * eigenVectors_v;
			eigenVectors_u.colwise().normalize();
			cout << "Finished computing eigenvectors!" << endl;

			for (int i = 0; i < eigenValues.rows(); i++)
			{
				EigenValVecPair pair;
				pair.eigenValue = eigenValues(i);
				pair.eigenVector = eigenVectors_u.col(i);
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
		}
		else if (inputString == "4") // Generate Omega based off threshold
		{
			double threshold, currentEigenValueNum = 0, totalEigenValueNum = 0;
			cout << "Select threshold value (0 to 1): ";
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
				for (unsigned i = 0; i < IMG_VEC_LEN; i++)
				{
					fin_vecs >> eigenVector(i);
				}
				
				topEigenVectors.push_back(eigenVector);
			}
			fin_vecs.close();

			for (unsigned i = 0; i < topEigenValues.size(); ++i)
			{
				totalEigenValueNum += topEigenValues[i];
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

			string omegaVectorsFile = "omegaVectors.txt";
			ofstream fout_omega_vecs(omegaVectorsFile.c_str());

			VectorXd phi_hat = VectorXd::Zero(IMG_VEC_LEN);

			for (unsigned i = 0; i < phi.size(); i++)
			{
				for (unsigned j = 0; j < K; j++)
				{
					MatrixXd u_t = ((MatrixXd)topEigenVectors[j]).transpose();
					MatrixXd phi_i = (phi[i]).cast<double>();
					double w = (u_t * phi_i)(0, 0);
					fout_omega_vecs << w;

					if (j < K - 1)	// if not last element
						fout_omega_vecs << " ";

					if (i == 0)
					{
						phi_hat += w * topEigenVectors[j];
					}
				}

				if (i < phi.size() - 1)	// if not last element
					fout_omega_vecs << endl;
			}

			phi_hat += avgFaceVector.cast<double>();
			VectorXi phi_hat_int = phi_hat.cast<int>();

			ImageType reconstructedFaceImg(IMG_H, IMG_W, 255);

			for (int i = 0; i < IMG_H; i++)
			{
				for (int j = 0; j < IMG_W; j++)
				{
					int val = phi_hat_int(i * IMG_W + j);
					reconstructedFaceImg.setPixelVal(i, j, val);
				}
			}
			
			writeImage((char*) "reconstructedFaceImg.pgm", reconstructedFaceImg);

			fout_omega_vecs.close();
		}
		else if (inputString == "5") // Generate 10 largest/smallest eigenvectors
		{
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
				for (unsigned i = 0; i < IMG_VEC_LEN; i++)
				{
					fin_vecs >> eigenVector(i);
				}
				
				topEigenVectors.push_back(eigenVector);
			}
			fin_vecs.close();


			for (int i = 0; i < 10; ++i)
			{
				double min,max;
				getMinMax(topEigenVectors[i],min,max);
				cout << "min: " << min << endl;
				cout << "max: " << max << endl;
				ImageType eigenFace(IMG_H, IMG_W, 255);
				for (int j = 0; j < IMG_H; j++)
				{
					for (int k = 0; k < IMG_W; k++)
					{
						
						double val = (topEigenVectors[i](j * IMG_W + k) - min) * (255 / (max - min));
						val += avgFaceVector(j * IMG_W + k);
						cout << val << endl;
						eigenFace.setPixelVal(j, k, (int)val);
					}
				}
				
				string folder = "largestEigenFaces/";
				cout << "folder->" << folder << endl;
				mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
				
				stringstream ss;
				ss << i;
				string str = ss.str();
				folder += str;
				folder += ".pgm";
				
				ofstream fout;
				fout.open(folder.c_str());
				fout.close();

				writeImage((char*) folder.c_str(), eigenFace);
			}
			
			vector<VectorXd> temp = topEigenVectors;
			temp.erase(temp.begin(), temp.end()-10);
			for (int i = 0; i < 10; ++i)
			{
				double min,max;
				getMinMax(temp[i],min,max);
				cout << "min: " << min << endl;
				cout << "max: " << max << endl;
				ImageType eigenFace(IMG_H, IMG_W, 255);
				for (int j = 0; j < IMG_H; j++)
				{
					for (int k = 0; k < IMG_W; k++)
					{
						
						double val = (temp[i](j * IMG_W + k) - min) * (255 / (max - min));
						val += avgFaceVector(j * IMG_W + k);
						cout << val << endl;
						eigenFace.setPixelVal(j, k, (int)val);
					}
				}
				
				string folder = "smallestEigenFaces/";
				cout << "folder->" << folder << endl;
				mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
				
				stringstream ss;
				ss << i;
				string str = ss.str();
				folder += str;
				folder += ".pgm";
				
				ofstream fout;
				fout.open(folder.c_str());
				fout.close();

				writeImage((char*) folder.c_str(), eigenFace);
			}
		}
		else if (inputString == "6")
		{
			string fileName;
			cout << "Input image filename you would like to check for in the database (format- nnnnn_yymmdd_xx.pgm OR nnnnn_yymmdd_xx_q.pgm):" << endl;
			cin >> fileName;

			string temp;
			int imageID, imageDate, i, N, M, Q;
			string imageDataset;
			bool imageGlasses, type;
			
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
			currentImage.SetFileName(fileName);
			currentImage.SetIdentifier(imageID);
			currentImage.SetDateTaken(imageDate);
			currentImage.SetWearingGlasses(imageGlasses);		
			
			//string currentFile = directory + '/' + ent->d_name;
			int val;

			if (fileName != "." && fileName != "..")
			{
				cout << fileName << endl;
				readImageHeader((char*) fileName.c_str(), N, M, Q, type);

				// allocate memory for the image array
				ImageType tempImage(N, M, Q);

 				readImage((char*) fileName.c_str(), tempImage);

 				VectorXi imageVector(N * M);
				
				for (int k = 0; k < N; k++)
				{
					for (int j = 0; j < M; j++)
					{
			            tempImage.getPixelVal(k, j, val);
			            imageVector(k * M + j) = val;
					}
				}
				//cout << "before...";
				currentImage.setFaceVector(imageVector);
				//cout << "imageVector->rows: "<< imageVector.rows() << " ImageVector->cols: "<< imageVector.cols();
				//cout << "after!" << endl;

				//returnVector.push_back(currentImage);
			}
			//cout << "avgFaceVector->rows: "<< avgFaceVector.rows() << " avgFaceVector->cols: "<< avgFaceVector.cols();
			MatrixXi imagePhi = currentImage.getFaceVector() - avgFaceVector;
			MatrixXd imagePhi_double = imagePhi.cast<double>();

			//cout << "about to extract" << endl;
			currentImage.ExtractEigenVectors(imagePhi_double, IMG_W, IMG_H);


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
				for (unsigned i = 0; i < IMG_VEC_LEN; i++)
				{
					fin_vecs >> eigenVector(i);
				}
				
				topEigenVectors.push_back(eigenVector);
			}
			fin_vecs.close();

			string omegaVectorsFile = "omegaVectors.txt";
			vector<VectorXd> omegaVectors;

			ifstream fin_omegas;
			fin_omegas.open(omegaVectorsFile.c_str());
			
			VectorXd readVector(K);

			for (int i = 0; i < NUM_SAMPLES; ++i)
			{

				for (unsigned j = 0; j < K; ++j)
				{
					fin_omegas >> readVector(j);
				}
				omegaVectors.push_back(readVector);
			}

			fin_omegas.close();

			VectorXd wVector(K);

			VectorXd phi_hat = VectorXd::Zero(IMG_VEC_LEN);
			for (unsigned j = 0; j < K; j++)
			{
				MatrixXd u_t = ((MatrixXd)topEigenVectors[j]).transpose();
				MatrixXd phi_i = (imagePhi).cast<double>();
				double w = (u_t * phi_i)(0, 0);
				wVector(j) = w;
		 		//cout << w << " ";
		 		//fout_omega_vecs << w;

				if (j < K - 1)	// if not last element
					//fout_omega_vecs << " ";

				phi_hat += w * topEigenVectors[j];
			}
			cout << endl;

			phi_hat += avgFaceVector.cast<double>();
			MatrixXi phi_hat_int = phi_hat.cast<int>();

			//cout << wVector << endl;

			cout << "calculating e_r..." << endl;
			double e_r;
			//Compare omegaVectors for norm. min of all those is e_r
			for (int i = 0; i < NUM_SAMPLES; ++i)
			{
				double diffSum=0;
				for (unsigned j = 0; j < K; ++j)
				{
					//(image w - database image w)^2
					diffSum += (wVector(j) - omegaVectors[i](j)) * (wVector(j) - omegaVectors[i](j));
				}

				if(i==0)
				{
					e_r = diffSum;
				}
				else if(e_r > diffSum)
				{
					//cout << "New min-> " << e_r << endl;
					e_r = diffSum;
				}
			}
			cout << "e_r: " << e_r << endl;
		}
		else if (inputString == "7") //Initialize images
		{
			TestImageVector = obtainTrainingFaces(testingDataset, IMG_W, IMG_H);
			cout << "ImageVector.size() = " << ImageVector.size() << endl;
		}

		cout << endl;
	} while (inputString != "-1");
}

void getMinMax(VectorXd m, double& min, double& max)
{
	min = 1, max = -1;
	for (int j = 0; j < IMG_H; j++)
	{
		for (int k = 0; k < IMG_W; k++)
		{
			if(m(j * IMG_W + k) < min)
			{
				min = m(j * IMG_W + k);
			}
			if(m(j * IMG_W + k) > max)
			{
				max = m(j * IMG_W + k);
			}
		}
	}
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