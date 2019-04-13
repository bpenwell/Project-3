#include <iostream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

#include "image.h"
#include "rgb.h"

ImageType::ImageType()
{
    N = 0;
    M = 0;
    Q = 0;

    pixelValue = NULL;
}

ImageType::ImageType(int tmpN, int tmpM, int tmpQ)
{
    int i, j;

    N = tmpN;
    M = tmpM;
    Q = tmpQ;
    
    pixelValue = new int* [N];
    RGBpixelValue = new RGB* [N];
    for(i=0; i<N; i++)
    {
        RGBpixelValue[i] = new RGB[M];
        pixelValue[i] = new int[M];
        for(j=0; j<M; j++)
            pixelValue[i][j] = 0;
    }
}

void ImageType::getImageInfo(int& rows, int& cols, int& levels)
{
    rows = N;
    cols = M;
    levels = Q;
} 

void ImageType::setImageInfo(int rows, int cols, int levels)
{
    N = rows;
    M = cols;
    Q = levels;
} 

void ImageType::setPixelVal(int i, int j, int val)
{
    pixelValue[i][j] = val;
}

void ImageType::getPixelVal(int i, int j, int& val)
{
    val = pixelValue[i][j];
}

void ImageType::setPixelVal(int i, int j, RGB val)
{
    //GOES WRONG HERE
    RGBpixelValue[i][j] = val;
}

void ImageType::getPixelVal(int i, int j, RGB& val)
{
    val = RGBpixelValue[i][j];
}