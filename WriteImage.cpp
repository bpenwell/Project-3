#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

#include "image.h"
#include "rgb.h"

void writeImage(char fname[], ImageType& image)
{
    int i, j;
    int N, M, Q;
    unsigned char *charImage;
    ofstream ofp;

    image.getImageInfo(N, M, Q);

    charImage = (unsigned char *) new unsigned char [3*M*N];

    int val;

    for(i=0; i<N; i++)
    {
        for(j=0; j<M; j++)
        {
            image.getPixelVal(i, j/3, val);
            charImage[i+j]=(unsigned char)val;
            //charImage[i*3*M+j+1]=(unsigned char)val.g;
            //charImage[i*3*M+j+2]=(unsigned char)val.b;
        }
    }

    ofp.open(fname, ios::out | ios::binary);

    if (!ofp)
    {
        cout << "Can't open file: " << fname << endl;
        exit(1);
    }

    ofp << "P5" << " " << M << " " << N << " " << 255;

    ofp.write( reinterpret_cast<char *>(charImage), (3*M*N)*sizeof(unsigned char));

    if (ofp.fail())
    {
        cout << "Can't write image " << fname << endl;
        exit(0);
    }

    ofp.close();

    delete [] charImage;
}
