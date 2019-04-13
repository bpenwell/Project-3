#ifndef RGB_H
#define RGB_H

struct RGB 
{
	RGB();
    RGB(int, int, int);
    RGB& operator=(RGB);
    bool operator==(RGB);
    int r, g, b;
};

#endif
