#include "rgb.h"

RGB::RGB() : r(0), g(0), b(0)
{}

RGB::RGB(int r, int g, int b)
{
    this->r = r;
    this->g = g;
    this->b = b;
}

RGB& RGB::operator=(RGB rhs)
{
    this->r = rhs.r;
    this->g = rhs.g;
    this->b = rhs.b;

    return *this;
}

bool RGB::operator==(RGB rhs)
{
	if (this->r == rhs.r && this->g == rhs.g && this->b == rhs.b)
    {
		return true;
	}
	else
    {
		return false;
	}
}