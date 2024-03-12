#include "SSS_data_structure.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#define scalProd(a, b) (a.x * b.x + a.y * b.y + a.z * b.z)
using namespace std;

// general geometry function
vec mkVec(double x, double y, double z)
{
	vec tmp;
	tmp.x = x;
	tmp.y = y;
	tmp.z = z;
	return tmp;
}
vec sum(vec a, vec b)
{
	vec res;
	res.x = a.x + b.x;
	res.y = a.y + b.y;
	res.z = a.z + b.z;
	return res;
}
vec sumK(vec a, vec b, double k)
{
	vec res;
	res.x = a.x + k * b.x;
	res.y = a.y + k * b.y;
	res.z = a.z + k * b.z;
	return res;
}
vec vecToScal(vec a, double k)
{
	vec res;
	res.x = k * a.x;
	res.y = k * a.y;
	res.z = k * a.z;
	return res;
}
vec vecProd(vec a, vec b)
{
	vec product;
	product.x = a.y * b.z - a.z * b.y;
	product.y = a.z * b.x - a.x * b.z;
	product.z = a.x * b.y - a.y * b.x;
	return product;
}
matrix3X3 mk3X3(double a00, double a01, double a02,
				double a10, double a11, double a12,
				double a20, double a21, double a22)
{
	matrix3X3 tmp;
	tmp.element[0][0] = a00;
	tmp.element[0][1] = a01;
	tmp.element[0][2] = a02;
	tmp.element[1][0] = a10;
	tmp.element[1][1] = a11;
	tmp.element[1][2] = a12;
	tmp.element[2][0] = a20;
	tmp.element[2][1] = a21;
	tmp.element[2][2] = a22;
	return tmp;
}
vec MatrixXVec(matrix3X3 M, vec a)
{
	vec tmp;
	tmp.x = M.element[0][0] * a.x + M.element[0][1] * a.y + M.element[0][2] * a.z;
	tmp.y = M.element[1][0] * a.x + M.element[1][1] * a.y + M.element[1][2] * a.z;
	tmp.z = M.element[2][0] * a.x + M.element[2][1] * a.y + M.element[2][2] * a.z;
	return tmp;
}
volume Translate(vec a, volume v0) // translate original volume v0 to new position,
{								   // translation vector is a
	volume tr = v0;
	for (int i = 0; i < 6; i++)
	{
		tr.fc[i].p -= scalProd(a, v0.fc[i].n);
	}
	tr.m0n0 = sum(a, v0.m0n0);
	tr.center = sum(a, v0.center);

	return tr;
}
volume Rotate(vec n, double ang, volume v0)
{
	matrix3X3 Rotation;
	// because rotation around zaxis as usual, better start from n.z=1
	if (n.x == 0 && n.y == 0 && n.z == 1) // rotation around Z axis
	{
		Rotation = mk3X3(cos(ang), -sin(ang), 0,
						 sin(ang), cos(ang), 0,
						 0, 0, 1);
	}
	else
	{
		if (n.x == 0 && n.y == 1 && n.z == 0) // rotation around Y axis
		{
			Rotation = mk3X3(cos(ang), 0, sin(ang),
							 0, 1, 0,
							 -sin(ang), 0, cos(ang));
		}
		else
		{
			if (n.x == 1 && n.y == 0 && n.z == 0) // rotation around X axis
			{
				Rotation = mk3X3(1, 0, 0,
								 0, cos(ang), -sin(ang),
								 0, sin(ang), cos(ang));
			}
			else
			{
				cout << "Wrong rotation axis! (" << n.x << "," << n.y << "," << n.z << ")" << endl;
				cout << "Volume can be rotated only around x,y or z axises only!!!" << endl;
				return v0;
			}
		}
	}
	// rotation around Z axis counter clockwise on angle ang
	volume Vr = v0;
	for (int i = 0; i < 6; i++)
	{
		Vr.fc[i].n = MatrixXVec(Rotation, v0.fc[i].n);
	}
	Vr.m0n0 = MatrixXVec(Rotation, v0.m0n0);
	Vr.center = MatrixXVec(Rotation, v0.center);

	Vr.translateX = MatrixXVec(Rotation, v0.translateX);
	Vr.translateZ = MatrixXVec(Rotation, v0.translateZ);
	return Vr;
}
matrix3X3 RotationMatrix(vec n, double ang)
{
	// simulate the matrix for rotation and then translation
	matrix3X3 Rotation;
	if (n.x == 1 && n.y == 0 && n.z == 0) // rotation around X axis
	{
		Rotation = mk3X3(1, 0, 0,
						 0, cos(ang), -sin(ang),
						 0, sin(ang), cos(ang));
	}
	else
	{
		if (n.x == 0 && n.y == 1 && n.z == 0) // rotation around Y axis
		{
			Rotation = mk3X3(cos(ang), 0, sin(ang),
							 0, 1, 0,
							 -sin(ang), 0, cos(ang));
		}
		else
		{
			if (n.x == 0 && n.y == 0 && n.z == 1) // rotation around Z axis
			{
				Rotation = mk3X3(cos(ang), -sin(ang), 0,
								 sin(ang), cos(ang), 0,
								 0, 0, 1);
			}
			else
			{
				cout << "!!!!!!!!Wrong rotation axis! (" << n.x << "," << n.y << "," << n.z << ")" << endl;
				cout << "Volume can be rotated only around x,y or z axises only!!!" << endl;
			}
		}
	}
	return Rotation;
}