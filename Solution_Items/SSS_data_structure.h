#ifndef SSS_DATA_STRUCTURE_H
#define SSS_DATA_STRUCTURE_H

#pragma pack(push)
#include <vector>
#include <math.h>
using namespace std;


// general structures
struct vec
{
	double x, y, z; // (x,y,z)
};
struct ray
{		   // ray:
	vec p; // origin in point out(x,y,z)
	vec a; // direction vector a(x,y,z)//normalize!!!
};
struct plain
{			  // plane
	vec n;	  // normal vector n(x,y,z)//normalize
	double p; // distance from (0,0,0)
};
struct matrix3X3
{
	double element[3][3];
};
struct cryst
{
	int V;		  // index of the volume in V[]
	int detector; // index of detector inListOfDetecor[]
	int sector;	  // sectror
	int moduleX;  // index of the module in X
	int moduleZ;  // index of the module in Z
	int nx;		  // index in module in X
	int mz;		  // index in module in Z
};
// Volume/Detector structure
struct volume
{				 // parallelepiped
	plain fc[6]; // 6 plane definition (check indexing!)
				 // indexing- pairing oposite side	[0],[1]: -Z,+Z
				 //					[2],[3]: -X,+X
				 //					[4],[5]: -Y,+Y
	short subtype; // what detector has the same type of crystal array for detection efficies
				   // this was used in prev version (check Is it true here?)
	//	int mat;		//index of the material in the list
	// delete if "materials" NOT used!!!
	vec m0n0;		// absolute coordinates of the crystall (0,0)
	vec translateX; // two orthogonal vectors to replicate crystals	in x and y
	vec translateZ;
	vec center; // center of the module

	int detectorID; // ID in a list of detectors
	int sector;		// sector of the detectors of type ..., for detectors ONLY!!! type>=0
	int moduleX;	// module index in x direction
	int moduleZ;	// module index in z direction
};
// cross sections
struct crossSection
{
	float a[7][24]; // cross sections for one cystall theta(7) from 0 to 90 and Fi (24) from 0 to 360
};
// Detector and repetition structures (used for volumes generation only)
struct detector
{
	vec Base;  // detector system center in global coordinates
	vec detAx; // axis for rotational replication-sectors
			   //!!! only rotation around x, y or z!!!!!!!
	vec transl;		   // this is preliminary translation of original volume for following sectors replicas
	double startAngle; // start angle for rotational replica
	double pitchAngle; // pith angle for rotational replica

	int sectors; // total number of replicas

	double sectorHalfX; // half size of sector Do I need this?
	double sectorHalfY;
	double sectorHalfZ;
	// NOT-logical! subdivision of sector to the detection units - modules
	int modulesX; // number of modules in X direction!!!!
	int modulesZ; // number of modules in Z direction!!!!
	int gapX;	  // gap between modules in x-direction
	int gapZ;	  // gap between modules in z-direction

	double moduleX; // module dimensions - realy used only moduleZ
	double moduleY;
	double moduleZ;
	// crystal - logical division of module
	double crystalx; // x,y,z size of crystal
	double crystaly;
	double crystalz;
	int repeatX; // number of cystal in x direction
	int repeatZ; // number of cystal in z direction

	int volumes;
	int N;		 // total number of crystals in this detector
	int Nsector; // crystals in sector
	int Nmodule; // crystals in module

	vector<crossSection> cs511;
	vector<crossSection> cs450;
	// looks like I do not need this!!!!
	// int firstID;		//ID of the first crystal (global)
};
struct genRep
{
	float time;
	float angle;
	vec rotAxis;
	vec trans;
};
// scatter structures

struct scatPoint
{
	vec r; // coordinates
	int i; // index in a grid
	int j;
	int k;
};

// scatter Pathes from center ov volume to scatter point
struct scaterPath
{
	float ro_L;	 // line integral ro*dL
	float act_L; // line integral activity*dL
	float L;
	float ax; // direction from crystal to scatter point!!!
	float ay;
	float az;
	float solidAngle511;
	float solidAngle450;
};
//-----------Functions--------------------------------------------------------
// geometry functions
vec mkVec(double x, double y, double z); // create vector from its component
vec sum(vec a, vec b);
// sum of two vectors sum=a+b
vec sumK(vec a, vec b, double k);
// sum of two vectors. whis second vector multiplied by a scalar k
// sum=a+k*b
volume Translate(vec a, volume v0); // translate original volume v0 to new position,
									// translation vector is a
volume Rotate(vec n, double ang, volume v0); // rotation around one of the axis only!!!!
											 // n can be (1,0,0),(0,1,0) or (0,0,1) only!!!!
matrix3X3 RotationMatrix(vec n, double ang);
matrix3X3 mk3X3(double a00, double a01, double a02, // create matrix 3X3;
				double a10, double a11, double a12,
				double a20, double a21, double a22);
vec MatrixXVec(matrix3X3 M, vec a); // multiplay matrix and vector
vec vecToScal(vec a, double k);		// multiplication of vector and scalar
vec vecProd(vec a, vec b);			// vector product






#pragma pack(pop)
#endif