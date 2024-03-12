#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <algorithm>
#include <time.h>
#include <stdlib.h>
#include "SSS_data_structure.h"
#include "GATE_data_structure.h"
#include "ImageArray.h"
#include "PET_data.h"
// ... (other includes and definitions)

#define SGN(a) (((a) < 0) ? -1 : 1)
#define scalProd(a, b) (a.x * b.x + a.y * b.y + a.z * b.z)
#define MIN(a, b) (((a) < (b)) ? a : b)
#define crystals_in_ring 760

// for ran1
#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836

// for both ran1 and ran2
#define NTAB 32
#define NDIV (1 + (IM - 1) / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

// for ran2

#define IM1 2147483563
#define IM2 2147483399
#define AM_2 (1.0 / IM1)
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791

class Single_Scatter_Simulation
{
public:
    Single_Scatter_Simulation();
    ~Single_Scatter_Simulation();
    // Configure the simulation parameters
    void setSimulationParameters(ImageArray<float>& Emission_Image,ImageArray<float>& Attenuation_Image);
    void runSimulation();
    void runSimulationOnGPU();
    void setEmissionAndAttImg(ImageArray<float>& Emission_Image,ImageArray<float>& Attenuation_Image);
    void setPromptData(std::vector<PET_LST_event>& data);
    void setProtocol(PET_protocol_type protocol);    
    void setImgPosition(float imgPosition);    
    void setChunkSize(float chunkSize);
    void setOutputFilename(string outputFilename);
    void setDataPath(string datapath);
    //if static no partition; if Step and shoot, partion based on bed position; if CBM, partion every 30mm
    void partition_data(); 
private:
    string _datapath;
    int TotalCrystals = 0; // here it is global!!!! //should be done befor path generation!!!!
    // Image dimension
    int Ddimx; // number of voxel on x direction
    int Ddimy; // number of voxel on y direction
    int Ddimz; // number of voxel on z direction
    int Dstep = 30; // voxel size in mm
    float Dorigine[3] = {-299.5, -299.5, -131};
    double profile[760][760];
    double ScaleCompt = 11300; // scaling coefficient to noralize ComptonCorr output
                               // should be calculated or measured
    double Ethreshold = 450;   // default. It should be taken from detector file! check!
    float _imgPosition;
	PET_protocol_type _protocol=STATIC; //static by default
    float _chunkSize = 30; //default:30mm
    // Member variables
    std::vector<volume> V;
    std::vector<PET_LST_event>* promptData;
    //start index: [], left close, right close
    std::vector<int> _startIndex;
    std::vector<int> _endIndex;
    std::vector<int> _bedposition; //bed position in mm
    std::vector<detector> ListOfDetectors;
    std::vector<genRep> GenericRepetition;
    std::vector<scatPoint> ScatPoints;
    scaterPath *PathesMatrix;
    //float Dimg[Ddimx][Ddimy][Ddimz];
    //float Dmumap[Ddimx][Ddimy][Ddimz];
    ImageArray<float> Dimg;
    ImageArray<float> Dmumap;
    string outputFilename="scatter.dat"; //default name
    long seed, seed1;
    long *mypointer;
    long *mypointer1;
    // Member functions
    bool GenerateDetectors(int shiftBedPosition);
    bool GenerateScatterPoints();
    bool GeneratePathes();
    bool FillCrossSections();
    double ComptonCorr(int crystal1, int crystal2, double TOF_dist);

    cryst CrystPosFromID(int crystalID);
    int CrystalInd(cryst IN);

    // Random number generators
    float ran1(long *idum);
    float ran2(long *idum);
    float gasdev(long *idum);

    // some constants used in calculations-----------------------------------
    // HUGE should be larger than size of the system...
    // TINY any reasonable small value
    double const HUGEv = 1.0e50;
    double const TINYv = 1.0e-5;
    // attenuation coeff for 511keV
    double const mwater511 = 9.55; // mm*mm/g  //0.00955;	//(1/mm)?
    double roWater = 0.001;        // g/mm^3
    // 0.001(g/mm^3)* mu(mm^2/g)=mwater511(1/mm)
    double MuToRo = 1. / mwater511; // 0.001/mwater511;//g/mm^2

    double const mLSO511 = 0.085;

    double const PI = 3.1415926535897932384626433832795;
    double const twoPI = 2 * PI;
    double const halfPI = PI / 2.;
    double const quaterPI = PI / 4.;
    double const angleStep = PI / 12.;
    double const HalfAngleStep = angleStep / 2.;
    double const GradInRad = PI / 180.; // step in theta and fi angles in precalculated crystal crass sections

    //////////////////////////////////////////
    float timeFWHM = 150;                                // psec/////////////////////////////////////////input value!!!!!!!!
    float timeBinSize = 50;                              // in psec
                                                         // this define the Inegration step!!!
    float speedOfLight = 0.3;                            // in mm/psec
    float timeSigma = timeFWHM / 2.355;                  // ps
    float distBinSize = timeBinSize * speedOfLight / 2.; // mm
    float distSigma = timeSigma * speedOfLight / 2.;     // mm
    float halfDistSize = -distSigma * (int)(distBinSize / 2);
    //////////////////////////////////////////
};
