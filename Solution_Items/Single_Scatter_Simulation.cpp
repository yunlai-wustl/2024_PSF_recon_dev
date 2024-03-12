#include "Single_Scatter_Simulation.h"
#include "SSS_data_structure.h"
#include <time.h>
Single_Scatter_Simulation::Single_Scatter_Simulation()
{
}
Single_Scatter_Simulation::~Single_Scatter_Simulation()
{
}

void Single_Scatter_Simulation::setPromptData(std::vector<PET_LST_event> &data)
{
	promptData = &data;
}
void Single_Scatter_Simulation::setProtocol(PET_protocol_type protocol)
{
	_protocol = protocol;
}

void Single_Scatter_Simulation::setChunkSize(float chunkSize)
{
	_chunkSize = chunkSize;
}

void Single_Scatter_Simulation::setImgPosition(float imgPosition)
{
	_imgPosition = imgPosition;
}
void Single_Scatter_Simulation::setOutputFilename(string filename)
{
	outputFilename = filename;
}
void Single_Scatter_Simulation::setDataPath(string datapath){
	_datapath = datapath;
}


void Single_Scatter_Simulation::partition_data()
{
	if (_protocol == STATIC)
	{
		std::cout << "Static dataset" << endl;
		_startIndex.push_back(_imgPosition);
		_endIndex.push_back(promptData->size());
		_bedposition.push_back(_imgPosition);
	}
	if (_protocol == STEP_AND_SHOOT)
	{
		std::cout << "Step and shoot dataset" << endl;
		int currentBedPosition;
		int ind = 0;
		_startIndex.push_back(0);
		// note that the first several lines of data have bed position = 0;
		while ((*promptData)[ind].t0 == 0.00)
			ind++;
		currentBedPosition = int((*promptData)[ind].bed_position);
		_bedposition.push_back(currentBedPosition);

		for (int i = ind; i < promptData->size(); i++)
		{
			if (int((*promptData)[i].bed_position) != currentBedPosition)
			{
				ind = i;
				while ((*promptData)[ind].t0 == 0.00)
					ind++;
				_endIndex.push_back(i - 1);
				currentBedPosition = int((*promptData)[ind].bed_position);
				_bedposition.push_back(currentBedPosition);
				_startIndex.push_back(i);
				i = ind;
			}
		}
		_endIndex.push_back(promptData->size() - 1);
	}
	if (_protocol == CBM)
	{
		std::cout << "CBM dataset" << endl;
		float currentBedPosition;
		int ind = 0;
		while ((*promptData)[ind].t0 == 0.00)
			ind++;
		currentBedPosition = (*promptData)[ind].bed_position;

		float currentMax = currentBedPosition;

		size_t startIndex = ind;

		while (startIndex < promptData->size())
		{
			currentMax += _chunkSize;
			auto it = std::upper_bound(promptData->begin() + startIndex, promptData->end(), currentMax,
									   [](const float &val, const PET_LST_event &event)
									   {
										   return val < event.bed_position;
									   });
			size_t endIndex = it - promptData->begin() - 1;
			startIndex == ind ? _startIndex.push_back(0) : _startIndex.push_back(startIndex);
			_endIndex.push_back(endIndex);
			_bedposition.push_back(((*promptData)[startIndex].bed_position + (*promptData)[endIndex].bed_position) / 2);
			startIndex = endIndex + 1;
		}
	}
	for (int i = 0; i < _startIndex.size(); i++)
	{
		std::cout << "group " << i << ", startIndex:" << _startIndex[i] << ", endIndex:" << _endIndex[i]
			 << ", current group bed position:" << _bedposition[i] << endl;
	}
}

void Single_Scatter_Simulation::setEmissionAndAttImg(ImageArray<float> &Emission_Image, ImageArray<float> &Attenuation_Image)
{	
	Ddimx = Emission_Image.getDimX();
	Ddimy = Emission_Image.getDimY();
	Ddimz = Emission_Image.getDimZ();
	Dimg.CopyFromImageArray(Emission_Image);
	Dmumap.CopyFromImageArray(Attenuation_Image);
	Dmumap.ScaledBy(MuToRo);
	Dimg.WriteToFile("Dimg.img");
	Dmumap.WriteToFile("Dmumap.img");
	std::cout << "Image Dimension for SSS:" << Ddimx << "," << Ddimy << "," << Ddimz << endl;
}

void Single_Scatter_Simulation::runSimulation()
{
	clock_t startTime,endTime;
	std::cout << sizeof(float) << endl;
	std::cout << "Size of scaterPath: " << sizeof(scaterPath) << " bytes" << endl;
	unsigned bunit = sizeof(scaterPath);
	// std::cout<<sizeof(scaterPath)<<" : "<<sizeof(scaterPath)*75801600.<<endl;
	std::cout << "Start" << endl;
	seed = -time(NULL); // 200;
	mypointer = &seed;
	seed1 = -200; //-time(NULL);//200;//not used
	std::cout << "Start seed " << seed;
	std::cout << "; Start seed1 " << seed1 << endl;
	mypointer1 = &seed1;
	std::cout << "grouping data ..." << endl;
	startTime = clock();
	partition_data();
	endTime = clock();
	std::cout<<"grouping data took "<<(endTime-startTime)/CLOCKS_PER_SEC<<"s"<<endl;
	if (_protocol == STATIC)
	{
		startTime = clock();
		bool flag = GenerateDetectors(0);
		endTime = clock();
		std::cout<<"generating detectors took "<<(endTime-startTime)/CLOCKS_PER_SEC<<"s"<<endl;

		if (!flag)
		{
			std::cout << "Detectors are not created!" << endl;
			return;
		}
		// calculate total number of crystall
		for (int i = 0; i < ListOfDetectors.size(); i++)
		{
			TotalCrystals += ListOfDetectors[i].N;
			std::cout << i << " " << ListOfDetectors[i].N << endl;
		}
		std::cout << "TotalCrystals= " << TotalCrystals << endl;
		std::cout << "Detectors " << ListOfDetectors.size() << endl;
		std::cout << Ethreshold << " keV" << endl;
		// system display: cycle can be deleted
		for (int i = 0; i < ListOfDetectors.size(); i++)
		{
			std::cout << "detector " << i << endl;
			std::cout << "modulesX x modulesZ " << ListOfDetectors[i].modulesX << "x" << ListOfDetectors[i].modulesZ << endl;
			std::cout << "repeatX x repeatZ " << ListOfDetectors[i].repeatX << "x" << ListOfDetectors[i].repeatZ << endl;
			std::cout << "volumes " << ListOfDetectors[i].volumes << endl;
			std::cout << "sectors " << ListOfDetectors[i].sectors << endl;
			std::cout << "crystals " << ListOfDetectors[i].N << endl;
			std::cout << "crystals in sector " << ListOfDetectors[i].Nsector << endl;
			std::cout << "crystals in module " << ListOfDetectors[i].Nmodule << endl
				 << endl;
		}
		startTime = clock();
		flag = GenerateScatterPoints(); // Generate scatter points and scatter pathes
		endTime = clock();
		std::cout<<"generating scatter points took "<<(endTime-startTime)/CLOCKS_PER_SEC<<"s"<<endl;

		if (!flag)
		{
			std::cout << "Scatter points are not created!" << endl;
			return;
		}
		else
			std::cout << flag << " scatter poins are generated!" << endl;
		
		startTime = clock();
		flag = FillCrossSections(); // fill cross sections from outside files
		endTime = clock();
		std::cout<<"filling cross sections took "<<(endTime-startTime)/CLOCKS_PER_SEC<<"s"<<endl;

		if (!flag)
		{
			std::cout << "Cross sections are not created!" << endl;
			return;
		}
		else
			std::cout << "Cross sections are IN" << endl;
		startTime = clock();
		flag = GeneratePathes(); // Generate scatter points and scatter pathes
		endTime = clock();
		std::cout<<"generating pathes took "<<(endTime-startTime)/CLOCKS_PER_SEC<<"s"<<endl;

		if (!flag)
		{
			std::cout << "Scatter Pathes are not created!" << endl;
			return;
		}
		std::cout << "Here" << endl;

		// Read Files
		int totalEvent = 0;
		char cname[100];

		vector<double> scatterFraction(promptData->size(), 0.0);
		
		startTime = clock();
		for (int i = 0; i < promptData->size(); i++)
		{
			auto event = (*promptData)[i];
			scatterFraction[i] = ComptonCorr(event.src_id, event.dest_id, event.TOF_dist);
		}
		endTime = clock();
		std::cout<<"calculating scatter fractions took "<<(endTime-startTime)/CLOCKS_PER_SEC<<"s"<<endl;

		ofstream outFile(outputFilename.c_str(), ios::binary);
		if (!outFile)
		{
			throw runtime_error("Unable to open file ");
		}
		// Write the entire vector to the file in one go
		if (!scatterFraction.empty())
		{
			outFile.write(reinterpret_cast<const char *>(scatterFraction.data()), scatterFraction.size() * sizeof(double));
		}
		outFile.close();
	}
	else
	{
		double *scatterFraction = (double *)malloc(promptData->size() * sizeof(double));
		int nChunk = _bedposition.size();
		for (int i = 0; i < nChunk; i++)
		{
			std::cout << "calculating scatter fraction for chunk " << i << endl;
			int detectorShiftDistance = _bedposition[i];
			startTime = clock();
			bool flag = GenerateDetectors(detectorShiftDistance);
			endTime = clock();
			std::cout<<"GenerateDetectors took "<<(endTime-startTime)/CLOCKS_PER_SEC<<"s"<<endl;

			std::cout<<"chunk "<<i<<", shift the whole system by "<<detectorShiftDistance<<"mm"<<endl;
			if (!flag)
			{
				std::cout << "Detectors are not created!" << endl;
				return;
			}
			// calculate total number of crystall
			for (int i = 0; i < ListOfDetectors.size(); i++)
			{
				TotalCrystals += ListOfDetectors[i].N;
				std::cout << i << " " << ListOfDetectors[i].N << endl;
			}
			std::cout << "TotalCrystals= " << TotalCrystals << endl;
			std::cout << "Detectors " << ListOfDetectors.size() << endl;
			std::cout << Ethreshold << " keV" << endl;
			// system display: cycle can be deleted
			for (int i = 0; i < ListOfDetectors.size(); i++)
			{
				std::cout << "detector " << i << endl;
				std::cout << "modulesX x modulesZ " << ListOfDetectors[i].modulesX << "x" << ListOfDetectors[i].modulesZ << endl;
				std::cout << "repeatX x repeatZ " << ListOfDetectors[i].repeatX << "x" << ListOfDetectors[i].repeatZ << endl;
				std::cout << "volumes " << ListOfDetectors[i].volumes << endl;
				std::cout << "sectors " << ListOfDetectors[i].sectors << endl;
				std::cout << "crystals " << ListOfDetectors[i].N << endl;
				std::cout << "crystals in sector " << ListOfDetectors[i].Nsector << endl;
				std::cout << "crystals in module " << ListOfDetectors[i].Nmodule << endl
					 << endl;
			}
	
			startTime = clock();
			flag = GenerateScatterPoints(); // Generate scatter points and scatter pathes
			if (!flag)
			{
				std::cout << "Scatter points are not created!" << endl;
				return;
			}
			else
				std::cout << flag << " scatter poins are generated!" << endl;
			endTime = clock();
			std::cout<<"GenerateScatterPoints took "<<(endTime-startTime)/CLOCKS_PER_SEC<<"s"<<endl;

			startTime = clock();
			flag = FillCrossSections(); // fill cross sections from outside files
			if (!flag)
			{
				std::cout << "Cross sections are not created!" << endl;
				return;
			}
			else
				std::cout << "Cross sections are IN" << endl;
			endTime = clock();
			std::cout<<"FillCrossSections took "<<(endTime-startTime)/CLOCKS_PER_SEC<<"s"<<endl;

			startTime = clock();
			flag = GeneratePathes(); // Generate scatter points and scatter pathes
			if (!flag)
			{
				std::cout << "Scatter Pathes are not created!" << endl;
				return;
			}
			endTime = clock();
			std::cout<<"GeneratePathes took "<<(endTime-startTime)/CLOCKS_PER_SEC<<"s"<<endl;

			startTime = clock();
			for (int j = _startIndex[i]; j <= _endIndex[i]; j++)
			{
				auto event = (*promptData)[j];
				scatterFraction[j] = ComptonCorr(event.src_id, event.dest_id, event.TOF_dist);
			}
			endTime = clock();
			std::cout<<"calculating scatter fractions took "<<(endTime-startTime)/CLOCKS_PER_SEC<<"s"<<endl;

			std::cout<<"releasing memory"<<endl;
			free(PathesMatrix);
			ListOfDetectors.clear(); //clear detectors from last run
			ScatPoints.clear();
			GenericRepetition.clear();
			V.clear();
			TotalCrystals = 0;
		}
		FILE *file = fopen(outputFilename.c_str(), "wb");
		if (file == NULL)
		{
			perror("Error opening file");
			free(scatterFraction);
			return;
		}

		// Write the entire array to the file in one operation
		size_t itemsWritten = fwrite(scatterFraction, sizeof(double), promptData->size(), file);
		if (itemsWritten != promptData->size())
		{
			perror("Error writing to file");
		}

		// Close the file
		fclose(file);

		// Free the allocated memory
		free(scatterFraction);
	}

	return;
}

void Single_Scatter_Simulation::runSimulationOnGPU()
{

	std::cout << sizeof(float) << endl;
	std::cout << "Size of scaterPath: " << sizeof(scaterPath) << " bytes" << endl;
	unsigned bunit = sizeof(scaterPath);
	// std::cout<<sizeof(scaterPath)<<" : "<<sizeof(scaterPath)*75801600.<<endl;
	std::cout << "Start" << endl;
	seed = -time(NULL); // 200;
	mypointer = &seed;
	seed1 = -200; //-time(NULL);//200;//not used
	std::cout << "Start seed " << seed;
	std::cout << "; Start seed1 " << seed1 << endl;
	mypointer1 = &seed1;
	std::cout << "grouping data ..." << endl;
	partition_data();
	if (_protocol == STATIC)
	{
		bool flag = GenerateDetectors(0);
		if (!flag)
		{
			std::cout << "Detectors are not created!" << endl;
			return;
		}
		// calculate total number of crystall
		for (int i = 0; i < ListOfDetectors.size(); i++)
		{
			TotalCrystals += ListOfDetectors[i].N;
			std::cout << i << " " << ListOfDetectors[i].N << endl;
		}
		std::cout << "TotalCrystals= " << TotalCrystals << endl;
		std::cout << "Detectors " << ListOfDetectors.size() << endl;
		std::cout << Ethreshold << " keV" << endl;
		// system display: cycle can be deleted
		for (int i = 0; i < ListOfDetectors.size(); i++)
		{
			std::cout << "detector " << i << endl;
			std::cout << "modulesX x modulesZ " << ListOfDetectors[i].modulesX << "x" << ListOfDetectors[i].modulesZ << endl;
			std::cout << "repeatX x repeatZ " << ListOfDetectors[i].repeatX << "x" << ListOfDetectors[i].repeatZ << endl;
			std::cout << "volumes " << ListOfDetectors[i].volumes << endl;
			std::cout << "sectors " << ListOfDetectors[i].sectors << endl;
			std::cout << "crystals " << ListOfDetectors[i].N << endl;
			std::cout << "crystals in sector " << ListOfDetectors[i].Nsector << endl;
			std::cout << "crystals in module " << ListOfDetectors[i].Nmodule << endl
				 << endl;
		}

		flag = GenerateScatterPoints(); // Generate scatter points and scatter pathes
		if (!flag)
		{
			std::cout << "Scatter points are not created!" << endl;
			return;
		}
		else
			std::cout << flag << " scatter poins are generated!" << endl;
		flag = FillCrossSections(); // fill cross sections from outside files
		if (!flag)
		{
			std::cout << "Cross sections are not created!" << endl;
			return;
		}
		else
			std::cout << "Cross sections are IN" << endl;

		flag = GeneratePathes(); // Generate scatter points and scatter pathes
		if (!flag)
		{
			std::cout << "Scatter Pathes are not created!" << endl;
			return;
		}
		std::cout << "Here" << endl;

		// Read Files
		int totalEvent = 0;
		char cname[100];

		vector<double> scatterFraction(promptData->size(), 0.0);
		for (int i = 0; i < promptData->size(); i++)
		{
			auto event = (*promptData)[i];
			scatterFraction[i] = ComptonCorr(event.src_id, event.dest_id, event.TOF_dist);
		}
		ofstream outFile(outputFilename.c_str(), ios::binary);
		if (!outFile)
		{
			throw runtime_error("Unable to open file ");
		}
		// Write the entire vector to the file in one go
		if (!scatterFraction.empty())
		{
			outFile.write(reinterpret_cast<const char *>(scatterFraction.data()), scatterFraction.size() * sizeof(double));
		}
		outFile.close();
	}
	/*
	ifstream dataFile;
	ofstream dataOut;
	mini_coinc_event event;
	for (int f = 0; f < 1; f++)
	{
		dataFile.open("prompt_data.lst", ios::in | ios::binary);
		if (!dataFile)
		{
			std::cout << "File prompt_data.lst do not exists! continue" << endl;
			continue;
		}
		std::cout << "File to read prompt_data.lst open!" << endl;
		// file to write//////////////////////////////////////////////////////////////////
		// sprintf(cname,"Out_eltub_prompt_SS_%dmin.dat",5);
		// dataOut.open(cname, ios::out | ios::binary);
		// was		dataOut.open("Out_smallcylprompt_SS.dat", ios::out | ios::binary);
		dataOut.open("scatter_sergey.dat", ios::out | ios::binary);
		if (!dataOut)
		{
			std::cout << "File scatter.dat do not open! continue" << endl;
			continue;
		}
		std::cout << "File to write scatter.dat open!" << endl;
		unsigned int counter = 0;
		do
		{
			dataFile.read((char *)&event, sizeof(mini_coinc_event)); // read data
			counter++;
			if (!(counter % 1000000000))
				std::cout << counter << " ";
			int DTime = (double)event.diff_time * 1.e12; // time inPS
			double ComCor = ComptonCorr(event.crystal_index_1, event.crystal_index_2, DTime);
			dataOut.write((char *)&ComCor, sizeof(double));
		} while (!dataFile.eof());
		dataFile.close();
		dataOut.close();
		std::cout << endl;
		std::cout << "Total events in this file=" << counter << endl;
	}
	*/
	return;
}

bool Single_Scatter_Simulation::GenerateDetectors(int shiftBedPosition)
{
	std::cout << ListOfDetectors.size() << " first detector ID" << endl;
	// read the detector structure------------------------------------------------------------
	ifstream inDet;
	inDet.open(_datapath+"detectorFromMac.txt");
	if (!inDet)
	{
		std::cout << "detectorFromMac.txt is not found! Exit!" << endl;
		return 0;
	}
	// line 1
	inDet >> Ethreshold;
	inDet.ignore(100, '&');
	std::cout << Ethreshold << endl;
	// line2
	string NewDetector;
	inDet >> NewDetector;
	if (NewDetector != "NEW")
	{
		std::cout << "This is not detectors file! Exit!" << endl;
		return 0;
	}
	detector dettmp;
	do
	{
		inDet.ignore(100, '&');
		double X0, Y0, Z0;
		// line3 - read the coordinates of detector center(==rotation axis origin)
		inDet >> X0 >> Y0 >> Z0;
		inDet.ignore(100, '&');
		vec Base = mkVec(X0, Y0, Z0);
		dettmp.Base = Base;
		int ax, ay, az;
		// line4 - read the axis of the rotation
		inDet >> ax >> ay >> az;
		inDet.ignore(100, '&');
		if ((ax * ax + ay * ay + az * az) != 1 || ax < 0 || ay < 0 || az < 0) //??? may be consider (-1) for negative rotation. Do I need this??
		{
			std::cout << "Wrong rotation axis! Can be only (100)(010)(001)! You have (" << ax << ay << az << ")" << endl;
			return 0;
		}
		vec n = mkVec(ax, ay, az); // axis direction
		dettmp.detAx = n;
		// line 5 - read the translation vector from X0,Y0,Z0
		// This will be rotation radius
		// can be any but for clarity it will be better to move only in (-Y)
		//  shift in X and Z will be better to do before: see "Detector center X0,Y0,Z0
		inDet >> X0 >> Y0 >> Z0;
		Z0 = Z0+shiftBedPosition;
		inDet.ignore(100, '&');
		vec Translation = mkVec(X0, Y0, Z0);
		dettmp.transl = Translation;
		// lines 6,7 and 8: crystal size h_x, h_y, h_z
		double h_y;
		inDet >> h_y;
		inDet.ignore(100, '&');
		double h_x;
		inDet >> h_x;
		inDet.ignore(100, '&');
		double h_z;
		inDet >> h_z;
		inDet.ignore(100, '&');
		dettmp.crystalx = h_x;
		dettmp.crystaly = h_y;
		dettmp.crystalz = h_z;
		// lines 9 and 10: number of the array crystals in x and z derections
		// there is no crystal replica in Y-direction
		int repeatX, repeatZ;
		inDet >> repeatX;
		inDet.ignore(100, '&');
		inDet >> repeatZ;
		inDet.ignore(100, '&');
		dettmp.repeatX = repeatX;
		dettmp.repeatZ = repeatZ;
		// MODULES IN SECTOR
		int modulesX, modulesZ;
		double gapX, gapZ;
		// lines 11 and 12: modules in X and gap between modules in this direction
		inDet >> modulesX;
		inDet.ignore(100, '&');
		inDet >> gapX;
		inDet.ignore(100, '&');
		dettmp.modulesX = modulesX;
		dettmp.gapX = gapX;
		// lines 13 and 14: modules in Z and gap between modules in this direction
		inDet >> modulesZ;
		inDet.ignore(100, '&');
		inDet >> gapZ;
		inDet.ignore(100, '&');
		dettmp.modulesZ = modulesZ;
		dettmp.gapZ = gapZ;
		// Repeater
		string Repeater;
		inDet >> Repeater;
		inDet.ignore(100, '&');

		int repeatR;
		double pitchAngle;
		double StartAngle;
		//	dettmp.firstID=0;
		//	for(int n=0;n<ListOfDetectors.size();n++)
		//	{
		//		dettmp.firstID+=ListOfDetectors[n].N;
		//	}
		if (Repeater == "ring")
		{
			// dettmp.firstID=0;
			// for(int n=0;n<ListOfDetectors.size();n++)
			//{
			//	dettmp.firstID+=ListOfDetectors[n].N;
			// }
			// SECTORS
			// line 15 - read start angle (in grad) for sector replica
			inDet >> StartAngle;
			inDet.ignore(100, '&');
			StartAngle = StartAngle * GradInRad; // transform to radian
			dettmp.startAngle = StartAngle;
			// line 16 - read span angle (in grad) for whole sector replicas
			double SpanAngle;
			inDet >> SpanAngle;
			inDet.ignore(100, '&');
			SpanAngle = SpanAngle * GradInRad; // transform to radian
			// line 17 - number of sectors
			inDet >> repeatR;
			inDet.ignore(100, '&');
			dettmp.sectors = repeatR;
			// here was a problem in old codes? check it was so!!!!!
			if (repeatR > 1)
				pitchAngle = SpanAngle / (repeatR - 1);
			else
				pitchAngle = 0;
			// possibly I do not need this - delete this structure element if true
			dettmp.N = (dettmp.repeatX * dettmp.repeatZ) * (dettmp.modulesX * dettmp.modulesZ) * dettmp.sectors;
			//	crystals in module	x	modules			x	sectors
		}
		else
		{
			if (Repeater != "genericRepeater")
			{
				std::cout << "Wrong Repeater " << Repeater << endl;
				return 0;
			}
			int sectors;
			inDet >> sectors;
			inDet.ignore(100, '&');
			dettmp.sectors = sectors;
			// possibly I do not need this - delete this structure element if true
			// the same as above -see prev coment
			dettmp.N = (dettmp.repeatX * dettmp.repeatZ) * (dettmp.modulesX * dettmp.modulesZ) * dettmp.sectors;
			for (int s = 0; s < sectors; s++)
			{
				genRep genTMP;
				inDet >> genTMP.angle;
				inDet >> genTMP.rotAxis.x >> genTMP.rotAxis.y >> genTMP.rotAxis.z;
				inDet >> genTMP.trans.x >> genTMP.trans.y >> genTMP.trans.z;
				inDet.ignore(100, '&');
				GenericRepetition.push_back(genTMP);
			}
		}
		// MATERIAL not used reed blanks
		// lines 18,19,20
		// name: one word. !!! Use the same name for identical materials!!!
		// dencity in g/sm3
		// number of elements that compose "material" - used for loop over all individual elements
		string TMP;
		inDet >> TMP;
		inDet.ignore(100, '&');
		float Density;
		inDet >> Density;
		inDet.ignore(100, '&');
		int maxElements;
		inDet >> maxElements;
		inDet.ignore(100, '&');
		// read elements maxElements time every line is: Z and number of atoms
		for (int i = 0; i < maxElements; i++)
		{
			float zz;
			float nn;
			// set Z and n for element
			inDet >> zz >> nn;
			inDet.ignore(100, '&');
		}
		// pause in reading to create volumes...
		// Attention!!! HERE is a change!!! modules are separate volumes!!!!
		// half sizes of volume/array;
		double Xmax = repeatX * h_x / 2.;
		double Ymax = h_y / 2.;
		double Zmax = repeatZ * h_z / 2.; // NO replica for module in other directions exept X and Z

		double stepModX = repeatX * h_x + gapX;
		double stepModZ = repeatZ * h_z + gapZ;

		dettmp.sectorHalfX = (repeatX * h_x * modulesX + gapX * (modulesX - 1)) / 2.;
		dettmp.sectorHalfY = Ymax;
		dettmp.sectorHalfZ = (repeatZ * h_z * modulesZ + gapZ * (modulesZ - 1)) / 2.;

		dettmp.moduleX = repeatX * h_x;
		dettmp.moduleY = h_y;
		dettmp.moduleZ = repeatZ * h_z;
		// new addition!!!!
		dettmp.Nmodule = repeatX * repeatZ;
		dettmp.Nsector = dettmp.Nmodule * modulesX * modulesZ;
		dettmp.volumes = modulesX * modulesZ * dettmp.sectors;
		//		dettmp.N //already calculated

		ListOfDetectors.push_back(dettmp);
		// base volume
		volume VI;
		VI.fc[0].n.x = 0.;
		VI.fc[0].n.y = 0.;
		VI.fc[0].n.z = -1.;
		VI.fc[0].p = -Zmax;
		VI.fc[1].n.x = 0.;
		VI.fc[1].n.y = 0.;
		VI.fc[1].n.z = 1.;
		VI.fc[1].p = -Zmax;
		VI.fc[2].n.x = -1.;
		VI.fc[2].n.y = 0.;
		VI.fc[2].n.z = 0.;
		VI.fc[2].p = -Xmax;
		VI.fc[3].n.x = 1.;
		VI.fc[3].n.y = 0.;
		VI.fc[3].n.z = 0.;
		VI.fc[3].p = -Xmax;
		VI.fc[4].n.x = 0.;
		VI.fc[4].n.y = -1.;
		VI.fc[4].n.z = 0.;
		VI.fc[4].p = -Ymax;
		VI.fc[5].n.x = 0.;
		VI.fc[5].n.y = 1.;
		VI.fc[5].n.z = 0.;
		VI.fc[5].p = -Ymax;
		VI.m0n0 = mkVec(-Xmax + h_x / 2., 0, -Zmax + h_z / 2.);
		VI.center = mkVec(0, 0, 0);

		VI.translateX = mkVec(h_x, 0, 0);
		VI.translateZ = mkVec(0, 0, h_z);

		VI.detectorID = ListOfDetectors.size() - 1;

		VI = Translate(Translation, VI);
		// translate to the module (0,0) position, than move to the needed position
		vec translateTo0_0 = mkVec(-dettmp.sectorHalfX + Xmax, 0, -dettmp.sectorHalfZ + Zmax);

		VI = Translate(translateTo0_0, VI);
		// replication part
		if (Repeater == "ring")
		{
			for (int s = 0; s < repeatR; s++)
			{
				double ang = StartAngle + s * pitchAngle;

				for (int mz = 0; mz < dettmp.modulesZ; mz++)
				{
					for (int nx = 0; nx < dettmp.modulesX; nx++)
					{
						// translate to nx,mz position
						vec translateNM = mkVec(nx * stepModX, 0, mz * stepModZ);
						volume tmp = Translate(translateNM, VI);

						tmp = Rotate(n, ang, tmp);
						tmp.sector = s;
						tmp.moduleX = nx;
						tmp.moduleZ = mz;
						// if translate is not 0 translate the whole detector structure
						if (Base.x != 0 || Base.y != 0 || Base.z != 0)
							tmp = Translate(Base, tmp);
						V.push_back(tmp);
					}
				}
			}
		}
		if (Repeater == "genericRepeater")
		{
			for (int ss = 0; ss < GenericRepetition.size(); ss++)
			{
				for (int mz = 0; mz < dettmp.modulesZ; mz++)
				{
					for (int nx = 0; nx < dettmp.modulesX; nx++)
					{
						// translate to nx,mz position
						vec translateNM = mkVec(nx * stepModX, 0, mz * stepModZ);
						volume tmp = Translate(translateNM, VI);
						tmp = Rotate(GenericRepetition[ss].rotAxis, GenericRepetition[ss].angle * GradInRad, tmp);
						tmp.sector = ss;
						tmp.moduleX = nx;
						tmp.moduleZ = mz;
						tmp = Translate(GenericRepetition[ss].trans, tmp);
						// if translate is not 0 translate the whole detector structure
						if (Base.x != 0 || Base.y != 0 || Base.z != 0)
							tmp = Translate(Base, tmp);
						V.push_back(tmp);
					}
				}
			}
		}
		GenericRepetition.clear();
		inDet >> NewDetector;
		std::cout << "new " << NewDetector << endl;
	} while (NewDetector == "NEW");
	inDet.close();

	std::cout << ListOfDetectors.size() << " last" << endl;
	std::cout << V.size() << " total number of volumes" << endl;

	if (ListOfDetectors.size() && V.size())
		return true;
	else
		return false;
}

// Generate scatter points and scatter pathes
bool Single_Scatter_Simulation::GenerateScatterPoints()
{
	// TFile *ff;
	// ff = new TFile("scatterPointsV2.root", "RECREATE");
	// TH2F *scatPoints = new TH2F("scatPoints", "scatPoints", 600, Dorigine[0], Dorigine[0] + 600, 600, Dorigine[1], Dorigine[1] + 600);
	for (int k = 0; k < Ddimz; k++)
	{
		for (int j = 0; j < Ddimy; j++)
		{
			for (int i = 0; i < Ddimx; i++)
			{
				if (Dmumap(i, j, k) > 0) // only for non empty voxels
				{
					scatPoint tmp;
					tmp.i = i;
					tmp.j = j;
					tmp.k = k;
					tmp.r.z = Dorigine[2] + k * Dstep + ran1(mypointer) * Dstep;
					tmp.r.x = Dorigine[0] + i * Dstep + ran1(mypointer) * Dstep; // it is supposed that there ia no border events;
					tmp.r.y = Dorigine[1] + j * Dstep + ran1(mypointer) * Dstep;
					ScatPoints.push_back(tmp);
				}
			}
		}
	}
	std::cout << ScatPoints.size() << " scatter poins are generated!!!" << endl;
	if (ScatPoints.size())
		return true;
	else
		return false;
}

bool Single_Scatter_Simulation::FillCrossSections()
{
	int E = 511;
	float aToRead[168]; // aToRead[7][24];
	crossSection tmp;
	for (int e = 0; e < 2; e++)
	{
		if (e)
			E = Ethreshold;
		for (int dd = 0; dd < ListOfDetectors.size(); dd++)
		{
			std::cout << "Start detector " << dd << " for energy " << E << endl;
			char cname[100];
			sprintf(cname, "crossSection_%d_%d.dat", dd, E);
			ifstream crosSection;
			crosSection.open(_datapath+cname, ios::in | ios::binary);
			if (!crosSection)
			{
				std::cout << " File " << cname << " do not open! Exit!" << endl;
				return 0;
			}
			else
				std::cout << "File " << cname << " open." << endl;
			// read data and fill crossection vector
			int totalCrystals = ListOfDetectors[dd].repeatX * ListOfDetectors[dd].repeatZ;
			for (int n = 0; n < totalCrystals; n++)
			{
				crosSection.read((char *)aToRead, sizeof(aToRead));
				// why tmp.a=aToRead doesn't work???
				// try this
				// for(int t=0;t<7;t++) for(int f=0;f<24;f++)
				for (int ind = 0; ind < 168; ind++)
				{
					// int ind=t+f*7;
					int t = ind % 7;
					int f = ind / 7;
					tmp.a[t][f] = aToRead[ind];
				}
				if (e == 0)
					ListOfDetectors[dd].cs511.push_back(tmp);
				if (e == 1)
					ListOfDetectors[dd].cs450.push_back(tmp);
			}
			crosSection.close();
		}
	}
	return true;
}

bool Single_Scatter_Simulation::GeneratePathes() // crystals-scatter
{
	std::cout << "Generate pathes" << endl;
	vec g0 = mkVec(Dorigine[0], Dorigine[1], Dorigine[2]);
	vec gStep = mkVec(Dstep, Dstep, Dstep); // grid steps(deltax,deltay,deltaz)
	int max[3] = {Ddimx, Ddimy, Ddimz};

	int Pathes = ScatPoints.size() * TotalCrystals;
	std::cout << "number of scatter points:" << ScatPoints.size() << ", number of crystals:" << TotalCrystals << endl;
	PathesMatrix = new (nothrow) scaterPath[Pathes];
	if (PathesMatrix == 0)
	{
		std::cout << "Error: memory could not be allocated for PathesMatrix data!" << endl;
		// return false;
	}
	else
		std::cout << "Memory alocated for PathesMatrix:" << Pathes << "*" << sizeof(scaterPath) << " b" << endl;

	// loop over scatter points
	unsigned int NNN = 0;
	for (int n = 0; n < ScatPoints.size(); n++)
	{
		unsigned int NN = 0;
		vec r_n = ScatPoints[n].r;
		// loop over crystals
		for (int i = 0; i < TotalCrystals; i++)
		{
			// loop over crystals
			cryst Crystal = CrystPosFromID(i);
			int v = Crystal.V;
			int nc = Crystal.nx;
			int mc = Crystal.mz;
			int crystalInd = nc + mc * ListOfDetectors[V[v].detectorID].repeatX;

			scaterPath tmp;
			tmp.ro_L = 0;
			tmp.act_L = 0;

			// a is vector from crystal center (r_i) to scatter pont (r_n)
			vec r_i = V[v].m0n0;
			r_i = sumK(r_i, V[v].translateX, nc);
			r_i = sumK(r_i, V[v].translateZ, mc);
			vec a = sumK(r_n, r_i, -1);
			double norm2 = a.x * a.x + a.y * a.y + a.z * a.z;
			double norm = sqrt(norm2); // distance between points->put in a structure
			tmp.L = norm;
			a.x = a.x / norm;
			a.y = a.y / norm;
			a.z = a.z / norm;
			tmp.ax = a.x;
			tmp.ay = a.y;
			tmp.az = a.z;
			double cosTheta = scalProd(V[v].fc[5].n, a);
			if (cosTheta < 0.01) // theta<0.6deg Almost normal!!!
			{
				tmp.solidAngle511 = ListOfDetectors[V[v].detectorID].cs511[crystalInd].a[0][0] / norm2;
				tmp.solidAngle450 = ListOfDetectors[V[v].detectorID].cs450[crystalInd].a[0][0] / norm2;
			}
			else
			{
				// test//remove it later if it doesn'happen
				if (abs(cosTheta) > 1.)
				{
					std::cout << "abs(cosTheta)>1.;" << cosTheta << endl;
					cosTheta = 1;
				}
				double Theta = acos(cosTheta);
				int intTheta = Theta / angleStep; // theta index in sensitivity file
				int intThetaNext = intTheta + 1;
				double xTheta = Theta / angleStep - intTheta;
				if (intThetaNext > 6) // 90 deg is a maximum angle in sencetivity file!!!! change it if it is not!!!
				{
					intThetaNext = 6;
					xTheta = 0.;
				}
				double z = scalProd(V[v].fc[1].n, a);
				double x = scalProd(V[v].fc[3].n, a);
				double Fi = atan2(x, z);
				if (Fi < 0)
					Fi = twoPI + Fi;
				int intFi = Fi / angleStep;
				int intFiNext = intFi + 1;
				double xFi = Fi / angleStep - intFi;
				if (intFiNext > 23) // means full circle done use Fi=0 as next!!!
				{
					intFiNext = 0;
				}
				// interpolation511
				double aj = ListOfDetectors[V[v].detectorID].cs511[crystalInd].a[intTheta][intFi] +
							(ListOfDetectors[V[v].detectorID].cs511[crystalInd].a[intThetaNext][intFi] -
							 ListOfDetectors[V[v].detectorID].cs511[crystalInd].a[intTheta][intFi]) *
								xTheta;

				double aj1 = ListOfDetectors[V[v].detectorID].cs511[crystalInd].a[intTheta][intFiNext] +
							 (ListOfDetectors[V[v].detectorID].cs511[crystalInd].a[intThetaNext][intFiNext] -
							  ListOfDetectors[V[v].detectorID].cs511[crystalInd].a[intTheta][intFiNext]) *
								 xTheta;
				tmp.solidAngle511 = (aj + (aj1 - aj) * xFi) / norm2;

				// interpolation450
				aj = ListOfDetectors[V[v].detectorID].cs450[crystalInd].a[intTheta][intFi] +
					 (ListOfDetectors[V[v].detectorID].cs450[crystalInd].a[intThetaNext][intFi] -
					  ListOfDetectors[V[v].detectorID].cs450[crystalInd].a[intTheta][intFi]) *
						 xTheta;

				aj1 = ListOfDetectors[V[v].detectorID].cs450[crystalInd].a[intTheta][intFiNext] +
					  (ListOfDetectors[V[v].detectorID].cs450[crystalInd].a[intThetaNext][intFiNext] -
					   ListOfDetectors[V[v].detectorID].cs450[crystalInd].a[intTheta][intFiNext]) *
						  xTheta;
				tmp.solidAngle450 = (aj + (aj1 - aj) * xFi) / norm2;
			}
			// Ray from scatter point to module center - I think we do not need this
			// here will be simulated path integrals - ro*dx and activity*dx
			// keep in mind - the Ray is opposite to the ray used in this part of code!!!!
			// trace the phantom grid
			// start and end points---------------------------------------------
			vec p = ScatPoints[n].r;	// mkVec(-2.6,15.1,-60);		//ray origine scatter point
										// should be always inside voxel structure!!!!
			vec pd = r_i;				// crystal center//mkVec(-20,-220,-110);		//end point, can be any point!=start point
			vec aa = vecToScal(a, -1.); // invert existing Ray.a//
			// find curren voxel indexes and decriments for them during the tracking--------
			int ind[3] = {ScatPoints[n].i, ScatPoints[n].j, ScatPoints[n].k}; // i,j,k
			int incr[3] = {SGN(aa.x), SGN(aa.y), SGN(aa.z)};				  // incremet/decriment for indexes
			double t[3];													  // distance between planes along ray(should be >0)
			double tnext[3];												  // distance to the next voxel border intersections
																			  // from the origin point
			bool doNot[3] = {0, 0, 0};										  // do not track flag for parallel directions
			if (abs(aa.x) > TINYv)
			{
				t[0] = incr[0] * gStep.x / aa.x;
				if (incr[0] > 0)
					tnext[0] = (g0.x + gStep.x * (ind[0] + 1.) - p.x) / aa.x;
				else
					tnext[0] = (g0.x + gStep.x * ind[0] - p.x) / aa.x;
			}
			else
			{
				// t[0]=HUGEv;//not used
				doNot[0] = 1;
				tnext[0] = HUGEv;
			}
			if (abs(aa.y) > TINYv)
			{
				t[1] = incr[1] * gStep.y / aa.y;
				if (incr[1] > 0)
					tnext[1] = (g0.y + gStep.y * (ind[1] + 1.) - p.y) / aa.y;
				else
					tnext[1] = (g0.y + gStep.y * ind[1] - p.y) / aa.y;
			}
			else
			{
				// t[1]=HUGEv;//not used
				doNot[1] = 1;
				tnext[1] = HUGEv;
			}
			if (abs(aa.z) > TINYv)
			{
				t[2] = incr[2] * gStep.z / aa.z;
				if (incr[2] > 0)
					tnext[2] = (g0.z + gStep.z * (ind[2] + 1.) - p.z) / aa.z;
				else
					tnext[2] = (g0.z + gStep.z * ind[2] - p.z) / aa.z;
			}
			else
			{
				// t[2]=HUGEv;//not used
				doNot[2] = 1;
				tnext[2] = HUGEv;
			}
			bool in = true;		 // flag that you stil inside the grid and between given point
			double range;		 // range inside current voxel
			double tcurrent = 0; // current position of tracking
			int in1, in2, in3;
			do
			{
				in1 = ind[0];
				in2 = ind[1];
				in3 = ind[2];
				// find min next
				int min = 0;
				if (tnext[0] > tnext[1])
					min = 1;
				if (tnext[min] > tnext[2])
					min = 2;
				if (tnext[min] > norm)
				{
					range = norm - tcurrent;
					tcurrent = norm;
					in = false;
					// Do What you need here!!!
					tmp.ro_L += Dmumap(ind[0], ind[1], ind[2]) * range;
					tmp.act_L += Dimg(ind[0], ind[1], ind[2]) * range;
				}
				else
				{
					range = tnext[min] - tcurrent;
					tcurrent = tnext[min];
					// Do What you need here!!!
					tmp.ro_L += Dmumap(ind[0], ind[1], ind[2]) * range;
					tmp.act_L += Dimg(ind[0], ind[1], ind[2]) * range;
					// increment index//it should be never incremented for parallel beams!!!
					ind[min] += incr[min];
					if (ind[min] >= max[min] || ind[min] < 0)
						in = false; // go outside voxel structure - stop tracking
					else
						tnext[min] += t[min];
				}
			} while (in);
			// PathesMatrix[n*V.size()+i]=tmp;
			PathesMatrix[n * TotalCrystals + NN] = tmp;
			NN++;
			NNN++;
		}
	}
	std::cout << NNN << " pathes" << endl;
	double pSize = NNN / 1024. / 1024. / 1024. * sizeof(scaterPath);
	std::cout << pSize << " Gb pases matrix size" << endl;

	return true;
}

// Compton cotribution
double Single_Scatter_Simulation::ComptonCorr(int crystal1, int crystal2, double TOF_dist)
{
	double value = 0;

	double cosLimit = 2. - 511. / Ethreshold; // if cos of scatter angle is less there is no detection!!!
	double cosHighLimit = 0.999;			  // here is the restriction ONLY LARGE SCATTER ANGLES CONSIDERED!!!!!!
	cryst First = CrystPosFromID(crystal1);
	cryst Second = CrystPosFromID(crystal2);

	// loop over scatter points/centers
	for (int s = 0; s < ScatPoints.size(); s++)
	{
		int ind1 = s * TotalCrystals + crystal1;
		int ind2 = s * TotalCrystals + crystal2;

		double act1_L = 0;
		double act2_L = 0;

		// cos of scatter angle is -scalar product of direction to the crystals!!!
		vec a1 = mkVec(PathesMatrix[ind1].ax, PathesMatrix[ind1].ay, PathesMatrix[ind1].az);
		vec a2 = mkVec(PathesMatrix[ind2].ax, PathesMatrix[ind2].ay, PathesMatrix[ind2].az);

		double cosOfScatter = -scalProd(a1, a2); //"-" because one vector should be in opposite direction

		double Escatter = 511. / (2. - cosOfScatter);
		// ONLY LARGE SCATTER ANGLES CONSIDERED!!!!!!
		if (cosOfScatter > cosLimit && cosOfScatter < cosHighLimit)
		{
			// find center of trajectory
			double LdT = (PathesMatrix[ind1].L + PathesMatrix[ind2].L) / 2. + TOF_dist;

			vec rs = ScatPoints[s].r;
			if (LdT <= PathesMatrix[ind1].L)
			{
				double d = PathesMatrix[ind1].L - LdT;
				double dX = (rs.x - d * a1.x) - Dorigine[0];
				double dY = (rs.y - d * a1.y) - Dorigine[1];
				double dZ = (rs.z - d * a1.z) - Dorigine[2];
				int i = dX / Dstep;
				int j = dY / Dstep;
				int k = dZ / Dstep;
				if (i >= 0 && i < Ddimx && j >= 0 && j < Ddimy && k >= 0 && k < Ddimz && Dimg(i, j, k) > 0)
					act1_L = Dimg(i, j, k) * Dstep;
			}
			else
			{
				double d = LdT - PathesMatrix[ind1].L;
				double dX = (rs.x - d * a2.x) - Dorigine[0];
				double dY = (rs.y - d * a2.y) - Dorigine[1];
				double dZ = (rs.z - d * a2.z) - Dorigine[2];
				int i = dX / Dstep;
				int j = dY / Dstep;
				int k = dZ / Dstep;
				if (i >= 0 && i < Ddimx && j >= 0 && j < Ddimy && k >= 0 && k < Ddimz && Dimg(i, j, k) > 0)
					act2_L = Dimg(i, j, k) * Dstep;
			}

			double Escatter = 511. / (2. - cosOfScatter);
			// compton cross section - without constant coefficiens - the result will be scaled anyway!!!!
			double Pcoeff = 1. / (2 - cosOfScatter);
			double Pcoeff2 = Pcoeff * Pcoeff;
			double dSigma_dOmega = Pcoeff - Pcoeff2 * (1 - cosOfScatter * cosOfScatter) + Pcoeff2 * Pcoeff;
			// double mwaterEscatter=(-9.1E-6)*Escatter+0.0142;//this is in muxro (1/mm) units: if it will be in mm2/g unit x1000
			double mwaterEscatter = (-9.1E-3) * Escatter + 14.2;

			// int crystalInd1=nz1*ListOfDetectors[V[indV1].detectorID].repeatX+nx1;
			// int crystalInd2=nz2*ListOfDetectors[V[indV2].detectorID].repeatX+nx2;

			double solidAngle511_1 = PathesMatrix[ind1].solidAngle511;
			double solidAngle511_2 = PathesMatrix[ind2].solidAngle511;

			double solidAngle450_1 = PathesMatrix[ind1].solidAngle450;
			double solidAngle450_2 = PathesMatrix[ind2].solidAngle450;

			double eff_first = solidAngle511_1;
			double eff_second = solidAngle450_2 + (solidAngle511_2 - solidAngle450_2) * (Escatter - Ethreshold) / (511. - Ethreshold);
			// double attenuat	=exp(- mwater511*PathesMatrix[ind1].ro_L/0.001- mwaterEscatter*PathesMatrix[ind2].ro_L/0.001 );
			// double activity = PathesMatrix[ind1].act_L;//*emmiter dencity in the object!!!!
			double attenuat = exp(-mwater511 * PathesMatrix[ind1].ro_L - mwaterEscatter * PathesMatrix[ind2].ro_L);
			// recalculate according TOF correction
			// double Ifirst=eff_first*eff_second*attenuat*PathesMatrix[ind1].act_L;//IactAS;//activity;
			double Ifirst = eff_first * eff_second * attenuat * act1_L; // IactAS;//activity;
			// second path
			eff_first = solidAngle450_1 + (solidAngle511_1 - solidAngle450_1) * (Escatter - Ethreshold) / (511. - Ethreshold);
			eff_second = solidAngle511_2;
			// eff_second=PathesMatrix[secondInd].SolidAngle511;
			// attenuat=exp(-mwater511*PathesMatrix[ind2].ro_L/0.001- mwaterEscatter*PathesMatrix[ind1].ro_L/0.001 );
			attenuat = exp(-mwater511 * PathesMatrix[ind2].ro_L - mwaterEscatter * PathesMatrix[ind1].ro_L);
			// activity = PathesMatrix[ind2].act_L;
			// double Isecond=eff_first*eff_second*attenuat*PathesMatrix[ind2].act_L;
			double Isecond = eff_first * eff_second * attenuat * act2_L;
			// value+=(Ifirst+Isecond)*dSigma_dOmega*Dmumap[ScatPoints[s].i][ScatPoints[s].j][ScatPoints[s].k]/MuToRo;//should be multiplied on mu!!!!!!
			value += (Ifirst + Isecond) * dSigma_dOmega * Dmumap(ScatPoints[s].i, ScatPoints[s].j, ScatPoints[s].k);
		}
	}
	// value=value*ScaleCompt;
	return value;
}

// random generators
float Single_Scatter_Simulation::ran1(long *idum)
//"Minimal" random number generator of Park and Miller with Bays-Durham shuffle and added safeguards
// returns (0,1) !!!NOT [0,1]!!!
// Call with idum a negative integer to initialize;
// thereafter, do not alter idum between successive deviates in a sequence.
// RNMX should approximate the largest floating value that is less than 1
//(from "Numerical recipes in C" p.280)
{
	int j;
	long k;
	static long iy = 0;
	static long iv[NTAB];
	float temp;
	if (*idum <= 0 || !iy) // initialaize
	{
		if (-(*idum) < 1)
			*idum = 1; // to be shure that idum!=0
		else
			*idum = -(*idum);
		for (j = NTAB + 7; j >= 0; j--) // load the shaffle table(after 8 warm-ups
		{
			k = (*idum) / IQ;
			*idum = IA * (*idum - k * IQ) - IR * k;
			if (*idum < 0)
				*idum += IM;
			if (j < NTAB)
				iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ;						// start here when not initializing
	*idum = IA * (*idum - k * IQ) - IR * k; // compute idum=(IA*idum)%IM without overflows by Schrage's method
	if (*idum < 0)
		*idum += IM;
	j = iy / NDIV; // will be in the range 0..NTAB-1
	iy = iv[j];	   // output previosly stored value and refill the shuffle table
	iv[j] = *idum;
	if ((temp = AM * iy) > RNMX)
		return RNMX; // becouse no endpoint values!!!
	else
		return temp;
}
float Single_Scatter_Simulation::ran2(long *idum)
// Long period random number generator of L'Ecuyer with Bays-Durham shuffle and added safeguards
// returns (0,1) !!!NOT [0,1]!!!
// Call with idum a negative integer to initialize;
// thereafter, do not alter idum between successive deviates in a sequence.
// RNMX should approximate the largest floating value that is less than 1
//(from "Numerical recipes in C" p.282)
{
	int j;
	long k;
	static long idum2 = 123456789;
	static long iy = 0;
	static long iv[NTAB];
	float temp;
	if (*idum <= 0) // initialaize
	{
		if (-(*idum) < 1)
			*idum = 1; // to be shure that idum!=0
		else
			*idum = -(*idum);
		idum2 = (*idum);
		for (j = NTAB + 7; j >= 0; j--) // load the shaffle table(after 8 warm-ups
		{
			k = (*idum) / IQ1;
			*idum = IA1 * (*idum - k * IQ1) - k * IR1;
			if (*idum < 0)
				*idum += IM1;
			if (j < NTAB)
				iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ1;						   // start here when not initializing
	*idum = IA1 * (*idum - k * IQ1) - k * IR1; // compute idum=(IA*idum)%IM without overflows by Schrage's method
	if (*idum < 0)
		*idum += IM1;
	k = idum2 / IQ2;
	idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
	if (idum2 < 0)
		idum2 += IM2;
	j = iy / NDIV;		// will be in the range 0..NTAB-1
	iy = iv[j] - idum2; // here idum is shuffled, idum andidum2 are combined to generate output
	iv[j] = *idum;
	if (iy < 1)
		iy += IMM1;
	if ((temp = AM_2 * iy) > RNMX)
		return RNMX; // becouse no endpoint values!!!
	else
		return temp;
}
float Single_Scatter_Simulation::gasdev(long *idum) // From Numerical recipes
{
	// float ran1(long *idum);
	static int iset = 0;
	static float gset;
	float fac, rsq, v1, v2;
	if (iset == 0)
	{
		do
		{
			v1 = 2.0 * ran1(idum) - 1.0;
			v2 = 2.0 * ran1(idum) - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= 1.0 || rsq == 0);
		fac = sqrt(-2.0 * log(rsq) / rsq);
		gset = v1 * fac;
		iset = 1;
		return v2 * fac;
	}
	else
	{
		iset = 0;
		return gset;
	}
}

int Single_Scatter_Simulation::CrystalInd(cryst IN) // Crystall ID from .detector, .sector, ...
{
	int crystalID =
		IN.sector * ListOfDetectors[IN.detector].Nsector +
		ListOfDetectors[IN.detector].Nmodule * ListOfDetectors[IN.detector].modulesX * IN.moduleZ +
		ListOfDetectors[IN.detector].Nmodule * IN.moduleX +
		IN.nx +
		IN.mz * ListOfDetectors[IN.detector].repeatX;
	for (int d = 0; d < IN.detector; d++)
		crystalID += ListOfDetectors[d].N;
	return crystalID;
}
cryst Single_Scatter_Simulation::CrystPosFromID(int crystalID)
{
	cryst out;
	// find detector
	int detector = 0;
	int volume = 0;
	for (int i = 0; i < ListOfDetectors.size(); i++)
	{
		if (crystalID < ListOfDetectors[i].N)
			break;
		crystalID -= ListOfDetectors[i].N;
		detector++;
		volume += ListOfDetectors[i].volumes;
	}
	volume += crystalID / ListOfDetectors[detector].Nmodule;
	out.V = volume;
	// find sector
	out.detector = detector;
	int sector = crystalID / ListOfDetectors[detector].Nsector;
	out.sector = sector;

	int IDinSector = crystalID % ListOfDetectors[detector].Nsector;
	int inRow = ListOfDetectors[detector].modulesX * ListOfDetectors[detector].Nmodule;

	int moduleZ = IDinSector / inRow;
	out.moduleZ = moduleZ;

	int IDinRow = IDinSector % inRow;
	int moduleX = IDinRow / ListOfDetectors[detector].Nmodule;
	out.moduleX = moduleX;

	int IDinModule = IDinRow % ListOfDetectors[detector].Nmodule;

	out.mz = IDinModule / ListOfDetectors[detector].repeatX;
	out.nx = IDinModule % ListOfDetectors[detector].repeatX;
	return out;
}
