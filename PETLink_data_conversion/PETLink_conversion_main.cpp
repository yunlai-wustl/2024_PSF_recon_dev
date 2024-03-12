/*
Decode siemens format and write out out own defined data format
change:
1. 	vector<*> coinc_event_list; definition of output data type
2. record_time_alignemnt_coinc, call right func;
3. output.write(reinterpret_cast<char *>(&coinc_event_list[0]), coinc_event_list.size()*sizeof(time_alignment_event));
4. std::vector<time_alignment_event>().swap(coinc_event_list); //free memory
*/




#include "../Solution_Items/GATE_data_structure.h"
#include <stdio.h>
#include <stdlib.h>
using namespace std;
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#define MAX_COINC 50000000
void record_mini_coinc(vector<mini_coinc_event>& mini_coinc_event_ptr, int event_count, int crystal1, int crystal2, double tof, float elapsed_time);
void record_cbm_coinc(vector<tof_coinc_event_cbm>& tof_coinc_event_cbm_list, int crystal1, int crystal2, float elapsed_time, double tof, float current_bed_position);
//global timing offset lookup table; 
//look up table; key: cyrstal index in GATE convention; value: timing offset value in ps; 
//correction = uncorrected + time_offset[crystal_2] - time_offset[crystal_1]
std::map<int, double> LUT_Timing_offset; 

int main(int argc, char * argv[])
{

	int crystal1, crystal2, block1, block2;
	int nRandom = 0;
	double tof;
	int event_count = 0;


	//char file_name[155];
	int  pt, tf;
	int x[2];
	// note that pair (XE) strats from 1
	unsigned char XE_lookup_A[155] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, \
		2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, \
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, \
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, \
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, \
		6, 6, 6, 6, 6, 6, 6, 6, 6, 6, \
		7, 7, 7, 7, 7, 7, 7, 7, 7, \
		8, 8, 8, 8, 8, 8, 8, 8, \
		9, 9, 9, 9, 9, 9, 9, \
		10, 10, 10, 10, 10, 10, \
		11, 11, 11, 11, 11, \
		12, 12, 12, 12, \
		13, 13, 13, \
		14, 14, \
		15, \
		0, 1, 2, 3, 4, \
		14, 15, 16, 17, 18, \
		4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, \
		19
	};
	unsigned char XE_lookup_B[155] = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, \
		4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, \
		5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, \
		6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, \
		7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, \
		8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, \
		9, 10, 11, 12, 13, 14, 15, 16, 17, 18, \
		10, 11, 12, 13, 14, 15, 16, 17, 18, \
		11, 12, 13, 14, 15, 16, 17, 18, \
		12, 13, 14, 15, 16, 17, 18, \
		13, 14, 15, 16, 17, 18, \
		14, 15, 16, 17, 18, \
		15, 16, 17, 18, \
		16, 17, 18, \
		17, 18, \
		18, \
		19, 19, 19, 19, 19, 19, 19, 19, 19, 19, \
		20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20 };

	int dea1, dea2, ax, bx, ay, by, xe;
	int tf1, tf2, tf3, tf4, tf5;
	int max_block_1=-1;
	int max_block_2=-1;
	unsigned int Xsize;
	unsigned int total_num_counts=0;
	int time_tag_num = 0;
	int cnt = 0;
	int sign = 0;
	int file_ind = 0;
	int sync, event_tag;	
	float time_start= 1000000.0f;//start time
	float time_end=-1.0f;
	float time_tag = 0.0f;//time of event 1
	float current_bed_position = 0.0; //bed position, in mm
	float min_bed_position=10000.0f;
	float max_bed_position=-10000.0f;
	int tag_word;
	FILE *fp;
	int method;
	//FILE *output;
	//size_t t;
	for (int i = 0; i < 60800 + 900 * 64; i++)
		LUT_Timing_offset[i] = 0.0f; //initialize LUT
	if (argc == 3){ 
		printf("using timing offset lookuptable; usage: %s inputdatafilename outputdatafilename lookuptable \n", argv[0]);
		ifstream if_T(argv[2]);
		int crystal_ind;
		float timing_offset;
		while (if_T >> crystal_ind >> timing_offset){
			LUT_Timing_offset[crystal_ind] = timing_offset;
			//cout << "crystal id:" << crystal_ind << ", time offset:" << LUT_Timing_offset[crystal_ind] << endl;;
		}
		cout << "timing offset look up table read successfully" << endl;
	} // Usage
	if (argc == 2){ 
		printf("no timing offset lookuptable; usage: %s inputdatafilename \n", argv[0]); 
		
	} // Usage

	if (NULL == (fp = fopen(argv[1], "rb"))){ exit(1); } // IF File name is not proper
	vector<tof_coinc_event_cbm> coinc_event_list;
	//vector<mini_coinc_event> coinc_event_list;
	//vector<time_alignment_event> coinc_event_list;
	string outputfileName = "converted_" + string(argv[1]);
	vector<vector<float>> random_matrix_block(336, vector<float>(336, 0.0f));
	coinc_event_list.reserve(MAX_COINC);
	ofstream output(outputfileName, ios_base::binary);
	while (1){
		size_t rc = fread(x, 4, 2, fp); // read "4 bytes (8bit * 4) " 2 time
		if (rc != 2)
		{
			printf("end of file \n");
			break;
		}
		//--------------------------------extract bit 30,31 of both (sync, prompt/Tag_56PL (pt) , Tag_64(event tag) )------//
		sync = ((1 - x[0] & 0x80000000) >> 31)&((x[1] & 0x80000000) >> 31);
		event_tag = (x[0] & 0x40000000) >> 30; //event tage, 1 indicates non-event, 0 indicates prompt/delayed events
		pt = (x[1] & 0x40000000) >> 30; //when  event_tag=1, pt=1 indicates 32-bit payload, 0 for 56-bit payload; 
		//when event_tag=0,  pt=1 is prompt event, pt=0 is delayed events


		//syc!=1: not synchronized
		//EVENT_TAG=1: non-event
		if (sync != 1) //not synchronized
		{
			continue;
		}
		//---------------------------------------------------process non event---------------------------------------------//

		if (event_tag == 1) //non-event, see if this is 23/56 bit payload
		{	//check if it is 32-bit tag packet payload (pt=0: 32-bit; pt=1: 56 bit)

			if (pt == 1) //56-bit packet payload, singles rate
				continue;

			//otherwise this is a 32bit payload, extract tag_word;
			tag_word = (x[0] & 0x0000ffff) + ((x[1] & 0x0000ffff) << 16);
			//time, see elpased time or dead time tracker
			int tag = (tag_word & 0xf0000000) >> 28;
			switch (tag){
			case 8:{ //elapsed time marker
				time_tag = float(tag_word & 0x1fffffff) / 1000.0f; //in s
				time_start = min(time_tag,time_start);
				time_end = max(time_end,time_tag);
				continue;
			}
			case 12:{ //gantry motions & position
				//bed position, see if this is horizontal bed position
				if (((tag_word & 0x0f000000) >> 24) == 4){
					current_bed_position = ((tag_word & 0x000fffff) - (1 << 20)) / 100.0f; //in mm
					min_bed_position=min(min_bed_position,current_bed_position);
					max_bed_position=max(max_bed_position,current_bed_position);
				}
				continue;
			}
			case 14: //patient monitoring
				continue;
			case 15: //control/acquisition parameters
				continue;
			}
			continue;
		}

		//else: it is an event; extract crystal/block id first, then check if it is delayed event (pt==0) or prompt event (pt=1)
		//-----------------------------------------------------decode event------------------------------------------------------------------//
		xe = ((x[1] & 0x00070000) >> 13) + ((x[0] & 0x00070000) >> 16) + ((x[0] & 0x00000040)) + ((x[1] & 0x00000040) << 1);  // module pair
		//xe<=133: Scaner-Scanner coincidence;
		if (xe <= 133){
			dea1 = XE_lookup_A[xe - 1];
			dea2 = XE_lookup_B[xe - 1];
			ax = x[0] & 0x0000003f;
			ay = (x[0] & 0x00007f00) >> 8;
			bx = x[1] & 0x0000003f;
			by = (x[1] & 0x00007f00) >> 8;

			if ((ax > 39) || (ay > 79) || (bx > 39) || (by > 79))
			{
				printf("error: tp:%d, dea1:%d, ax:%d, ay:%d, dea2:%d, bx:%d, by:%d\n", pt, dea1, ax, ay, dea2, bx, by);
				continue;
			}
			/*
			crystal1 = (23 - dea1) % 19 * 3200 + (1 - ax / 20) * 1600 + ay * 20 + (39 - ax) % 20;
			crystal2 = (23 - dea2) % 19 * 3200 + (1 - bx / 20) * 1600 + by * 20 + (39 - bx) % 20;
			*/
			crystal1 = (dea1-14+19) % 19 * 3200 + (ax / 20) * 1600 + (79 - ay) * 20 + ax % 20;
			crystal2 = (dea2-14+19) % 19 * 3200 + (bx / 20) * 1600 + (79 - by) * 20 + bx % 20;


			if (crystal1 < 0 || crystal1>60799)
			{
				printf("error: tp:%d,crystal1:%d, dea1:%d, ax:%d, ay:%d\n", pt, crystal1, dea1, ax, ay);
				continue;
			}if (crystal2 < 0 || crystal2>60799)
			{
				printf("error: crystal2:%d, dea2:%d, bx:%d, by:%d\n", crystal2, dea2, bx, by);
				continue;
			}
		}
		//134=<xe<=154: Outsert-Scanner coincdence
		else if (xe <= 154)
		{
			dea1 = XE_lookup_A[xe - 1];
			dea2 = XE_lookup_B[xe - 1];
			ax = x[0] & 0x0000003f;
			ay = (x[0] & 0x00007f00) >> 8;
			bx = (x[1] & 0x0000003f) + ((x[1] & 0x00600000) >> 15);
			by = ((x[1] & 0x00007f00) >> 8) + ((x[1] & 0x01800000) >> 16);
			//crystal1 = (23 - dea1) % 19 * 3200 + (1 - ax / 20) * 1600 + ay * 20 + (39 - ax) % 20;
			crystal1 = (5 + dea1) % 19 * 3200 + (ax / 20) * 1600 + (79 - ay) * 20 + ax % 20;
			crystal2 = 60800 + 120*240*(dea2 - 19) + bx/30*7200+by/30*900+bx%30*30+by%30; //bx:0-119, by:0-239
			

		}
		//xe==155: Outsert-Outsert coincdence
		else if (xe == 155)
		{
			dea1 = XE_lookup_A[xe - 1];
			dea2 = XE_lookup_B[xe - 1];
			ax = (x[0] & 0x0000003f) + ((x[0] & 0x00600000) >> 15);
			ay = ((x[0] & 0x00007f00) >> 8) + ((x[0] & 0x01800000) >> 16);
			bx = (x[1] & 0x0000003f) + ((x[1] & 0x00600000) >> 15);
			by = ((x[1] & 0x00007f00) >> 8) + ((x[1] & 0x01800000) >> 16);
			crystal1 = 60800 + 120 * 240 * (dea1 - 19) + ax / 30 * 7200 + ay / 30 * 900 + ax % 30 * 30 + ay % 30;
			crystal2 = 60800 + 120 * 240 * (dea2 - 19) + bx / 30 * 7200 + by / 30 * 900 + bx % 30 * 30 + by % 30;

			//crystal1 = 60800 + by + (127 - bx) * 255;
			//crystal2 = 60800 + 128 * 256 + by + (127 - bx) * 255;
		}
		else
		{
			continue;
			printf("error: xe=%d\n", xe);
		}

		if (pt == 0) //delayed coincdence, record this to random matrix and move to next; get the rate when write out
		{	
			//SS
			if(xe <= 133){
			block1 = (5+dea1) % 19 * 16 + (1 - ax / 20) * 8 + ay / 10; //Gate definition
			block2 = (5+dea2) % 19 * 16 + (1 - bx / 20) * 8 + by / 10;
			}
			/*
			block1 = (23 - dea1) % 19 * 16 + (1 - ax / 20) * 8 + ay / 10; //Gate definition
			block2 = (23 - dea2) % 19 * 16 + (1 - bx / 20) * 8 + by / 10;
			*/
			//block2 = dea2 * 16 + by / 10 * 2 + bx / 20;
			else if (xe <= 154){
				block1 = (5+dea1) % 19 * 16 + (1 - ax / 20) * 8 + ay / 10; //Gate definition
				block2 = (crystal2-60800)/1800+304;
			}
			else if (xe == 155){	
				block1= (crystal1-60800)/1800+304;	
				block2 = (crystal2-60800)/1800+304;

			}
			max_block_1=max(max_block_1,block1);
			max_block_2=max(max_block_2,block2);

			random_matrix_block[block1][block2] += 1.0f;
			random_matrix_block[block2][block1] += 1.0f;
			nRandom++;
			continue;
		}
		//prompt event, continue decoding
		tf1 = (x[0] & 0x0E000000) >> 25;
		tf2 = (x[1] & 0x0E000000) >> 22;  // 3-25
		tf3 = (x[0] & 0x10000000) >> 22; // 6-28
		tf4 = (x[1] & 0x10000000) >> 21; // 7-28
		tf5 = (x[0] & 0x20000000) >> 21; //8-29
		tf = tf1 + tf2 + tf3 + tf4 + tf5;
		if (tf <= 255)
		{
			tof = 13.021 * tf;
		}
		else if (tf == 256)
		{
			tof = 0;
		}
		else
		{
			tof = 13.021 * (tf - 512); //in ps
		}
		//---------------------------------write out prompt when it reaches the max_coinc---------------------------------------------------//
		Xsize = coinc_event_list.size();
		if (int(current_bed_position)!=0 && Xsize < MAX_COINC){
			record_cbm_coinc(coinc_event_list, crystal1, crystal2, time_tag, tof, current_bed_position);
			//record_mini_coinc(coinc_event_list, event_count, crystal1, crystal2, (double)tof, (double)time_tag);
			//record_time_alignemnt_coinc(coinc_event_list, crystal1, crystal2, short int(tof / 13.021));
			event_count++;
		}
		//if the count reaches max count, write out buffer and reset the vector
		else
		{
			cout << "writing out " << Xsize << " coincidence events" << endl;
			output.write(reinterpret_cast<char *>(&coinc_event_list[0]), coinc_event_list.size()*sizeof(tof_coinc_event_cbm));
			//output.write(reinterpret_cast<char *>(&coinc_event_list[0]), coinc_event_list.size()*sizeof(mini_coinc_event));
			//output.write(reinterpret_cast<char *>(&coinc_event_list[0]), coinc_event_list.size()*sizeof(time_alignment_event));
			total_num_counts += MAX_COINC;
			coinc_event_list.clear();
			event_count = 0;
		}
	}

	//------------------------------------write out  files ------------------------------------------------------------//

	Xsize = coinc_event_list.size();
	total_num_counts += Xsize;
	cout << "writing out " << Xsize << " coincidence events" << endl;
	cout <<"total number of coincidence:"<<total_num_counts<< ", filesize should be:"<< total_num_counts*sizeof(tof_coinc_event_cbm) / 1024 / 1024 << "MB" << endl;
	cout << "number of delayed coincidence:" << nRandom << ",random rate:" << nRandom*100 / total_num_counts <<"%"<< endl;
	cout<<"max block 1:" <<max_block_1<<", max block 2:" <<max_block_2<<endl;
	output.write(reinterpret_cast<char *>(&coinc_event_list[0]), coinc_event_list.size()*sizeof(tof_coinc_event_cbm));
	//output.write(reinterpret_cast<char *>(&coinc_event_list[0]), coinc_event_list.size()*sizeof(mini_coinc_event));
	//output.write(reinterpret_cast<char *>(&coinc_event_list[0]), coinc_event_list.size()*sizeof(time_alignment_event));

	output.close();
	std::vector<tof_coinc_event_cbm>().swap(coinc_event_list); //free memory
	//std::vector<mini_coinc_event>().swap(coinc_event_list); //free memory
	//std::vector<time_alignment_event>().swap(coinc_event_list); //free memory

	cout << "scan start time " << time_start <<"s, end time:"<<time_end<<"s, scan duration:"<< time_end-time_start <<"s"<<endl;
	cout << "min bed position: "<<min_bed_position<<"mm, max bed position: "<< max_bed_position<< "mm" <<endl;
	// convert count to count rate and wirte out//	
	std::ofstream file_random_matrix;
	string filename_random_matrix = "random_"+string(argv[1])+".txt";

	file_random_matrix.open(filename_random_matrix);
	for (int i = 0; i < random_matrix_block.size(); i++){
		for (int j = 0; j < random_matrix_block[i].size(); j++)
		{
			file_random_matrix << random_matrix_block[i][j] << "	";
		}
		file_random_matrix << std::endl;
	}
	file_random_matrix.close();


	return 0;
}

void record_mini_coinc(vector<mini_coinc_event>& mini_coinc_event_ptr, int event_count, int crystal1, int crystal2, double tof, float elapsed_time)
{
	mini_coinc_event new_event;
	new_event.crystal_index_1 = crystal1;
	new_event.crystal_index_2 = crystal2;
	new_event.time_1 = elapsed_time;
	new_event.diff_time = tof*(1e-12) + (LUT_Timing_offset[crystal2] - LUT_Timing_offset[crystal1])*(1e-12); //in s 
	mini_coinc_event_ptr.push_back(new_event);
}


void record_cbm_coinc(vector<tof_coinc_event_cbm> &tof_coinc_event_cbm_list, int crystal1, int crystal2, float elapsed_time, double tof, float current_bed_position)
{
	tof_coinc_event_cbm new_event;
	new_event.crystal_index_1 = crystal1;
	new_event.crystal_index_2 = crystal2;
	new_event.time_1 = elapsed_time;
	new_event.diff_time = tof*(1e-12) + (LUT_Timing_offset[crystal2] - LUT_Timing_offset[crystal1])*(1e-12); //in s 
	new_event.bed_position = current_bed_position;
	tof_coinc_event_cbm_list.push_back(new_event);
}
