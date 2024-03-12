// The coordinate of this sorting code is from the robot side
// the direction to the right hand side is positive x, to up is positive Y, to yourself is positive Z
// for the Scanner itself, ax goes counterclock direction, ay goes to the negtive Z
// by Jianyong Jiang

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>

#include "../Solution_Items/global.h"
#include "../Solution_Items/GATE_data_structure.h"
#include "../Solution_Items/command_line.h"
#include "../Solution_Items/config.h"
#include "../Solution_Items/time_period.h"

#define MAX_COINC 2000000000
#define MAX_WRITE_BYTES 1000000000

void Usage(char argv0[])
{
	string program_name = CommandLine::GetProgramName(argv0);
	printf("Usage: %s <data_list_file>\n", program_name.c_str());
}


int sort_biograph_with_insert_non_tof(string inputfilename, vector<nontof_coinc_event> &list_data, float time_start, float time_end, float add_time_tag_offset)
{
	int tp, tg, pt;
	unsigned char x[8];
	unsigned char XE_lookup_SS_II_A[42] = { 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7, 7, 8 };
	unsigned char XE_lookup_SS_II_B[42] = { 3, 4, 5, 6, 7, 8, 9, 4, 5, 6, 7, 8, 9, 10, 5, 6, 7, 8, 9, 10, 11, 6, 7, 8, 9, 10, 11, 7, 8, 9, 10, 11, 8, 9, 10, 11, 9, 10, 11, 10, 11, 11 };
	unsigned char XE_lookup_IS_A[18] = { 0, 1, 2, 3, 4, 8, 9, 10, 11, 0, 1, 2, 3, 4, 8, 9, 10, 11 };
	unsigned char XE_lookup_IS_B[18] = { 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13 };
	int SS_SIZE = 42;
	int IS_SIZE = 18;
	int MAX_SIZE = 64;
	int Num_cry_per_Block_Scanner = 169;
	int Num_cry_per_array_Scanner = 13;
	int Num_cry_per_Block_Insert = 256;
	int Num_cry_per_array_Insert = 16;
	int Num_cry;
	int dea1, dea2, ax, bx, ay, by;
	int blk1, blk2, scr1, scr2, cryID1, cryID2;
	float time;
	int time_tag_num = 0;
	int MP;
	int EMP;
	int cry_index1, cry_index2;
	int eventID = 0;
	int eventRecord = 0;
	FILE *fid;

	int index_start;
	int index_end;


	fid = fopen(inputfilename.c_str(), "rb");
	if (fid == NULL){
		std::cout << inputfilename << " can NOT be openned. Continued to process next file." << endl;
		return 0;
	}

	std::cout << "File: " << inputfilename << endl;
	index_start = list_data.size();

	while (!feof(fid)){
		fread(x, 1, 8, fid); // read "1 byte * 8" one time
		tp = x[3] & 0x40;  // Tag_64 (0/1)
		tp = tp / 64;      // Tag == (0/1)?

		// According to the Tag_64, we perform different data decoding

		if (tp == 1){ // TYPE 1: non-event Tag
			//	printf("Time marker");
			tg = x[5] & 0xe0; //Tag word, 100x xxxx Time Marker
			if (tg == 128){
				time_tag_num++;
				time = (float)((x[5] & 0x1f) * 16777216 + x[4] * 65536 + x[1] * 256 + x[0]) / 1000.0;
			}

		}

		else if (tp == 0){ // TYPE 0: EVENT DATA DECODING 
			pt = (x[7] & 0x40) / 64; //Prompt tag

			MP = ((x[6] & 0x07) * 8) + (x[2] & 0x07);  // module pair
			EMP = (x[6] & 0x38) + ((x[2] & 0x38) >> 3);  // module pair


			// determine if the ring MP is non-zero
			if ((MP != 0) && (MP <= SS_SIZE))
			{
				dea1 = XE_lookup_SS_II_A[MP - 1];   // Event S
				dea2 = XE_lookup_SS_II_B[MP - 1];   // Event S
				Num_cry = Num_cry_per_array_Scanner;
				ax = x[0];
				ay = x[1];
				bx = x[4];
				by = x[5];

				//switch ay from negtive Z to positive Z

				blk1 = 3 - ay / Num_cry;
				blk2 = 3 - by / Num_cry;

				scr1 = (dea1 * 4 + ax / Num_cry + 10) % 48;
				scr2 = (dea2 * 4 + bx / Num_cry + 10) % 48;
				// rotate the scanner in order to let the starting crystal ID same with the gate geometry

				cryID1 = (12 - ay % 13) * 13 + ax % 13;
				cryID2 = (12 - by % 13) * 13 + bx % 13;

				cry_index1 = 169 * (scr1 * 4 + blk1) + cryID1;
				cry_index2 = 169 * (scr2 * 4 + blk2) + cryID2;
			}

			// determine if the insert module pair is non-zero
			else if ((MP == 0) && (EMP != 0) && (EMP <= IS_SIZE))
			{
				dea1 = XE_lookup_IS_A[EMP - 1];    //Event S
				dea2 = XE_lookup_IS_B[EMP - 1];    //Event I
				Num_cry = Num_cry_per_array_Insert;
				bx = x[0];
				by = x[1];
				ax = x[4];
				ay = x[5];


				blk1 = 3 - ay / 13;
				scr1 = (dea1 * 4 + ax / 13 + 10) % 48;
				// rotate the scanner in order to let the starting crystal ID same with the gate geometry
				cryID1 = (12 - ay % 13) * 13 + ax % 13;
				cry_index1 = 169 * (scr1 * 4 + blk1) + cryID1;             // find the cry_index1 for event in the scanner

				//if ((bx == 25 || bx == 26) && (by == 25 || by == 26)){
				//	if ((datout = fopen("IS_datout.dat", "at+")) == NULL){ printf("Cannot open file"); exit(1); }
				//	fprintf(datout, "%d %d %d %d %d %d\n", dea2, ax, ay, dea1, bx, by);
				//	fclose(datout);
				//}
				//if (ay > 51){ printf("%d %d %d %d %d %d %d %d %d\n", dea1, ax, ay, blk1, scr1, cryID1, dea2, bx, by); }


				if (dea2 == 13 && (bx > 15 && bx < 32))
				{
					if (by < 16) { bx = 47 - bx; by = 47 - by; }
					else if (by >15 && by < 32) { bx = 47 - bx; by = 79 - by; }
					else if (by > 31 && by < 48) { bx = 47 - bx; by = 47 - by; }
					else if (by > 47 && by < 64) { bx = 47 - bx; by = 79 - by; }
				}
				else if (dea2 == 13 && (bx > 47 && bx < 64))
				{
					if (by < 16) { bx = 111 - bx; by = 47 - by; }
					else if (by >15 && by < 32) { bx = 111 - bx; by = 79 - by; }
					else if (by > 31 && by < 48) { bx = 111 - bx; by = 47 - by; }
					else if (by > 47 && by < 64) { bx = 111 - bx; by = 79 - by; }
				}                          // correct the error connection in junction board for dea13


				blk2 = 3 - bx / Num_cry;
				scr2 = 3 - by / Num_cry;
				cryID2 = (15 - by % 16) * 16 + (15 - bx % 16);

				cry_index2 = 256 * (scr2 * 8 + blk2 + (13 - dea2) * 4) + cryID2 + 126336;
				if (cry_index2 >(126336 + 8192 - 1)){ printf("error!"); }

			}  // end of "else if ((MP == 0) && (EMP != 0) && (EMP <= IS_SIZE))"


			if (eventID < MAX_COINC && time<time_end && time>time_start){

				nontof_coinc_event coinc;
				coinc.crystal_index_1 = cry_index1;
				coinc.crystal_index_2 = cry_index2;
				coinc.time_1 = time + add_time_tag_offset;

				list_data.push_back(coinc);
				eventRecord++;

			}

			eventID++;
		} //end of else



	} // end of while
	fclose(fid);

	index_end = list_data.size() - 1;


	std::cout << eventID << " coincidence events with time span of " << time << " seconds." << endl;
	std::cout << eventRecord << " coincidence events from " << time_start << " to " << time_end << " are recorded." << endl;
	std::cout << "Event index in data structure for this file ranges from " << index_start << " to " << index_end << "." << endl;
	std::cout << "Time stamp in data structure for this file ranges from " << list_data[index_start].time_1 << " to " << list_data[index_end].time_1 << "." << endl;



	return 0;
}


int sort_compact_with_insert_non_tof(string inputfilename, vector<nontof_coinc_event> &list_data, float time_start, float time_end, float add_time_tag_offset){
	struct stat64  file_status;
	mini_coinc_event *compact_coincidence_data;
	int num_coincidence;
	FILE *inputfile;
	size_t t;
	__int64  size;


	printf("Size of mini_coinc_event structure is %d\n", (int)sizeof(mini_coinc_event));
	printf("Following coincidence files are being processed:\n");
	if (stat64(inputfilename.c_str(), &file_status) < 0) {
		printf("Read file error: %s \n", inputfilename.c_str());
		return 0;
	}

	num_coincidence = 0;

	size = file_status.st_size;
	printf("%I64d File size in byte \n", size);
	size /= 52;
	if (size >= INT_MAX){
		printf("File \"%s\" length is tooooooooooo long. This program only supports coincidence events fewer than %d. Current file has %lld. Program stopped\n", inputfilename.c_str(), INT_MAX, size);
		return 0;
	}


	/*
	if ((inputfile = fopen(inputfilename.c_str(), "rb")) == NULL) {
	printf("Error: Could not read input file. \n");
	return 0;
	}
	_fseeki64(inputfile, 0, SEEK_END);   // non-portable
	size = _ftelli64(inputfile);
	fclose(inputfile);
	printf("%I64d File size in byte \n", size);

	size /= 52;
	*/

	num_coincidence = (int)size;
	printf("%d coincidence events in file: %s\n", num_coincidence, inputfilename.c_str());




	//Coincidence data memory
	compact_coincidence_data = (mini_coinc_event*)malloc(num_coincidence*sizeof(mini_coinc_event));

	//Read in coincidence data
	printf("Opening list mode file, %s \n", inputfilename.c_str());
	if ((inputfile = fopen(inputfilename.c_str(), "rb")) == NULL) {
		printf("Error: Could not read input file. \n");
		return 0;
	}

	t = fread(compact_coincidence_data, sizeof(mini_coinc_event), num_coincidence, inputfile);
	if (num_coincidence == (int)t){
		printf("%s read successfully.\n", inputfilename.c_str());
	}
	else{
		printf("%s read error. %d elements were read\n", inputfilename.c_str(), (int)t);
		return 0;
	}
	fclose(inputfile);

	int eventID;
	float time_tag;
	for (eventID = 0; eventID < num_coincidence; eventID++){
		time_tag = (float)compact_coincidence_data[eventID].time_1;

		if (eventID < MAX_COINC && time_tag<time_end && time_tag>time_start){

			nontof_coinc_event coinc;
			coinc.crystal_index_1 = compact_coincidence_data[eventID].crystal_index_1;
			coinc.crystal_index_2 = compact_coincidence_data[eventID].crystal_index_2;
			coinc.time_1 = time_tag + add_time_tag_offset;

			list_data.push_back(coinc);
		}
	}
	free(compact_coincidence_data);

}

int encode_crystal_mCT(volume_id crystal_id){
	int crystal_index;
	int block_index;

	block_index = crystal_id.module + 4 * crystal_id.sector;
	crystal_index = 169 * block_index + crystal_id.crystal;

	return crystal_index;
}

int encode_crystal_mCT_Insert(volume_id crystal_id){
	int crystal_index;
	int block_index;

	if (crystal_id.system == 0){//belongs to Scanner
		block_index = crystal_id.module + 4 * crystal_id.sector;
		crystal_index = 169 * block_index + crystal_id.crystal;
	}
	else{//belongs to insert
		block_index = crystal_id.module;
		crystal_index = 256 * block_index + crystal_id.crystal + 126336;
	}

	return crystal_index;
}


int record_mini_coinc(GATE_singles_event *d_ptr1, GATE_singles_event *d_ptr2, mini_coinc_event *coinc_ptr){
	//record this coinc
	coinc_ptr->event_id = d_ptr1->eventID;
	coinc_ptr->crystal_index_1 = encode_crystal_mCT_Insert(d_ptr1->single_ID);
	coinc_ptr->crystal_index_2 = encode_crystal_mCT_Insert(d_ptr2->single_ID);
	coinc_ptr->time_1 = d_ptr1->time;
	coinc_ptr->diff_time = d_ptr2->time - d_ptr1->time;
	coinc_ptr->source_pos = d_ptr1->source;
	//recorded
	return 0;
}

int record_compact_nontof_coinc(GATE_coincidence_event *d_ptr, nontof_coinc_event *coinc_ptr){
	//record this coinc
	coinc_ptr->crystal_index_1 = encode_crystal_mCT_Insert(d_ptr->single_ID_1);
	coinc_ptr->crystal_index_2 = encode_crystal_mCT_Insert(d_ptr->single_ID_2);
	coinc_ptr->time_1 = d_ptr->time_1;
	//recorded
	return 0;
}

int record_compact_tof_coinc(GATE_coincidence_event *d_ptr, tof_coinc_event *coinc_ptr){
	//record this coinc
	coinc_ptr->crystal_index_1 = encode_crystal_mCT_Insert(d_ptr->single_ID_1);
	coinc_ptr->crystal_index_2 = encode_crystal_mCT_Insert(d_ptr->single_ID_2);
	coinc_ptr->time_1 = d_ptr->time_1;
	coinc_ptr->diff_time = (d_ptr->time_2 - d_ptr->time_1)*1.5*100000000;
	//recorded
	return 0;
}

int sort_gate_coincidence_with_insert_non_tof(string inputfilename, vector<nontof_coinc_event> &list_data, float time_start, float time_end, float add_time_tag_offset){
	struct stat64  file_status;
	GATE_coincidence_event *gate_coincidence_data;
	int num_coincidence;
	FILE *inputfile;
	size_t t;
	__int64  size;


	printf("Size of gate_coinc_event structure is %d\n", (int)sizeof(GATE_coincidence_event));
	printf("Following coincidence files are being processed:\n");
	if (stat64(inputfilename.c_str(), &file_status) < 0) {
		printf("Read file error: %s \n", inputfilename.c_str());
		return 0;
	}

	num_coincidence = 0;

	size = file_status.st_size;
	printf("%I64d File size in byte \n", size);
	size /= 264;
	if (size >= INT_MAX){
		printf("File \"%s\" length is tooooooooooo long. This program only supports coincidence events fewer than %d. Current file has %lld. Program stopped\n", inputfilename.c_str(), INT_MAX, size);
		return 0;
	}

	num_coincidence = (int)size;
	printf("%d coincidence events in file: %s\n", num_coincidence, inputfilename.c_str());


	//Coincidence data memory
	gate_coincidence_data = (GATE_coincidence_event*)malloc(num_coincidence*sizeof(GATE_coincidence_event));

	//Read in coincidence data
	printf("Opening list mode file, %s \n", inputfilename.c_str());
	if ((inputfile = fopen(inputfilename.c_str(), "rb")) == NULL) {
		printf("Error: Could not read input file. \n");
		return 0;
	}

	t = fread(gate_coincidence_data, sizeof(GATE_coincidence_event), num_coincidence, inputfile);
	if (num_coincidence == (int)t){
		printf("%s read successfully.\n", inputfilename.c_str());
	}
	else{
		printf("%s read error. %d elements were read\n", inputfilename.c_str(), (int)t);
		return 0;
	}
	fclose(inputfile);

	int eventID;
	float time_tag;
	for (eventID = 0; eventID < num_coincidence; eventID++){
		time_tag = (float)gate_coincidence_data[eventID].time_1;

		if (eventID < MAX_COINC && time_tag<time_end && time_tag>time_start){

			nontof_coinc_event coinc;
			record_compact_nontof_coinc(&gate_coincidence_data[eventID], &coinc);

			list_data.push_back(coinc);
		}
	}
	free(gate_coincidence_data);

}


int sort_gate_coincidence_with_insert_tof(string inputfilename, vector<tof_coinc_event> &list_data, float time_start, float time_end, float add_time_tag_offset){
	struct stat64  file_status;
	GATE_coincidence_event *gate_coincidence_data;
	int num_coincidence;
	FILE *inputfile;
	size_t t;
	__int64  size;


	printf("Size of gate_coinc_event structure is %d\n", (int)sizeof(GATE_coincidence_event));
	printf("Following coincidence files are being processed:\n");
	if (stat64(inputfilename.c_str(), &file_status) < 0) {
		printf("Read file error: %s \n", inputfilename.c_str());
		return 0;
	}

	num_coincidence = 0;

	size = file_status.st_size;
	printf("%I64d File size in byte \n", size);
	size /= 264;
	if (size >= INT_MAX){
		printf("File \"%s\" length is tooooooooooo long. This program only supports coincidence events fewer than %d. Current file has %lld. Program stopped\n", inputfilename.c_str(), INT_MAX, size);
		return 0;
	}

	num_coincidence = (int)size;
	printf("%d coincidence events in file: %s\n", num_coincidence, inputfilename.c_str());


	//Coincidence data memory
	gate_coincidence_data = (GATE_coincidence_event*)malloc(num_coincidence*sizeof(GATE_coincidence_event));

	//Read in coincidence data
	printf("Opening list mode file, %s \n", inputfilename.c_str());
	if ((inputfile = fopen(inputfilename.c_str(), "rb")) == NULL) {
		printf("Error: Could not read input file. \n");
		return 0;
	}

	t = fread(gate_coincidence_data, sizeof(GATE_coincidence_event), num_coincidence, inputfile);
	if (num_coincidence == (int)t){
		printf("%s read successfully.\n", inputfilename.c_str());
	}
	else{
		printf("%s read error. %d elements were read\n", inputfilename.c_str(), (int)t);
		return 0;
	}
	fclose(inputfile);

	int eventID;
	float time_tag;
	for (eventID = 0; eventID < num_coincidence; eventID++){
		time_tag = (float)gate_coincidence_data[eventID].time_1;

		if (eventID < MAX_COINC && time_tag<time_end && time_tag>time_start){

			tof_coinc_event coinc;
			record_compact_tof_coinc(&gate_coincidence_data[eventID], &coinc);

			list_data.push_back(coinc);
		}
	}
	free(gate_coincidence_data);

}



int main(int argc, char * argv[]){

	FILE* fid;
	unsigned int total_num_events;

	string config_filename;
	string type_of_data;
	string data_path;
	int num_files;
	float time_seperation;
	string data_list_filename;
	string outputfilename;

	vector<string> filenames;
	vector<time_period> time_periods;

	if (argc != 2)
	{
		Usage(argv[0]);
		return 0;
	}

	// ============================= read in config file ============================= 
	config_filename = argv[1];
	Config config(config_filename);
	config.GetValue<string>("DAT_PATH", data_path);
	config.GetValue<int>("Number of files", num_files);
	config.GetValue<string>("Type of data", type_of_data);
	config.GetValue<string>("Data list file", data_list_filename);
	config.GetValue<float>("Time seperation", time_seperation);
	config.GetValue<string>("List mode data output file name", outputfilename);

	// ============================= read in data list file ==========================
	ifstream myfile(data_path + data_list_filename);
	if (!myfile.is_open()){
		std::cout << "Unable to open data list file: " << data_list_filename << endl;
	}
	else{
		string line;
		string name;
		time_period tp;
		for (int i = 0; i < num_files; i++){
			getline(myfile, line);
			stringstream ss(line);
			ss >> name >> tp.t_start >> tp.t_end;
			std::cout << name << " list mode events from " << tp.t_start << " to " << tp.t_end << endl;
			filenames.push_back(name);
			time_periods.push_back(tp);
		}
		myfile.close();
	}


	vector<tof_coinc_event> tof_list_data;

	for (int i = 0; i < num_files; i++){
		std::cout << endl;
		std::cout << "Now processing file " << i + 1 << " of " << num_files << endl;
		sort_gate_coincidence_with_insert_tof(data_path + filenames[i], tof_list_data, time_periods[i].t_start, time_periods[i].t_end, i*time_seperation);

	}

	total_num_events = tof_list_data.size();
	std::cout << "\nTotal number of coincidence events recorded: " << total_num_events << endl;

	fid = fopen((data_path + outputfilename).c_str(), "wb");

	int max_write_events = MAX_WRITE_BYTES / sizeof(tof_coinc_event);
	size_t t = 0;
	unsigned int events_written = 0;
	while (total_num_events - events_written > max_write_events){
		t += fwrite(&tof_list_data[events_written], sizeof(tof_coinc_event), max_write_events, fid);
		events_written += max_write_events;
	}
	t += fwrite(&tof_list_data[events_written], sizeof(tof_coinc_event), total_num_events - events_written, fid);

	if (total_num_events == t){
		printf("%d MB data was written successfully.\n", sizeof(tof_coinc_event)*total_num_events / (1024 * 1024));
	}
	else{
		printf("%s written with error. %d elements were written.\n", (data_path + outputfilename).c_str(), t);
	}
	fclose(fid);

}






