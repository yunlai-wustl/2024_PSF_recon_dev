#include "../Solution_Items/GATE_data_structure.h"
#include <stdio.h>
#include <stdlib.h>

#define MAX_COINC 5000000000

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

int main(int argc, char * argv[])
{

	int head1, module1, crystal1, head2, module2, crystal2;
	volume_id vol_id;
	float toa;
	int event_count = 0;

	int i, j, tp, tg, pt;
	unsigned char x[8];
	unsigned char XE_lookup_A[42] = { 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 7, 7, 8 };
	unsigned char XE_lookup_B[42] = { 3, 4, 5, 6, 7, 8, 9, 4, 5, 6, 7, 8, 9, 10, 5, 6, 7, 8, 9, 10, 11, 6, 7, 8, 9, 10, 11, 7, 8, 9, 10, 11, 8, 9, 10, 11, 9, 10, 11, 10, 11, 11 };
	int dea1, dea2, ax, bx, ay, by, xe;
	int blk1, blk2, scr1, scr2, cryID1, cryID2;
	int lm_event_ID[6];
	int profile[23][52];
	int profile1[13];
	int time2step;
	int step;
	long time = 0, tim;
	int time_tag_num = 0;
	int cnt = 0;
	int sign = 0;
	FILE *fp;
	FILE *output;
	size_t t;
	for (i = 0; i<23; i++){
		for (j = 0; j<52; j++){
			profile[i][j] = 0;
		}
	}

	if (argc != 3){ printf("usage: %s inputdatafilename outputdatafilename\n", argv[0]); exit(1); } // Usage
	if (NULL == (fp = fopen(argv[1], "rb"))){ exit(1); } // IF File name is not proper

	mini_coinc_event* mini_coinc_event_ptr = (mini_coinc_event*)malloc(MAX_COINC*sizeof(mini_coinc_event));

	while (!feof(fp) && (event_count <MAX_COINC)){
		fread(x, 1, 8, fp); // read "1 byte * 8" one time
		tp = x[3] & 0x40;  // Tag_64 (0/1)
		tp = tp / 64;      // Tag == (0/1)?

		// According to the Tag_64, we perform different data decoding

		if (tp == 1){ // TYPE 1: non-event Tag
			//	printf("Time marker");
			sign = 0;
			tg = x[5] & 0xe0; //Tag word, 100x xxxx Time Marker
			if (tg == 128){
				time_tag_num++;
				tim = ((x[5] & 0x1f) * 16777216 + x[4] * 65536 + x[1] * 256 + x[0]);
				time = tim;
			}  //end of if

		}
		else if (tp == 0){ // TYPE 0: EVENT DATA DECODING 
			pt = (x[7] & 0x40) / 64; //Prompt tag
			xe = (x[6] & 0x07) * 8 + (x[2] & 0x07);  // module pair
			dea1 = XE_lookup_A[xe - 1];
			dea2 = XE_lookup_B[xe - 1];
			ax = x[0];
			ay = x[1];
			bx = x[4];
			by = x[5];


			blk1 = dea1 * 4 + ax / 13;
			blk2 = dea2 * 4 + bx / 13;

			head1 = blk1 + 10;
			head2 = blk2 + 10;
			module1 = 3 - ay/13;
			module2 = 3 - by/13;
			
			crystal1 = (12 - ay % 13) * 13 + ax % 13;
			crystal2 = (12 - by % 13) * 13 + bx % 13;

			toa = ((float)time)*0.001f;
	
			if (event_count < MAX_COINC){
				vol_id.system = 0;
				vol_id.sector = head1;
				vol_id.module = module1;
				vol_id.submodule = 0;
				vol_id.crystal = crystal1;
				vol_id.layer = 0;
				mini_coinc_event_ptr[event_count].crystal_index_1 = encode_crystal_mCT_Insert(vol_id);
				
				vol_id.system = 0;
				vol_id.sector = head2;
				vol_id.module = module2;
				vol_id.submodule = 0;
				vol_id.crystal = crystal2;
				vol_id.layer = 0;
				mini_coinc_event_ptr[event_count].crystal_index_2 = encode_crystal_mCT_Insert(vol_id);

				mini_coinc_event_ptr[event_count].time_1 = toa;
				mini_coinc_event_ptr[event_count].event_id = event_count;
				mini_coinc_event_ptr[event_count].diff_time = 0.0f;
				mini_coinc_event_ptr[event_count].source_pos.x = 0.0f;
				mini_coinc_event_ptr[event_count].source_pos.y = 0.0f;
				mini_coinc_event_ptr[event_count].source_pos.z = 0.0f;

			}

			event_count++;
		}

	}

	printf("......%d coincidence events were returned.\n", event_count);
	
	
	output = fopen(argv[2], "wb");
	t = fwrite(mini_coinc_event_ptr, sizeof(mini_coinc_event), event_count, output);
	if (event_count == (int)t){
		//printf("%s containing %d coincidence events was written successfully.\n",coin_filename, num_coincidence);
	}
	else{
		printf("%s written with error. %d elements were written.\n", argv[2], (int)t);
	}



	return 0;
}