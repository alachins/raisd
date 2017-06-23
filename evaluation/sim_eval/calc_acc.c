#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char ** argv)
{
	// Read selection locations
	FILE * fp = fopen("sorted_list.txt", "r");
	int size = 0;
	fscanf(fp, "%d", &size);
	float target = 0.0;
	fscanf(fp, "%f", &target);
	float length = 0.0;
	fscanf(fp, "%f", &length);
	float threshold = 0.0;
	fscanf(fp, "%f", &threshold);
	
	int success_cnt=0;
	float value = 0.0;
	double sum_dist = 0.0;
	int i;
	for(i=0;i<size;i++)
	{
		fscanf(fp, "%f", &value);
		
		value = value - target;
		value = abs(value);
		sum_dist += (double) value;

		if(value<=threshold)
			success_cnt++;		
	}

	double average_distance = sum_dist / size;
	fprintf(stdout, "Average reported distance to Target: %.3f (%.5f%%)\n", average_distance, (average_distance/length)*100.0);
	fprintf(stdout, "Error distance threshold: %.0f (%.5f%%)\n", threshold, (threshold/length)*100.0);
	fprintf(stdout, "Success rate (reported locations at distance to selection target < threshold): %.3f (%.4f%%)\n", ((float)success_cnt)/((float)size), (((float)success_cnt)/((float)size))*100.0); 		fclose(fp);

	// Read neutral scores
	fp = fopen("sorted_list2.txt", "r");
	size = 0;
	fscanf(fp, "%d", &size);
	int fpr_index = 0.0;
	fscanf(fp, "%d", &fpr_index);
	float fpr_neutral_score = 0.0;
	for(i=0;i<size;i++)
	{
		fscanf(fp, "%f", &value);
		if(i==size-fpr_index)
		{
			fpr_neutral_score = value;
			break;
		}
	}
	fclose(fp);

	// Read selection scores
	fp = fopen("sorted_list3.txt", "r");
	size = 0;
	fscanf(fp, "%d", &size);
	int tpr_cnt=0;
	for(i=0;i<size;i++)
	{
		fscanf(fp, "%f", &value);
		if(value>=fpr_neutral_score)
		{
			tpr_cnt++;
		}
	}
	fprintf(stdout, "True Positive Rate for 5%% False Positive Rate: %f%%\n", ((float)tpr_cnt/(float)size)*100.0);
	fclose(fp);

	fp = fopen(argv[1], "a");

	if(!strcmp(argv[2], "SweeD"))
		fprintf(fp, "[%s]\t\tDistanceError: %f %% | SuccessRate: %f %% | TPR: %f %% | ", argv[2], (double)((average_distance/length)*100.0), (double)((((float)success_cnt)/((float)size))*100.0), (double)(((float)tpr_cnt/(float)size)*100.0));
	else
		fprintf(fp, "[%s]\tDistanceError: %f %% | SuccessRate: %f %% | TPR: %f %% | ", argv[2], (double)((average_distance/length)*100.0), (double)((((float)success_cnt)/((float)size))*100.0), (double)(((float)tpr_cnt/(float)size)*100.0));

	fclose (fp);            			
}
