#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

int comp (const void * elem1, const void * elem2) 
{
    double f = *((double*)elem1);
    double s = *((double*)elem2);
    if (f >= s) return  1;
    if (f < s) return -1;
    return 0;
}

int dis_comp (const void * elem1, const void * elem2) 
{
    double f = *((double*)elem1);
    double s = *((double*)elem2);
    if (f >= s) return  -1;
    if (f < s) return 1;
    return 0;
}

int main(int argc, char ** argv)
{
	FILE * fp = fopen("sorted_list.txt", "r");

	float FPR_RANGE = 1.0;
	float FPR_STEP = 0.001;

	float * fpr_vec = (float*)malloc(sizeof(float)*(FPR_RANGE/FPR_STEP));
	float * tpr_vec_raisd = (float*)malloc(sizeof(float)*(FPR_RANGE/FPR_STEP));
	float * tpr_vec_omegaplus = (float*)malloc(sizeof(float)*(FPR_RANGE/FPR_STEP));
	float * tpr_vec_sweed = (float*)malloc(sizeof(float)*(FPR_RANGE/FPR_STEP));
	float * tpr_vec_sweepfinder2 = (float*)malloc(sizeof(float)*(FPR_RANGE/FPR_STEP));

	{
		int i, size1 = 0;
		fscanf(fp, "%d", &size1);

		double * vec1 = (double *) malloc(sizeof(double)*size1);
		for(i=0;i<size1;i++)
		{
			char string[1000];
			fscanf(fp, "%s", string);
			vec1[i] = (double)atof(string);
		}
	
		//qsort (vec1, size1, sizeof(double), comp);

		int size2 = 0;
		fscanf(fp, "%d", &size2);

		double * vec2 = (double *) malloc(sizeof(double)*size2);
		for(i=0;i<size2;i++)
		{
			char string[1000];
			fscanf(fp, "%s", string);
			vec2[i] = (double)atof(string);
		}
		//qsort (vec2, size2, sizeof(double), comp);
		float fpr_step = FPR_STEP;

		float fpr_check = 0.0;
	
		int index=0;
		while(fpr_check<=FPR_RANGE)
		{
			int fpr_index = (int)(fpr_check * size1);
			double fpr_thres = vec1[size1-fpr_index-1];

			//printf("index %d: %.7f ->", fpr_index, fpr_thres);

			int z=0, zstop=0;
			for(z=0;z<size2;z++)
			{
				zstop = z;
				
				if(vec2[z]>fpr_thres)
				{
					//zstop = z-1;
					break;
				}
			}
			fpr_vec[index] = fpr_check;
			tpr_vec_raisd[index++] = ((float)(size2-zstop)/size2);
			//printf("%f\t%f\n", fpr_check, ((float)(size2-zstop)/size2));
			fpr_check += fpr_step;	
		}
	}

	//printf("\n\n");
	{
		int i, size1 = 0;
		fscanf(fp, "%d", &size1);

		double * vec1 = (double *) malloc(sizeof(double)*size1);
		for(i=0;i<size1;i++)
		{
			char string[1000];
			fscanf(fp, "%s", string);
			vec1[i] = (double)atof(string);
		}

		//qsort (vec1, size1, sizeof(double), comp);

		int size2 = 0;
		fscanf(fp, "%d", &size2);

		double * vec2 = (double *) malloc(sizeof(double)*size2);
		for(i=0;i<size2;i++)
		{
			char string[1000];
			fscanf(fp, "%s", string);
			vec2[i] = (double)atof(string);
		}
		//qsort (vec2, size2, sizeof(double), comp);
		float fpr_step = FPR_STEP;

		float fpr_check = 0.0;
		int index=0;
		while(fpr_check<=FPR_RANGE)
		{
			int fpr_index = (int)(fpr_check * size1);
			double fpr_thres = vec1[size1-fpr_index-1];

			//printf("index %d: %.7f ->", fpr_index, fpr_thres);

			int z=0, zstop=0;
			for(z=0;z<size2;z++)
			{
				zstop = z;
				
				if(vec2[z]>fpr_thres)
				{
					//zstop = z-1;
					break;
				}
			}
		
			fpr_vec[index] = fpr_check;
			tpr_vec_omegaplus[index++] = ((float)(size2-zstop)/size2);
			//printf("%f\t%f\n", fpr_check, ((float)(size2-zstop)/size2));
			fpr_check += fpr_step;	
		}
	}
	
	//printf("\n\n");
	{
		int i, size1 = 0;
		fscanf(fp, "%d", &size1);

		double * vec1 = (double *) malloc(sizeof(double)*size1);
		for(i=0;i<size1;i++)
		{
			char string[1000];
			fscanf(fp, "%s", string);
			vec1[i] = (double)atof(string);
		}

		//qsort (vec1, size1, sizeof(double), comp);

		int size2 = 0;
		fscanf(fp, "%d", &size2);

		double * vec2 = (double *) malloc(sizeof(double)*size2);
		for(i=0;i<size2;i++)
		{
			char string[1000];
			fscanf(fp, "%s", string);
			vec2[i] = (double)atof(string);
		}
		//qsort (vec2, size2, sizeof(double), comp);
		float fpr_step = FPR_STEP;

		float fpr_check = 0.0;
		int index=0;
		while(fpr_check<=FPR_RANGE)
		{
			int fpr_index = (int)(fpr_check * size1);
			double fpr_thres = vec1[size1-fpr_index-1];

			//printf("index %d: %.7f ->", fpr_index, fpr_thres);

			int z=0, zstop=0;
			for(z=0;z<size2;z++)
			{
				zstop = z;
				
				if(vec2[z]>fpr_thres)
				{
					//zstop = z-1;
					break;
				}
			}
			fpr_vec[index] = fpr_check;
			tpr_vec_sweed[index++] = ((float)(size2-zstop)/size2);
			//printf("%f\t%f\n", fpr_check, ((float)(size2-zstop)/size2));
			fpr_check += fpr_step;	
		}
	}
	

	//printf("\n\n");
	{
		int i, size1 = 0;
		fscanf(fp, "%d", &size1);

		double * vec1 = (double *) malloc(sizeof(double)*size1);
		for(i=0;i<size1;i++)
		{
			char string[1000];
			fscanf(fp, "%s", string);
			vec1[i] = (double)atof(string);
		}

		//qsort (vec1, size1, sizeof(double), comp);

		int size2 = 0;
		fscanf(fp, "%d", &size2);

		double * vec2 = (double *) malloc(sizeof(double)*size2);
		for(i=0;i<size2;i++)
		{
			char string[1000];
			fscanf(fp, "%s", string);
			vec2[i] = (double)atof(string);
		}
		//qsort (vec2, size2, sizeof(double), comp);
		float fpr_step = FPR_STEP;

		float fpr_check = 0.0;
		int index=0;
		while(fpr_check<=FPR_RANGE)
		{
			int fpr_index = (int)(fpr_check * size1);
			double fpr_thres = vec1[size1-fpr_index-1];

			//printf("index %d: %.7f ->", fpr_index, fpr_thres);

			int z=0, zstop=0;
			for(z=0;z<size2;z++)
			{
				zstop = z;
				
				if(vec2[z]>fpr_thres)
				{
					//zstop = z-1;
					break;
				}
			}
			fpr_vec[index] = fpr_check;
			tpr_vec_sweepfinder2[index++] = ((float)(size2-zstop)/size2);
			//printf("%f\t%f\n", fpr_check, ((float)(size2-zstop)/size2));
			fpr_check += fpr_step;	
		}
	}
	
	fclose (fp);

	int i;
	fp = fopen("tpr_table.txt", "w");
	fprintf(fp, "%f\t%f\t%f\t%f\t\%f\n", 0.0, 0.0, 0.0, 0.0, 0.0);
	for(i=0;i<(FPR_RANGE/FPR_STEP);i++)
		fprintf(fp, "%f\t%f\t%f\t%f\t\%f\n", fpr_vec[i], tpr_vec_raisd[i], tpr_vec_omegaplus[i], tpr_vec_sweed[i], tpr_vec_sweepfinder2[i]);
	
	fprintf(fp, "%f\t%f\t%f\t%f\t\%f\n", 1.0, 1.0, 1.0, 1.0, 1.0);
	fclose(fp);

	//fscanf(fp, "%d", &size);
	//printf("next size %d\n", size);
	

	            			
}
