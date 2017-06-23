#include <stdlib.h>
#include <stdio.h>
#include <string.h>


int main(int argc, char ** argv)
{
	FILE * fp = fopen("command_line.txt", "r");
	FILE * fp2 = fopen("command_line2.txt", "w");

	char tstring[1000], tstring2[1000];
	int string_cnt=0;
	while(fscanf(fp, "%s", tstring)!=EOF)
	{
		//strncpy(tstring2, tstring, strlen(tstring));
		//printf("%s\n", tstring);
		if(strlen(tstring)>=10)
		{
			//printf("%s\n", tstring);
			tstring[10] = 0;

		}

		//string_cnt++;
		//if(string_cnt==15)
		//	fprintf(fp2, "\n");

		fprintf(fp2, "%s ", tstring);
	}


	fclose(fp);
	fclose(fp2);
}
