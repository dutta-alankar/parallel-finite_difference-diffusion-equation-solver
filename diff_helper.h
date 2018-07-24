#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>

double** array2D (int r, int c)
{
	/*Creates a 2D array*/
	double* B = (double *)malloc(r*c*sizeof(double));
	assert(NULL != B);
	double** A = (double **)malloc(sizeof(double *)*r);
	assert(A);
	A[0] = B;
	for (int i=1; i<r; i++) A[i] = A[i-1] + c;
	return A;
}



double* linspace(double start, double stop, int n)
{
	/*Similar to numpy linspace*/
	double delta = (stop-start)/(n-1);
	double* linsp = (double *)malloc(n*sizeof(double)); 
	assert(NULL != linsp);
	for (int i=0; i<n; i++) linsp[i] = start + i*delta;
	return linsp;
}

void make_dir()
{
	struct stat st = {0};

	if (stat("diff_out", &st) == -1) 
	mkdir("diff_out", 0700);
}


void save_array(double** a, int r, int c, int n)
{
	FILE *file;
	char str[10];
	char loc[50] = "diff_out/data.";
	sprintf(str, "%06d", n);
	strcat(loc,str);
	file=fopen(loc,"w");     

	for(int i=0; i<r; i++) 
	{
		for (int j=0; j<c; j++) fprintf(file,"%f ",a[i][j]);
		fprintf(file,"\n");
	}	
	fclose(file);
	//if (n==0)	for(int m=0; m<sizeof(data); m++) printf("%c",data[m]);
}

