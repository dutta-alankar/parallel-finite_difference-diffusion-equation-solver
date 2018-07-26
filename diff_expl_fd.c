#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diff_helper.h"
#define  MASTER		0

double** Initialize(int r, int c, double* x, double* y)
{
	/*Sets up the Initial and Boundary Conditions*/
	double** u = array2D(r, c);
	double radius;
	for (int i=0; i<r; i++)
	{
		for(int j=0; j<c; j++)
		{
			radius = pow(pow(x[j],2)+pow(y[i],2),0.5);
			if (radius>=0.5 && radius<=0.8) u[i][j] = 1.0;
			else u[i][j] = 0.0;
		}
	}
	make_dir();
	//FILE *f = fopen("diff_data/0.data", "wb");
	//fwrite(u, sizeof(double**), sizeof(u), f);
	//fclose(f);
	return u;
}

int main (int argc, char *argv[])
{
	int quit = 0, r = atoi(argv[1]), c = atoi(argv[2]); //rXc discrete points	
	
	double** u; 
	int   numtasks, taskid; 
	MPI_Status status;

	/***** MPI Initialization *****/
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

	if (r % numtasks != 0 || c % numtasks != 0) 
	{
		if (taskid == MASTER) printf("Quitting...\nNumber of rows & columns tasks must be divisible by the number of processors deployed.\n");
		MPI_Abort(MPI_COMM_WORLD, quit);
		exit(0);
	}

	int rend = (int)r/numtasks; //For Domian Decomposition
	if (numtasks==1) rend--;
	//printf("rend=%d\n",rend);

	/*Setting up the Computational Domain*/
	double xmin=-1, xmax=1, ymin=-1, ymax=1;
	double* x = linspace(xmin,xmax,c);
	double* y = linspace(ymin,ymax,r);
	//double xyd = 0.;
	//if (taskid==MASTER)
	//for (int m=0; m<r; m++){for(int n=0; n<c; n++){printf("(%f,%f)",x[n],y[m]);}printf("\n");}
	printf("Proc %d: Computation Domain successfully set.\n",taskid);
	MPI_Barrier(MPI_COMM_WORLD);

	double** u_old; 
	double** u_new;
	if (taskid == MASTER)
	{	
		u = Initialize(r,c,x,y);//Initialize and set the Boundary for the entire computational domain (only stored in MASTER)
		save_array(u, r, c, 0);
		//double **u_new = Initialize(r,c,x,y);//?
		
		/*Vertical (row) domain decomposition scheme (send to worker)*/
		int seek = 0;
		if (numtasks!=1) printf("Proc %d: Sending Domain Decomposition....\n",taskid);
		for (int splitter=0; splitter<numtasks; splitter++)
		{	
			//printf("%d\n",splitter);
			seek = rend*splitter-1;
			if(splitter==0) 
			{
				u_old = array2D(rend+1,c);
				u_new = array2D(rend+1,c);
				for (int i=0; i<(rend+1); i++)
				{ 
					for (int j=0; j<c; j++)
					{
					     u_old[i][j]=u[i][j];
					     u_new[i][j]=u[i][j];
					}
				}
			}
			else if(splitter==numtasks-1) MPI_Send(&u[seek][0], (rend+1)*c, MPI_DOUBLE, splitter, 75+splitter, MPI_COMM_WORLD);
			else if(splitter!=MASTER && splitter!=numtasks-1 ) MPI_Send(&u[seek][0], (rend+2)*c, MPI_DOUBLE, splitter, splitter+75, MPI_COMM_WORLD); 
				
		}
		// tags used till now: 0 to numtasks-1
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	/*Initialize u for each process and Vertical (row) domain decomposition scheme (received from master)*/
	else if (taskid==numtasks-1)
	{
		u_old = array2D(rend+1,c);
		u_new = array2D(rend+1,c);
		MPI_Recv(&u_old[0][0], (rend+1)*c, MPI_DOUBLE, MASTER, 75+taskid, MPI_COMM_WORLD, &status);
		for (int i=0; i<(rend+1); i++) for (int j=0; j<c; j++)	     u_new[i][j]=u_old[i][j];
		printf("Proc %d: Receiving Domain Decomposition....\n",taskid);
	}	
	else if(taskid!=MASTER && taskid!=numtasks-1)
	{
		u_old = array2D(rend+2,c);
		u_new = array2D(rend+2,c);
		MPI_Recv(&u_old[0][0], (rend+2)*c, MPI_DOUBLE, MASTER, 75+taskid, MPI_COMM_WORLD, &status);
		for (int i=0; i<(rend+2); i++) for (int j=0; j<c; j++)	     u_new[i][j]=u_old[i][j];
		printf("Proc %d: Receiving Domain Decomposition....\n",taskid);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	/*Computation starts from here*/
	printf("Proc %d: Computation starts...\n",taskid);
	double uxx, uyy, a = 1.0;
	double tstop = 1.0;
	double delx = x[1]-x[0], dely = y[1]-y[0]; 
	double dt = 0.5 *(1/(2*a))*(pow(delx,2)*pow(dely,2))/(pow(delx,2)+pow(dely,2));
	int nt = (int)(1/dt), stoprow, startrow;
	int collectfreq = (int) (nt/10);
	stoprow = rend+1;
	startrow = 1;
	if (taskid == MASTER || taskid == numtasks-1) stoprow=rend;
	//if (taskid == numtasks-1) stoprow = rend*numtasks -1;
	int n; double t;
	char pcnt = '%';
	double disp;
	for (n=1, t=0.0; n<=nt; n++, t=t+dt)
	{
		for (int i=startrow; i<stoprow; i++)
		{
			for (int j=1; j<c-1; j++)
			{
				disp = u_old[i][j];
				uxx = (u_old[i+1][j] - 2*u_old[i][j] + u_old[i-1][j])/pow((x[1]-x[0]),2);
				uyy = (u_old[i][j+1] - 2*u_old[i][j] + u_old[i][j-1])/pow((y[1]-y[0]),2);
				disp = u_old[i][j+1];
				u_new[i][j] = dt*a*(uxx + uyy) + u_old[i][j];
				disp = u_old[i+1][j];
				//printf("%.2f\n",u_new[i][j]);
			}
		}
		/*Swap*/
		for (int i=startrow; i<stoprow; i++)
		{
			for (int j=1; j<c-1; j++)
			{
				disp = u_old[i][j];
				u_old[i][j] = u_old[i][j]+u_new[i][j];
				u_new[i][j] = u_old[i][j]-u_new[i][j];
				u_old[i][j] = u_old[i][j]-u_new[i][j];
				disp = u_old[i][j];
				if (disp>1) printf("Error! Check Convergence!");
			}
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		/*Communicate Domain edge values*/
		if (taskid==MASTER && numtasks!=1)
		{
			//printf("taskid inside compute: %d\n",taskid);
			MPI_Send(&u_old[rend-1][0], c, MPI_DOUBLE, taskid+1, 30, MPI_COMM_WORLD);
			MPI_Recv(&u_old[rend][0], c, MPI_DOUBLE, taskid+1, 30, MPI_COMM_WORLD, &status);
		}
		else if (taskid==numtasks-1 && numtasks!=1)
		{
			//printf("taskid inside compute:%d\n",taskid);
			MPI_Send(&u_old[1][0], c, MPI_DOUBLE, taskid-1, 30, MPI_COMM_WORLD);
			MPI_Recv(&u_old[0][0], c, MPI_DOUBLE, taskid-1, 30, MPI_COMM_WORLD, &status);
		}
		else if(taskid!=MASTER && taskid!=numtasks-1 && numtasks!=1)
		{
			//printf("taskid inside compute:%d\n",taskid);
			MPI_Send(&u_old[1][0], c, MPI_DOUBLE, taskid-1, 30, MPI_COMM_WORLD);
			MPI_Send(&u_old[rend][0], c, MPI_DOUBLE, taskid+1, 30, MPI_COMM_WORLD);
			MPI_Recv(&u_old[0][0], c, MPI_DOUBLE, taskid-1, 30, MPI_COMM_WORLD, &status);//printf("mark\n");
			MPI_Recv(&u_old[rend+1][0], c, MPI_DOUBLE, taskid+1, 30, MPI_COMM_WORLD, &status);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		/*Collect computed result to MASTER*/
		if(n%collectfreq==0 || n==nt)
		{
			if (taskid==MASTER) printf("t=%f (%.2f%c done)\n",t,(t/tstop)*100,pcnt);
			if (taskid==MASTER)
			{
				//MPI_Recv(&u[0][0], rend*c, MPI_DOUBLE, MASTER, 100, MPI_COMM_WORLD, &status);
				//MPI_Recv(&u[(numtasks-1)*rend][0], rend*c, MPI_DOUBLE, numtasks-1, 100, MPI_COMM_WORLD, &status);
				if (numtasks!=1)
				{ 
					for (int receiver=1; receiver<numtasks; receiver++)
						MPI_Recv(&u[receiver*rend][0], rend*c, MPI_DOUBLE, receiver, 111+receiver, MPI_COMM_WORLD, &status);
				}
				for (int i=1; i<rend; i++)
				{
					for (int j=1; j<c-1; j++) u[i][j] = u_old[i][j];
				}
			}
			else if (numtasks!=1)
			{
				//if (taskid==MASTER) MPI_Send(&u_old[0][0], rend*c, MPI_DOUBLE, MASTER, taskid, MPI_COMM_WORLD);
				MPI_Send(&u_old[1][0], rend*c, MPI_DOUBLE, MASTER, 111+taskid, MPI_COMM_WORLD);
			}
		MPI_Barrier(MPI_COMM_WORLD);
		if (taskid==MASTER) save_array(u, r, c, n);
		//MPI_Barrier(MPI_COMM_WORLD);
		}
	}

/*	
		double** u_old = array2D(r, c);

	
	double** u_old = array2D(r, c);
	for (int i=0; i<r; i++)
	{
		for(int j=0; j<c; j++)
		{
			printf("u_old[%d][%d]=%f   ",i,j,u_old[i][j]);
		}
	printf("\n");
	}
*/
	MPI_Finalize(); 
	return 0;
} 
	
