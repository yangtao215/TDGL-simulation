#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(int argc, char** argv)
{
	FILE *fp,*fpw,*fpn;
	int i,j,k,Nx,Ny,Nz;
	long ijk;
	double pha,phb,wa,wb,dx,dy,dz;
	char in[100], out[100], *fname, comment[201];
	double *phasav;

	fpn=fopen("para_out.txt","r");
	fscanf(fpn,"%d,%d,%d,%lf,%lf,%lf\n",&Nx,&Ny,&Nz,&dx,&dy,&dz);
//	fscanf(fpn,"%f,%f,%f",&dx,&dy,&dz);
//	printf("%d, %d\n",Nx,Ny);
printf("%d,%d,%d,%lf,%lf,%lf\n",Nx,Ny,Nz,dx,dy,dz);
	fclose(fpn);

	printf("%d\n", argc);

	if(argc==1)
	{
		sprintf(in, "pha.dat");
		sprintf(out, "out.vtk");
	}
	else if(argc==2)
	{
		fname = argv[1];
		sprintf(in, "%s", fname);
		sprintf(out, "out.vtk");
	}
	else if(argc==3)
	{
		fname = argv[1];
		sprintf(in, "%s", fname);
		fname = argv[2];
		sprintf(out, "%s", fname);
	}
	
	if((fp=fopen(in,"r"))==NULL)
	{
		fprintf(stderr,"ERROR: Cannot open file %s\n",in);
                exit(1);
	}

//	fscanf(fp, "Nx=%d, Ny=%d, Nz=%d\n",&Nx,&Ny,&Nz);
//	fscanf(fp, "dx=%lf, dy=%lf, dz=%lf\n",&dx,&dy,&dz);
//	printf("Nx=%d, Ny=%d, Nz=%d\n",Nx,Ny,Nz);
//	printf("dx=%lf, dy=%lf, dz=%lf\n",dx,dy,dz);

	if((fpw=fopen(out,"w"))==NULL)
	{
                fprintf(stderr,"ERROR: Cannot open file %s\n",out);
                exit(1);
	}

	fprintf(fpw,"# vtk DataFile Version 3.0\n");
	fprintf(fpw,"# comment: converted from file %s\n", in);
	fprintf(fpw,"ASCII\n");
	fprintf(fpw,"DATASET STRUCTURED_POINTS\n");
	fprintf(fpw,"DIMENSIONS %d %d %d\n",Nx,Ny,Nz); //don't exchange Nx and Nz
	fprintf(fpw,"ORIGIN 0.000000 0.000000 0.000000\n");
	fprintf(fpw,"SPACING %.6lf %.6lf %.6lf\n\n",dx,dy,dz);//don't exchange dx and dz
	fprintf(fpw,"POINT_DATA %ld\n",Nx*Ny*Nz);
	fprintf(fpw,"SCALARS scalars float\n");
	fprintf(fpw,"LOOKUP_TABLE default\n\n");

	phasav = (double *)malloc(sizeof(double)*Nx*Ny*Nz);

	for(i=0;i<Nx;i++)for(j=0;j<Ny;j++)for(k=0;k<Nz;k++)
	{
		ijk=(long)((k*Ny+j)*Nx+i);
		fscanf(fp, "%lf", &pha);
		phasav[ijk]=pha;
	}

	for(k=0;k<Nz;k++)
	for(j=0;j<Ny;j++)
		{
		for(i=0;i<Nx;i++)
			{
			ijk=(long)((k*Ny+j)*Nx+i);
			//fscanf(fp,"%lf %lf %lf %lf",&pha,&phb,&wa,&wb);
			fprintf(fpw,"%7.4lf ",phasav[ijk]);
			}
		fprintf(fpw,"\n");
		}
	printf("%f\n",pha);		

	fclose(fp);
	fclose(fpw);
	free(phasav);
}

