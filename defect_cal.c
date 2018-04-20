#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define Pai (4.0*atan(1.0))	//3.141592653589
#define L0 9.3842
//#define Nx 289 
//#define Ny 473

#define Pi 3.1415926
#define Ley 1.732/2*L0+4.0
#define Lex 0.5*L0+4.0
#define Lnex 1.5*L0+4.0
#define Ns 80

/*data under this line may change,it should be redefine or read form file*/ 
#define Nsx 30.50  /*equal to Lx given in para.txt*/
#define dx 0.50	 /*equal to dx in tdgl.c*/
#define oneangle 72.0  /*72.0 equal to 360/n,n indicate n-side poly */

int Nx, Ny, Np, Ntri;
int *bd;

double *xpn, *ypn, Lx, Ly;


//int Locate_triag(double xi, double yi, int n1, 
//	int n2, int n3, double *h1, double *h2, double *h3);

//double tri_ang(double xa, double ya, double xb, 
//	double yb, double xp, double yp, double *hh);

//void ocf(double *psi);
//void ocf_bd();
//int change_n(int n);
void defect_N();
void surf2d();
int *pt;

main()
{
	double x0, y0, xi, yi, h1, h2, h3, deltaL;
//        double rx,ry,rxn,ryn;
	FILE *fp,*fpa;
	int nnew, i, j, k, n, n1, n2, n3, tag;
	double psiTot, psiR, psiI, dri, dxi, dyi, tht;
	double *Psi, *PsiB, psi_ij, psiB_ij;
	int *triag;
//	double nxpn,nypn;


//	Np = 289;
//	Npnew = 1076;
//	Ntri = 2132;

//	Nx = 360;
//	Ny = 468;
//	Np = 810;
//	Npnew = 1182;
//	Ntri = 2344;
//	printf("Input Nx, Ny, Np, Npnew, Ntri: \n");
	
	
	fp = fopen("paraNxyp.dat", "r");
	fscanf(fp, "%d, %d, %d", &Nx, &Ny, &Np);	
	fclose(fp);
	fp = fopen("paraN.txt","r");
	fscanf(fp,"%d",&Ntri);
	fclose(fp);
	printf("\n Nx = %d, Ny = %d\n", Nx, Ny);
	printf(" Np = %d,\n", Np);
	printf(" Ntri = %d\n", Ntri);
	Lx = Nx*1.0;
	Ly = Ny*1.0;


        pt = (int *)malloc(sizeof(int)*Np);
	triag = (int *)malloc(sizeof(int)*Ntri*3);
	Psi = (double *)malloc(sizeof(double)*Np);
	PsiB = (double *)malloc(sizeof(double)*Np);
	xpn = (double *)malloc(sizeof(double)*Np);
	ypn = (double *)malloc(sizeof(double)*Np);
	bd = (int *)malloc(sizeof(int)*Np*Ns);	//bd[Np][9], bd[i][0] save the number
					// of bonds on grid point of i
	for(i=0; i<Np; i++)
	{
		bd[i*Ns] = 0;	//initial:  bd[i][0]=0;
	}
	fp = fopen("points.txt", "r");
	for(i=0; i<Np; i++)
	{
		fscanf(fp, "%lf  %lf", &xpn[i], &ypn[i]);

	}
	fclose(fp);

	fp = fopen("points.delaunay.txt", "r");

	/**** Determine bonds on each point ****/
	for(i=0; i<Ntri; i++)
	{
		fscanf(fp, "%d  %d  %d", &n1, &n2, &n3);
		n1--;
		n2--;
		n3--;
		//k=7939;
		//if(n1==k||n2==k||n3==k){printf("%d  %d  %d\n", n1, n2, n3);getchar();}
		triag[i*3] = n1;
		triag[i*3+1] = n2;
		triag[i*3+2] = n3;
	}
	fclose(fp);
	for(i=0; i<Ntri; i++)
	{
		n1 = triag[i*3];
		n2 = triag[i*3+1];
		n3 = triag[i*3+2];

		if(n1<Np)	//only consider the points inside the box of (0, Nx, 0, Ny)
		{
			tag = 0;
			for(k=0; k<bd[n1*Ns]; k++)
			{
				if(n2==bd[n1*Ns+k+1])tag = 1;
			}
			if(tag==0)
			{
				k = bd[n1*Ns];
				bd[n1*Ns+k+1] = n2;
				bd[n1*Ns]++;
			}
                	tag = 0;
                	for(k=0; k<bd[n1*Ns]; k++)
                	{
                        if(n3==bd[n1*Ns+k+1])tag = 1;
                	}
                	if(tag==0)
                	{
                        	k = bd[n1*Ns];
                        	bd[n1*Ns+k+1] = n3;
                        	bd[n1*Ns]++;
                	}
		}

		if(n2<Np)
		{
                	tag = 0;
                	for(k=0; k<bd[n2*Ns]; k++)
                	{
                        	if(n1==bd[n2*Ns+k+1])tag = 1;
                	}
                	if(tag==0)
                	{
                        	k = bd[n2*Ns];
                        	bd[n2*Ns+k+1] = n1;
                        	bd[n2*Ns]++;
                	}
                	tag = 0;
                	for(k=0; k<bd[n2*Ns]; k++)
                	{
                        	if(n3==bd[n2*Ns+k+1])tag = 1;
                	}
                	if(tag==0)
                	{
                        	k = bd[n2*Ns];
                        	bd[n2*Ns+k+1] = n3;
                        	bd[n2*Ns]++;
                	}
		}
		if(n3<Np)
		{
                	tag = 0;
                	for(k=0; k<bd[n3*Ns]; k++)
                	{
                        	if(n1==bd[n3*Ns+k+1])tag = 1;
                	}
                	if(tag==0)
                	{
                        	k = bd[n3*Ns];
                        	bd[n3*Ns+k+1] = n1;
                        	bd[n3*Ns]++;
                	}
                	tag = 0;
                	for(k=0; k<bd[n3*Ns]; k++)
                	{
                        	if(n2==bd[n3*Ns+k+1])tag = 1;
                	}
                	if(tag==0)
                	{
                        	k = bd[n3*Ns];
                        	bd[n3*Ns+k+1] = n2;
                        	bd[n3*Ns]++;
                	}
		}

//fpa=fopen("bd.txt","a");
//if(k>7){
//fprintf(fpa,"i=%d,k=%d,xpn=%lf,ypn=%lf\n",i,k,xpn[i],ypn[i]);
//fprintf(fpa,"%lf\n",xpn[bd[i*Ns+1]]);
////for(j=1;j<Ns;j++){
////fprintf(fpa,"%d\n",i*Ns+j);
////fprintf(fpa,"%lf,%lf",xpn[bd[i*Ns+j]],ypn[bd[i*Ns+j]]);
////}
//}
//fclose(fpa);

	}


//fpa=fopen("bd.txt","a");
//for(i=0;i<Np;i++)
//{
//if(k==13)
//fprintf(fpa,"%d\n",bd[i*Ns]);
//}
//fclose(fpa);

/*	
fpa=fopen("bd.dat","a");
for(i=1;i<Np;i++)
{
          	 rx=fabs(xpn[i]);
          	 rxn=fabs(Nx-rx);
          	 ry=fabs(ypn[i]);
         	 ryn=fabs(Ny-ry);
      		if((bd[i*Ns]>7)&&(rx>Lex)&&(ry>Ley)&&(rxn>Lex)&&(ryn>Ley))
		{
		fprintf(fpa,"i=%d,k=%d\n%lf,%lf\n",i,bd[i*Ns],xpn[i],ypn[i]);
		for(j=1;j<bd[i*Ns];j++)
		fprintf(fpa,"%lf,%lf\n",xpn[bd[i*Ns+j]],ypn[bd[i*Ns+j]]);
		}
	
}
fclose(fpa);
*/	
	
	
	for(i=0; i<Np; i++)
	{

		if(bd[i*Ns]>=Ns){
		printf("Error: %3d( %d ): ", i, bd[i*Ns]);
		for(k=1; k<=bd[i*Ns]; k++)
			printf(" %3d (%f, %f) \n", bd[i*Ns+k], xpn[bd[i*Ns+k]], ypn[bd[i*Ns+k]]);
		printf("\n");
		getchar();}
}

	/******* Compute defect numbers **********************/
	
	surf2d();
	
	defect_N();	


//	fclose(fp);

	free(triag);
	free(Psi);
	free(PsiB);
	free(bd);
	free(xpn);
	free(ypn);
        free(pt);
}

void defect_N()
{
        int i, k;
        int df5, df7, dftot;
	int nnedge;
        FILE *fpt;
        double rx1,rxn1,ry1,ryn1;
//        double Lex,Ley,Lney;
        int Nedge;
//        Lex=8/3*L0;
//        Ley=5/3*L0*1.732;
//        Lney=5/3*1.5*L0;

        df5 = 0;
        df7 = 0;
        dftot = 0;
        Nedge=0;
	nnedge=0;
	
	

	fpt=fopen("defect.dat","w");
        for(i=0; i<Np; i++)
        {
           rx1=fabs(xpn[i]);
           rxn1=fabs(Nx-xpn[i]);
           ry1=fabs(ypn[i]);
           ryn1=fabs(Ny-ypn[i]);
                 k=bd[i*Ns];
/*                 if((rx1<Lex)||(rxn1<Lex)||(ry1<Ley)||(ryn1<Ley))Nedge++;
                        else if((rx1<Lnex)||(rxn1<Lnex))
                        {
				nnedge++;
				printf("%d,%lf,%lf,%lf\n",i,xpn[i],Lnex,ypn[i]);
                                if(k!=5)
				{
				dftot++;
				fprintf(fpt,"edge defect is %d,%f,%f\n",i,rx1,ry1);
				printf("edge defect is %d,%f,%f,%f\n",i,rx1,Lnex,ry1);
				printf("nnedge is %d\n",nnedge);
				}
                        }*/
				
		if(pt[i]==2)Nedge++;					
//		printf("%d\n",nnedge);
		if(pt[i]==1){
                        if(k!=6)
				{
                 	               dftot++;
                        	       if(k==5)		{
								df5++;
//								printf("df5 %d,%f,%f\n\n",i,xpn[i],ypn[i]);
							}
			else if(k==7)
							{
								df7++;
//								printf("df7 %d,%f,%f\n\n",i,xpn[i],ypn[i]);
							}

                          }
					}
		}
//        fpt = fopen("defect.dat", "w");

	fprintf(fpt, "Nedge=%d dftot=%d df5=%d df7=%d Np=%d  rdf=%.3f\n",Nedge,dftot,df5,df7,Np,100.0*dftot/Np);
	fclose(fpt);	
	printf(" Nedge=%d\n dftot=%d\n df5=%d\n df7=%d\n df5+7=%d\n Np=%d,\n rdf57=%.3f\n rdftot=%.3f\n",Nedge,dftot,df5,df7,df5+df7,Np,100.0*(df5+df7)/Np,100.0*dftot/Np);
	
//	fprintf(fpt, "Nedge=%d dftot=%d df5=%d df7=%d Np=%d  rdf=%.3f\n",Nedge,dftot,df5,df7,Np,100.0*dftot/Np);
//	fclose(fpt);	
//	printf(" Nedge=%d,dftot=%d,df5=%d,df7=%d,Np=%d,rdf=%.3f\n",Nedge,dftot,df5,df7,Np,100.0*dftot/Np);
//      printf("\nThe total defects are %d with 5- and 7- defects:  %d,  %d,  %.3f\n", dftot, df5, df7, 100.0*dftot/Np);


// free(pt);

}

void surf2d()
{
	
	int i,ipn,n;
	double dist,dist0,xipn,yipn,xx,yy;
	double theta,r,rt;
	FILE *fp;

        for(i=0;i<Np;i++)pt[i]=0;

        dist0=(cos(oneangle/2.0*Pi/180))*Nsx*L0/dx;
        printf("dist=%lf\n",dist0);


        printf("oneangle=%lf\n",oneangle);

	for(ipn=0;ipn<Np;ipn++)
	{	
	xx=xpn[ipn]-Nsx*L0/dx;
	yy=ypn[ipn]-Nsx*L0/dx;

	xipn=xx*cos((90.0-oneangle/2.0)*Pi/180)-yy*sin((90.0-oneangle/2.0)*Pi/180);
	yipn=xx*sin((90.0-oneangle/2.0)*Pi/180)+yy*cos((90.0-oneangle/2.0)*Pi/180);


	r=sqrt(xipn*xipn+yipn*yipn);	
	theta = acos(xipn/r)*180/Pi;
	
	n = (int)(theta/oneangle);
	dist = r*cos((oneangle/2.0+oneangle*n-theta)*Pi/180);

	dist=fabs(dist);


	rt= dist - dist0;
	if(rt<0.0)
		{
			pt[ipn]=1;
			if(fabs(rt)<0.86602540378*L0*3/dx/4)pt[ipn]=2;
		}
	}
	
	
	
	fp=fopen("pt.txt","w");
	for(i=0;i<Np;i++)
		{
		fprintf(fp,"%lf,%lf,%d\n",xpn[i],ypn[i],pt[i]);
		}

fclose(fp);
}
