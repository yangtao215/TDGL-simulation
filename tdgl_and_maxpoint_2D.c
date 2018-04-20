//Phase separation dynamics based upon TDGL (Model B)--------------------
//i
//Homo and copolymer blend and expand to the 4th order
//calculate bond length and theta automatically 

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

#define Nstep 2000000
#define thetaCritical 0.001

#define Fluct 0.040
#define Pi 3.1415926535898
#define N 3000
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

#define PIOZ (atan(1.0)*2.0)
#define PI   (2.0*PIOZ)

#define ITMAX 1000
#define EPS 1.0e-8
#define TOL 2.0e-8

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

#define CGOLD 0.3819660
#define ZEPS 1.0e-10

#define NDIMM 5

int nside;
int Nx, Ny, *pt;
long NxNy, npt,NstepCritical=0;
double fA, fB, fC,  alpha, beta, gamm, Ac, Ac1, Arand;
double *phiA, *phiAold, *mu, PHIA0, *phiH, *phiHold, *muH, PHIH0, PHIHcritical;
double b1,b2,b3,b4;
int *ing, *ips, *jng, *jps;
double Nsx;
double L0, Ls, rs;
double Dx, Dy, Dz, DeltaT;
double V0, *Hs, *indexx, nLs;
double C1, C2;

void surf2d();
void Init();
void cal_mu();		
		//period boundary conditions
void outDen(long step);

int findmaxpoint(long step,FILE *fp0);
int caltheta(long step,FILE *fp0);
double Laplc(int i0, int j0 ,int kk);
void search_edge();

int ncom;
int Npoint=0;
double pcom[NDIMM],xicom[NDIMM];
double tau,gamma0,xisqr;

void frprmn(double p[],int n,double ftol,int *iter,double *fret);
void linmin(double p[],double xi[],int n,double *fret);
void mnbrak(double *ax,double *bx,double *cx,double *fa,double *fb,double *fc);
double brent(double ax,double bx,double cx,double tol,double *xmin);
double fldim(double x);

int maxf();
void interP(double *xmax,double *ymax,double *phmax);

double func(double x[]);
void dfunc(double x[],double dx[]);
double func1(double a,double b,double c,double x0,double y0);

double xt[9],yt[9],fv[9],fva[4],fvb[5];
double f, ph[N][N];

double ph_nb(int j, int k);
int Nx, Ny, Nz;

struct Pt{
	double x;
	double y;
};

int on_caltheta=1;
int switch_theta = 1;
int switch_Npoint = 1;


FILE *fp0;
//fp0=fopen("cal_theta_step.dat","w");


main()
{
	int i, j, iseed;
	int in, ip, jn, jp, swch;
	long ij, istep;
	double phav, phhv;
	double Nxf, Nyf;
	double *etax, *etay;
	time_t ts,ts_new;
	FILE *fp;
//	FILE *fp0;

	iseed = time(&ts);
	srand48(iseed);
//	fp0 = fopen("cal_theta_step.dat","a");
        fp = fopen("para.txt", "r");

        fscanf(fp,"%lf, %lf, %lf, %lf",&fA,&fC,&Ac,&Ac1);
        fscanf(fp,"%lf, %lf, %lf, %lf",&b1,&b2,&b3,&b4);
        fscanf(fp,"%lf, %lf",&alpha,&PHIHcritical);
        fscanf(fp,"%lf, %lf",&C1,&C2);
        fscanf(fp,"%lf",&DeltaT);
        fscanf(fp,"%lf, %lf",&V0,&Arand);
        fscanf(fp,"%lf",&Nsx);
        fscanf(fp,"%lf",&L0);
	fscanf(fp,"%d",&nside);
        fclose(fp);
	rs = 0.150;
	Dx = 0.50;
	Dz = 1.0;
	Nxf = 2*Nsx * L0 / Dx;
	Nx = (int)(Nxf);
	if(Nxf-Nx*1.0>0.50)
		Nx++;
	Nx += 4;

/*        Nyf = sqrt(3)/2*Nsx * L0 / Dx;
        Ny = (int)(Nyf);
        if(Nyf-Ny*1.0>0.50)
                Ny++;
        Ny += 4;*/

	Ny=Nx;


	NxNy = (long) (Nx*Ny);

	Nz=1;

	fB=(1-fA)*(1-fC);
	fA=fA*(1-fC);

	beta = alpha*(fB - fA)/(fB + fA);
	gamm = beta*(fB - fA)/(fB + fA);
	PHIA0 = fA - fB;
	PHIH0 = fA+fB-PHIHcritical;

	fp = fopen("para_out.txt", "w");
	fprintf(fp, "\n  %d,  %d,  %d,  %lf,  %lf,  %lf\n", Nx, Ny, Nz,Dx,Dx,Dz);
	fprintf(fp, "\n fA = %f, fB = %f, fC = %f, Ac = %f, Ac1 = %f\n", fA, fB, fC, Ac,Ac1);
	fprintf(fp, "\n b1 = %f, b2 = %f, b3 = %f, b4 = %f\n", b1, b2, b3, b4);
	fprintf(fp, "\n alpha = %f,beta = %f, gamm = %f, PHIHcritical = %f\n", alpha, beta, gamm, PHIHcritical);
	fprintf(fp, "\n phiA0 = %f, phiH0 = %f\n", PHIA0, PHIH0);
	fprintf(fp, "\n C1 = %f , C2 = %f\n", C1, C2);
	fprintf(fp, "\n DeltaT = %f\n", DeltaT);
	fprintf(fp, "\n V0 = %f\n", V0);
	fprintf(fp, "\n NsX = %f\n", Nsx);
	fprintf(fp, "\n L0 = %.1f\n", L0);
	fprintf(fp, "\n rs = %.3f*L0\n", rs);
	fprintf(fp, "\n nside = %d\n",nside);
	fprintf(fp, "\n iseed = %d\n",iseed);
	fclose(fp);


	ing = (int *)malloc(sizeof(int)*Nx);
	ips = (int *)malloc(sizeof(int)*Nx);
	jng = (int *)malloc(sizeof(int)*Ny);
	jps = (int *)malloc(sizeof(int)*Ny);
	etax = (double *)malloc(sizeof(double)*NxNy);  //noise
	etay = (double *)malloc(sizeof(double)*NxNy);
	phiA = (double *)malloc(sizeof(double)*NxNy);
	phiAold = (double *)malloc(sizeof(double)*NxNy);
	mu = (double *)malloc(sizeof(double)*NxNy);	//chemical potential
	phiH = (double *)malloc(sizeof(double)*NxNy);
	phiHold = (double *)malloc(sizeof(double)*NxNy);
	muH = (double *)malloc(sizeof(double)*NxNy);	//chemical potential
	Hs = (double *)malloc(sizeof(double)*NxNy);
	pt = (int *)malloc(sizeof(int)*NxNy);
	indexx = (double *)malloc(sizeof(double)*NxNy);

	Dx = 0.50;
	Dy = 0.50;
	Dz=1;
	for(i=0; i<Nx; i++)
	for(j=0; j<Ny; j++){
		ij=(long)(i*Ny+j);
		Hs[ij]=0.0;
		pt[ij]=0;
	}

	surf2d();
	search_edge();
	/******** Producing Periodic Boundary Conditions ***********/

	ing[0] = Nx - 1;
	ips[0] = 1;
	ing[Nx-1] = Nx - 2;
	ips[Nx-1] = 0;
	for(i=1; i<=Nx-2; i++)
	{
		ing[i] = i - 1;
		ips[i] = i + 1;
	}

	jng[0] = Ny - 1;
	jps[0] = 1;
	jng[Ny-1] = Ny - 2;
	jps[Ny-1] = 0;
	for(j=1; j<=Ny-2; j++)
	{
		jng[j] = j - 1;
		jps[j] = j + 1;
	}

	Init();		//initializing densities with random thermal fluctuation

	/*********** Dynamic evolution ***********/

	for(istep = 1; istep<= Nstep; istep ++)
	{
		
		/*if(istep%1000==0)
		{
			printf(" %6d : \n", istep);
			fp = fopen("fet.dat", "a");
			fprintf(fp," Step %6d : \n", istep);
			fclose(fp);
		}*/

		for(i=0; i<Nx; i++)
			for(j=0; j<Ny; j++)
			{
				etax[i*Ny+j] = 2.0*(drand48()-0.50);
				etay[i*Ny+j] = 2.0*(drand48()-0.50);
			}

		cal_mu();

		phav = 0.0;
		phhv = 0.0;

		for(i=0; i<Nx; i++)
		for(j=0; j<Ny; j++){
			ij = (long) (i*Ny + j);
			if(pt[ij]){
				ip = ips[i];
				jp = jps[j];
				if((istep%10000)>9900||(istep%10000)==0||((istep<=10000)&&(istep%1000>900))||((istep%1000)==0)) 
					swch = 0;
				else 
					swch = 1;
		
				phiA[ij] += DeltaT*((Laplc(i, j, 0) - alpha*(phiAold[ij]-PHIA0)-beta*(phiHold[ij]-PHIH0))
					+swch*Arand*(etax[ip*Ny+j]-etax[i*Ny+j]+etay[i*Ny+jp]-etay[i*Ny+j]));
				phav += phiA[ij];

				phiH[ij] += DeltaT*((Laplc(i, j, 1)-gamm*(phiHold[ij]-PHIH0)-beta*(phiAold[ij]-PHIA0))
					+swch*Arand*(etax[ip*Ny+j]-etax[i*Ny+j]+etay[i*Ny+jp]-etay[i*Ny+j]));
				phhv += phiH[ij];
			}
		}
			phav /= npt;
			phav -= PHIA0;
			phhv /= npt;
			phhv -= PHIH0;

		for(i=0; i<Nx; i++)
		for(j=0; j<Ny; j++){
			ij = (long) (i*Ny + j);
			if(pt[ij]){
				phiA[ij] -= phav;
				phiH[ij] -= phhv;
			}
			phiAold[ij] = phiA[ij];
			phiHold[ij] = phiH[ij];
		}

//		if(istep<10000)
//		if(istep%1000==0)
//			outDen(istep);
		if(istep%10000==0)
			{
				outDen(istep);
				fp0=fopen("cal_theta_step.dat","a");
				fprintf(fp0,"\n%ld\t",istep);
				printf("\n%ld\t",istep);
				on_caltheta = findmaxpoint(istep,fp0);
				if (on_caltheta==1) {switch_theta=caltheta(istep,fp0);} else {fclose(fp0);continue;}
				fclose (fp0);
			}
		if ((istep==NstepCritical)&&(switch_theta==0)&&(switch_Npoint==0)) return 0;			
	}

	free(ing);
	free(ips);
	free(jng);
	free(jps);
	free(etax);
	free(etay);
	free(Hs);
	free(pt);
	free(indexx);
	free(mu);
	free(phiA);
	free(phiAold);
	free(muH);
	free(phiH);
	free(phiHold);
	fclose(fp0);
}

void surf2d()
{
	int i, j, ij, n;
	double rsigma;	//rsigma = 1.0/sigma
	double rt, xt, yt, xp, yp, r, dist, rt1, rt2;
	double x, y, theta,xx,yy;
	double lx, ly;
	FILE *fp,*fp1;

	double oneangle;

	rsigma = 2.0;
	npt = 0;

	lx = Nx*1.0;		//Nsx*Ls;
	ly = Ny*1.0;		//Nsy*sqrt(3.0)*Ls;
	rs = rs*L0;

	for(ij=0; ij<NxNy; ij++){
		mu[ij] = 0.0;
		muH[ij] = 0.0;
		Hs[ij] = 0.0;
	}

	oneangle = 360/nside;

	for(i=0; i<Nx; i++)
	for(j=0; j<Ny; j++){
		
//		i=(int)(ii*cos(18*Pi/180)-jj*sin(18*Pi/180));
//		j=(int)(ii*sin(18*Pi/180)+jj*cos(18*Pi/180));
		
		ij = i * Ny + j;
		xx = i - Nx/2 + 0.50;
		yy = j - Ny/2 + 0.50;

                x=xx*cos((90-oneangle/2)*Pi/180)-yy*sin((90-oneangle/2)*Pi/180);
                y=xx*sin((90-oneangle/2)*Pi/180)+yy*cos((90-oneangle/2)*Pi/180);



//		x = i-1.0;
//		y = j-1.0;

		r = sqrt(x*x+y*y);
		theta = acos(x/r)*180/Pi;
		n = (int)(theta/oneangle);
		dist = r*cos(((oneangle/2.0)+(oneangle)*n-theta)*Pi/180);
//		dist = r*cos((36-theta)*Pi/180);
		dist = fabs(dist);

		rt = (cos((oneangle/2.0)*Pi/180))*Nsx*L0 / Dx - dist;
//		rt1 = r*sin((60-theta)*Pi/180);
//		rt2 = r*sin(theta*Pi/180);

		if(rt > 0){
			pt[ij] = 1;
			npt++;
			if(rt<=2*rs)Hs[ij] = 0.50*V0*(tanh(rsigma*(-rt+rs))+1.0);
//			if(rt2<=2*rs)Hs[ij] = 0.50*V0*(tanh(rsigma*(-rt2+rs))+1.0);
//			if(rt1<=2*rs)Hs[ij] = 0.50*V0*(tanh(rsigma*(-rt1+rs))+1.0);
				
		}
//		printf("%lf\n", Hs[ij]);
//		printf("%ld\n", pt[ij]);
	}
	fp1=fopen("Hs.txt","w");	
	fp=fopen("pt.txt","w");
	for(i=0;i<Nx;i++)
	for(j=0;j<Ny;j++)
	{
                ij=i*Ny+j;
		fprintf(fp,"%d   ",pt[ij]);
		fprintf(fp1,"%lf   ",Hs[ij]);
		if(j==Ny-1)
		{
			fprintf(fp,"\n");
			fprintf(fp1,"\n");
		}
	}
	fclose(fp);
	fclose(fp1);
}

void search_edge()
{
	int i, j, ij, ij1, sign;
	int m, n, num;
	int dexx;
	FILE *fp;
	
	for(i = 1; i < Nx - 1; i++)
	for(j = 1; j < Ny - 1; j++){
		ij = i * Ny + j;
		indexx[ij] = 1.0;
		if(pt[ij]){
			for(m = -1; m < 2; m++)
			for(n = -1; n < 2; n++){
				ij1 = ij + m* Ny + n;
				if(!pt[ij1]){
					dexx = m*m + n*n;
					if(dexx==2)
						indexx[ij] -= 1.0/12.0;
					if(dexx==1)
						indexx[ij] -= 1.0/6.0;
				}		
			}
		}
	}
}			

void Init()
{
	int i, j;
	long ij;
	double ran, phiav, phihv;

	phiav = 0.0;
	phihv = 0.0;

	for(j=0; j<Ny; j++)
	for(i=0; i<Nx; i++){
		ij = (long) (i*Ny + j);
		if(pt[ij]){
			ran = drand48() - 0.50;
			phiA[ij] = ran*Fluct;
			phiav += phiA[ij];
			ran = drand48() - 0.50;
			phiH[ij] = ran*Fluct;
			phihv += phiH[ij];
		}
		else{
			phiA[ij] = 0.0;
			phiH[ij] = 0.0;
		}
	}

	phiav /= npt;
	phihv /= npt;

	for(i=0; i<Nx ;i++)
	for(j=0; j<Ny; j++){
		ij = (long)(i*Ny + j);
		if(pt[ij]){
			phiA[ij] = phiA[ij] - phiav + PHIA0;
			phiH[ij] = phiH[ij] - phihv + PHIH0;
		}
		phiAold[ij] = phiA[ij];
		phiHold[ij] = phiH[ij];
	}
}

void cal_mu()
{
	int i, j, in, ip, jn, jp;
	long ij;
	double sum;


	for(i=0; i<Nx; i++)
	for(j=0; j<Ny; j++){
		ij = (long) (i*Ny + j);

		if(pt[ij]){
			in = ing[i];
			ip = ips[i];
			jn = jng[j];
			jp = jps[j];

			mu[ij] = - Ac*tanh(phiAold[ij]) + phiAold[ij]
					+b1*phiHold[ij]
					-b2*phiHold[ij]*phiAold[ij]
					-b3*phiHold[ij]*(phiHold[ij]+phiHold[ij]*phiHold[ij]+3.0*phiAold[ij]*phiAold[ij])
					+b4*phiHold[ij]*phiHold[ij]*phiAold[ij];
	
			muH[ij] = - Ac1*tanh(phiHold[ij]) + phiHold[ij]
					+b1*phiAold[ij]
					-0.5*b2*phiAold[ij]*phiAold[ij]
					-b3*phiAold[ij]*(2.0*phiHold[ij]+3.0*phiHold[ij]*phiHold[ij]+phiAold[ij]*phiAold[ij])
					+b4*phiHold[ij]*phiAold[ij]*phiAold[ij];

			/************** Ohta expression for Laplacian ***************/
		
			//assume Dx=Dy=1

			sum = (phiAold[in*Ny+j]+phiAold[ip*Ny+j]+phiAold[i*Ny+jn]+phiAold[i*Ny+jp])/6.0;
			sum += (phiAold[in*Ny+jn]+phiAold[ip*Ny+jn]+phiAold[ip*Ny+jp]+phiAold[in*Ny+jp])/12.0;
			sum -= indexx[ij]*phiAold[ij];
			sum /= (Dx*Dx);

			mu[ij] -= C1*sum;

			sum=0;
			sum = (phiHold[in*Ny+j]+phiHold[ip*Ny+j]+phiHold[i*Ny+jn]+phiHold[i*Ny+jp])/6.0;
			sum += (phiHold[in*Ny+jn]+phiHold[ip*Ny+jn]+phiHold[ip*Ny+jp]+phiHold[in*Ny+jp])/12.0;
			sum -= indexx[ij]*phiHold[ij];
			sum /= (Dx*Dx);

			muH[ij] -= C2*sum;

			mu[ij] += Hs[ij];
			muH[ij] += Hs[ij];
		}
		else{
			mu[ij] = 0.0;	
			muH[ij] = 0.0;	
		}
	}
}

double Laplc(int i0, int j0, int kk)
{
	int in, ip, jn, jp;
	long ij;
	double lp;

	ij = i0*Ny + j0;

	in = ing[i0];
	ip = ips[i0];
	jn = jng[j0];
	jp = jps[j0];

	if(kk==0){
		lp = (mu[in*Ny+j0]+mu[ip*Ny+j0]+mu[i0*Ny+jn]+mu[i0*Ny+jp])/6.0;
		lp += (mu[in*Ny+jn]+mu[in*Ny+jp]+mu[ip*Ny+jn]+mu[ip*Ny+jp])/12.0;
		lp -= indexx[ij]*mu[i0*Ny+j0];
		lp /= (Dx*Dx);
	}else
	{
		lp = (muH[in*Ny+j0]+muH[ip*Ny+j0]+muH[i0*Ny+jn]+muH[i0*Ny+jp])/6.0;
		lp += (muH[in*Ny+jn]+muH[in*Ny+jp]+muH[ip*Ny+jn]+muH[ip*Ny+jp])/12.0;
		lp -= indexx[ij]*muH[i0*Ny+j0];
		lp /= (Dx*Dx);
	}

	return lp;
}

void outDen(long step)
{
	int i, j;
	long ij;
	FILE *fp;
	char fname[50];
	double phiAmax, phiHmax, phiAmin, phiHmin;
	double sum, sumH;
	double phit;

	phiAmax = 0.0;
	phiHmax = 0.0;
	phiAmin = 2.0;
	phiHmin = 2.0;

	fp = fopen("fname.txt", "w");
	fprintf(fp, "pha%ld.dat", step);
	fclose(fp);

	fp = fopen("fname.txt", "r");
	fscanf(fp, "%s", fname);
	fclose(fp);

	sum = 0;
	sumH = 0;
	fp = fopen(fname, "w");
	for(i=0; i<Nx; i++){
		for(j=0; j<Ny; j++){
			ij = (long) (i*Ny + j);
			fprintf(fp, "%.8e\n", 0.5*(phiA[ij]+phiH[ij]+PHIHcritical));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}



int findmaxpoint(long step,FILE *fp0)
{
	int i, j, k, ng, n, nnew,Npointold;
	int in, ip, jn, jp, kn, kp, jn2, jp2, kn2, kp2;
	double e1,e2,pha,phb,wa,wb;
	double xmax,ymax,phmax;
	double xi, yi, Lx, Ly, deltaL;
	double xpn[50000], ypn[50000];
	FILE *fp, *fp2;
	char fname[100];


	fp=fopen("para_out.txt","r");
	fscanf(fp,"%d,%d",&Nx,&Ny);
	printf("%d, %d\n",Nx,Ny);
	fclose(fp);


	fp = fopen("fname.txt", "w");
	fprintf(fp, "pha%ld.dat", step);
	fclose(fp);

	fp = fopen("fname.txt", "r");
	fscanf(fp, "%s", fname);
	fclose(fp);
	fp=fopen(fname,"r");

	fp2 = fopen("fname.txt", "w");
	fprintf(fp2, "paraN%ld.dat", step);
	fclose(fp2);

	fp2 = fopen("fname.txt", "r");
	fscanf(fp2, "%s", fname);
	fclose(fp2);

	fp2 = fopen(fname, "w");

	Lx = Nx*1.0;
	Ly = Ny*1.0;

	for(i=0;i<Nx;i++)
	for(j=0;j<Ny;j++)
	{
		fscanf(fp,"%lf ", &pha);
			ph[i][j]=pha;
	}
	fclose(fp);

	Npointold = Npoint;
	Npoint=0;

	fp = fopen("fname.txt", "w");
	fprintf(fp, "max2d%ld.dat", step);
	fclose(fp);

	fp = fopen("fname.txt", "r");
	fscanf(fp, "%s", fname);
	fclose(fp);
	fp=fopen(fname,"w");
	

	for(i=0;i<Nx;i++)
	for(j=0;j<Ny;j++){
		fv[0]=ph[i][j];
		in=i-1;
		ip=i+1;
		if(in<0)
			in+=Nx;
		else if(ip>=Nx)
			ip-=Nx;
		jn=j-1;
		jp=j+1;
		if(jn<0)
			jn+=Ny;
		else if(jp>=Ny)
			jp-=Ny;
		fv[1]=ph[in][j];
		fv[2]=ph[ip][j];
		fv[3]=ph[i][jn];
		fv[4]=ph[i][jp];

		fv[5] = ph[in][jn];
		fv[6] = ph[ip][jn];
		fv[7] = ph[ip][jp];
		fv[8] = ph[in][jp];

		fva[0]=ph[in][jn];
		fva[1]=ph[in][jp];
		fva[2]=ph[ip][jn];
		fva[3]=ph[ip][jp];

		fvb[0] = ph_nb(i, j);
		fvb[1] = ph_nb(in, j);
		fvb[2] = ph_nb(ip, j);
		fvb[3] = ph_nb(i, jn);
		fvb[4] = ph_nb(i, jp);

		ng = maxf();	

		if(ng){
			interP(&xmax,&ymax,&phmax);
			xpn[Npoint]=i*1.0+xmax;
			ypn[Npoint]=j*1.0+ymax;
			fprintf(fp, "%.8e  %.8e\n", xpn[Npoint], ypn[Npoint]);
//			if((j>=15&&j<=18&&k>=9&&k<=12)||(xpn[Npoint]>16.0&&xpn[Npoint]<17.0&&ypn[Npoint]>10.0&&ypn[Npoint]<11.0))
//				printf("(%d, %d) %f   %f\n", i, j, xpn[Npoint], ypn[Npoint]);
			Npoint++;
		}
	}
//	printf("!!!\n");

	fprintf(fp2, "%d,  %d,  %d,  ", Nx, Ny, Npoint);
	printf("\n Npoint = %d\n", Npoint);
	fprintf(fp0, "%d\t", Npoint);

	if ((Npoint>=100)&&(abs(Npoint-Npointold)==0))
	{
		if (switch_Npoint==1) {NstepCritical=step+50000;}
		switch_Npoint=0;
		fprintf(fp2,"END");
	}
	else switch_Npoint = 1;

	fclose(fp);
	fclose(fp2);


	if (Npoint<=12) return 0; else return 1;

	
}


void frprmn(double p[],int n,double ftol,int *iter,double *fret)
{
int i,j,its;
double gg,gam,dgg,fp,fp1,xi[NDIMM],xisum;
double g[NDIMM],h[NDIMM];
  
fp=func(p);

dfunc(p,xi);
/*
printf("\n %20.18e",fp);
for(i=0;i<NDIMM;i++){fp+=0.001*xi[i]*xi[i]; p[i]+=0.001*xi[i];}
printf("\n %20.18e \n %20.18e",fp,func(p));
for(i=0;i<NDIMM;i++)p[i]-=0.002*xi[i];
printf("\n %20.18e",func(p));
//exit(1);

xisum=0.0;
for(j=0;j<NDIMM;j++)xisum+=xi[j]*xi[j];
*/
for(j=0;j<n;j++)
    	{
      	g[j]=-xi[j];
    	xi[j]=h[j]=g[j];
    	}


for(its=1;its<=ITMAX;its++)
    	{
      	*iter=its;


        if(its%10==0)
		{
	  	//wtden(p);
	  	dfunc(p,xi);

	  	for(j=0;j<n;j++)
	    			{
	      			g[j]=-xi[j];
	      			xi[j]=h[j]=g[j];
	    			}
		}

      	linmin(p,xi,n,fret);

      	//printf("\n %e, %e, %e",func(p),fp,fabs(*fret-fp));
      	if((2.0*fabs(*fret-fp))<=(ftol*(fabs(*fret)+fabs(fp)+EPS)))
	//if(xisum<1.0e-4)
		{
	  	return;
		}
         
      	fp=func(p);
      	dfunc(p,xi);
        //xisum=0.0;
        //for(i=0;i<NDIMM;i++)xisum+=xi[i]*xi[i]; 
      	dgg=gg=0.0;

      	for(j=0;j<n;j++)
		{
	  	gg+=g[j]*g[j];
	  	dgg+=(xi[j]+g[j])*xi[j];
		}

      	if(gg==0.0) return;

      	gam=dgg/gg;

      	for(j=0;j<n;j++)
		{
	  	g[j]=-xi[j];
	  	xi[j]=h[j]=g[j]+gam*h[j];
		}
    	}
printf("Numerical Recipes run-time error.../n");
printf("Too many iterations in frprmn");
printf("...now exiting to system...\n");
}

void linmin(double *p,double *xi,int n,double *fret)
{
int j;
double xx,xmin,fx,fb,fa,bx,ax;

ncom=n;

for(j=0;j<n;j++)
    	{
      	pcom[j]=p[j];
      	xicom[j]=xi[j];
    	}

ax=0.0;
xx=1.0;

mnbrak(&ax,&xx,&bx,&fa,&fx,&fb);

*fret=brent(ax,xx,bx,TOL,&xmin);
 
//printf("\n xmin=%e,*fret=%e",xmin,*fret);

for(j=0;j<n;j++)
    	{
      	xi[j]*=xmin;
      	p[j]+=xi[j];
    	}
}

void mnbrak(double *ax,double *bx,double *cx,double *fa,double *fb,double *fc)
{
double ulim,u,r,q,fu,dum;

*fa=fldim(*ax);

*fb=fldim(*bx);

if(*fb>*fa)
    	{
      	SHFT(dum,*ax,*bx,dum)
      	SHFT(dum,*fb,*fa,dum)
    	}
*cx=(*bx)+GOLD*(*bx-*ax);
*fc=fldim(*cx);

while(*fb>*fc)
    	{
      	r=(*bx-*ax)*(*fb-*fc);
      	q=(*bx-*cx)*(*fb-*fa);
      	dum=fabs(q-r);

      	if(dum<TINY) dum=TINY;
      	if((q-r)<0.0) dum=-dum;

      	u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*dum);
      	ulim=(*bx)+GLIMIT*(*cx-*bx);

      	if(((*bx-u)*(u-*cx))>0.0)
		{
	  	fu=fldim(u);
	  	if(fu<*fc)
	    			{
	      			*ax=(*bx);
	      			*bx=u;
	      			*fa=(*fb);
	      			*fb=fu;

	      			return;
	    			}
	  	else if(fu>(*fb))
	    			{
	      			*cx=u;
	      			*fc=fu;

	      			return;
	    			}
	  	u=(*cx)+GOLD*(*cx-*bx);
	  	fu=fldim(u);
		}
      	else if(((*cx-u)*(u-ulim))>0.0)
		{
	  	fu=fldim(u);
	  	if(fu<(*fc))
	    			{
	      			SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
	      			SHFT(*fb,*fc,fu,fldim(u))
	    			}
		}
      	else if(((u-ulim)*(ulim-*cx))>=0.0)
		{
	  	u=ulim;
	  	fu=fldim(u);
		}
      	else 
		{
	     	u=(*cx)+GOLD*(*cx-*bx);
	     	fu=fldim(u);
	   	}
      	SHFT(*ax,*bx,*cx,u)
      	SHFT(*fa,*fb,*fc,fu)
    	}
}

double brent(double ax,double bx,double cx,double tol,double *xmin)
{
int iter;
double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm,stemp;
double e=0.0;

a=(ax<cx?ax:cx);
b=(ax>cx?ax:cx);
x=w=v=bx;
fw=fv=fx=fldim(x);

for(iter=1;iter<=ITMAX;iter++)
    	{
      	xm=0.5*(a+b);

      	tol2=2.0*(tol1=(tol*fabs(x))+ZEPS);
      	if(fabs(x-xm)<=(tol2-0.5*(b-a)))
		{
		//printf("\n %d",iter);
		//getchar();
	  	*xmin=x;
	  	return fx;
		}
      	if(fabs(e)>tol1)
		{
	  	r=(x-w)*(fx-fv);
	  	q=(x-v)*(fx-fw);
	  	p=(x-v)*q-(x-w)*r;
	  	q=2.0*(q-r);

	  	if(q>0.0) p=-p;

	  	q=fabs(q);
	  	etemp=e;
	  	e=d;

	  	if(fabs(p)>=(fabs(0.5*q*etemp))||p<=(q*(a-x))||p>=(q*(b-x)))
	    	                                    d=CGOLD*(e=(x>=xm?a-x:b-x));
	  	else
	    			{
	      			d=p/q;
	      			u=x+d;
	      			if((u-a)<tol2||(b-u)<tol2)
						{
		  				if((xm-x)>0.0) d=fabs(tol1);
		  				else d=-fabs(tol1);
		  				//d=SIGN(tol1,xm-x);
						}
	    			}
		}
      	else
		{
	  	d=CGOLD*(e=(x>xm?(a-x):(b-x)));
		}
      	if(d>0.0) stemp=fabs(tol1);
      	else stemp=-fabs(tol1);
      	//u=(fabs(d)>=tol1?(x+d):(x+SIGN(tol1,d)));
      	u=(fabs(d)>=tol1?(x+d):(x+stemp));
      	fu=fldim(u);
      	if(fu<=fx)
		{
	  	if(u>=x) a=x;
	  	else b=x;
	  	SHFT(v,w,x,u)
	  	SHFT(fv,fw,fx,fu)
		}
      	else
		{
	  	if(u<x) a=u;
	  	else b=u;
	  	if(fu<=fw||w==x)
	    			{
	      			v=w;
	     		 	w=u;
	      			fv=fw;
	      			fw=fu;
	    			}
	  	else if(fu<=fv||v==x||v==w)
	    			{
	     			v=u;
	      			fv=fu;
	    			}
		}
    	}
printf("\nToomany iteration in brent\n");

*xmin=x;
return fx;
}

double fldim(double x)
{
int j;
double f,xt[NDIMM];

for(j=0;j<ncom;j++)
    	{
      	xt[j]=pcom[j]+x*xicom[j];
    	}

f=func(xt);
return f;
}



int maxf()
{
	
	//if(f<0.50){
	if(fv[0]>0.550&&fv[0]>fv[1]&&fv[0]>fv[2]&&fv[0]>fv[3]&&fv[0]>fv[4]&&
		fv[0]>fva[0]&&fv[0]>fva[1]&&fv[0]>fva[2]&&fv[0]>fva[3])//&&
		//fvb[0]>fvb[1]&&fvb[0]>fvb[2]&&fvb[0]>fvb[3]&&fvb[0]>fvb[4])
		return 1;
	else return 0;
	//else return 0;}
	/*else 
		{
		if(fv[0]<0.50&&fv[0]<fv[1]&&fv[0]<fv[2]&&fv[0]<fv[3]&&
		fv[0]<fv[4]&&fv[0]<fva[0]&&fv[0]<fva[1]&&fv[0]<fva[2]&&fv[0]<fva[3]
		&&fv[0]<fvb[0]&&fv[0]<fvb[1]&&fv[0]<fvb[2]&&fv[0]<fvb[3])return 1;
		else return 0;
		}*/	
}


void interP(double *xmax,double *ymax,double *phmax)
	{
	/**** ax*(x-x0)^2+ay*(y-y0)^2+c *****/
	double p[5],dp[5];
	int Nstep2;
	xt[0]=0.0; 	yt[0]=0.0;
	xt[1]=-1.0;	yt[1]=0.0;
	xt[2]=1.0;	yt[2]=0.0;
	xt[3]=0.0;	yt[3]=-1.0;
	xt[4]=0.0;	yt[4]=1.0;
	xt[5]=-1.0;	yt[5]=-1.0;
	xt[6]=1.0;	yt[6]=-1.0;
	xt[7]=1.0;	yt[7]=1.0;
	xt[8]=-1.0;	yt[8]=1.0;

	p[0]=1.0;		//a=b
	p[1]=fv[0];		//c
	p[2]=0.0;		//x0
	p[3]=0.0;		//y0

	frprmn(p,4,TOL,&Nstep2,phmax);

	dfunc(p,dp);

	//printf("%.8e, %.8e, %.8e, %.8e, %.8e\n",dp[0],dp[1],dp[2],dp[3],dp[4],dp[5]);

	*xmax=p[2];
	*ymax=p[3];

	*phmax=p[1];	
	}

double func1(double a,double b,double c,double x0,double y0)
	{
	double fret,ft;
	int i;

	fret=0.0;
	for(i=0;i<9;i++)
		{
		ft=(fv[i]-a*(xt[i]-x0)*(xt[i]-x0)-b*(yt[i]-y0)*(yt[i]-y0)-c);
		fret=fret+ft*ft;
		}

	return fret;
	}

void dfunc(double x[],double dx[])
	{
	double ft0,ft,a,b,c,x0,y0;
	int i;
	
	a=x[0];	b=x[0];	c=x[1];	x0=x[2]; y0=x[3];

	for(i=0;i<4;i++)dx[i]=0.0;

	for(i=0;i<9;i++)
		{
		ft0=fv[i]-a*(xt[i]-x0)*(xt[i]-x0)-b*(yt[i]-y0)*(yt[i]-y0)-c;
		ft=-ft0*2.0*((xt[i]-x0)*(xt[i]-x0)+(yt[i]-y0)*(yt[i]-y0));

		dx[0]=dx[0]+ft;

		ft=-ft0*2.0;
		dx[1]=dx[1]+ft;

		ft=-ft0*2.0*2.0*a*(x0-xt[i]);
		dx[2]=dx[2]+ft;

		ft=-ft0*2.0*2.0*b*(y0-yt[i]);

		dx[3]=dx[3]+ft;
		
		}
	}

double func(double x[])
	{
	return func1(x[0],x[0],x[1],x[2],x[3]);
	}

double ph_nb(int j, int k)
{
	int jn, jp, kn, kp;
	double phav_nb;
	jn = j - 1;
	jp = j + 1;
	kn = k - 1;
	kp = k + 1;

	if(jn<0)jn+=Nx;
	else if(jp>=Nx)jp-=Nx;

	if(kn<0)kn+=Ny;
	else if(kp>=Ny)kp-=Ny;

	phav_nb = ph[j][k];
	phav_nb += ph[jn][k];
	phav_nb += ph[jp][k];
	phav_nb += ph[j][kn];
	phav_nb += ph[j][kp];

	phav_nb += (ph[jn][kn]+ph[jn][kp]+ph[jp][kn]+ph[jp][kp]);

	return phav_nb;
}








int caltheta(long step,FILE *fp0)
{
	int i, j, num;
	int Nx, Ny, Npoint;
	double e1, e2, dist, min_dist, max_dist, radiu, sum, theta, mean, nring, length[100000];
	char Fname[50];
	char fname[50];
	FILE *fp,*fpw;
	struct Pt *pt;

	//scanf("%lf",&nring); 
	nring = 15;

	fp = fopen("fname.txt", "w");
	fprintf(fp, "paraN%ld.dat", step);
	fclose(fp);

	fp = fopen("fname.txt", "r");
	fscanf(fp, "%s", fname);
	fclose(fp);

	fp=fopen(fname,"r");
	fscanf(fp,"%d, %d, %d",&Nx, &Ny, &Npoint);
	printf("%d , %d, %d\n",Nx, Ny, Npoint);
	fclose(fp);

	pt = (struct Pt *)malloc(sizeof(struct Pt)*Npoint);

	fp = fopen("fname.txt", "w");
	fprintf(fp, "max2d%ld.dat", step);
	fclose(fp);

	fp = fopen("fname.txt", "r");
	fscanf(fp, "%s", fname);
	fclose(fp);

	fp=fopen(fname,"r");


	for(i=0;i<Npoint;i++){
		fscanf(fp,"%lf%lf",&e1, &e2);
		pt[i].x = e1;
		pt[i].y = e2;
	}
	fclose(fp);

	min_dist = 1000000000.0;
	max_dist = 0.0;

	for(i=0;i<Npoint;i++)
	for(j=0;j<i;j++){
		if(i!=j){
			dist = sqrt(pow(pt[i].x - pt[j].x,2)+pow(pt[i].y - pt[j].y,2));
			if(dist < min_dist)
				min_dist = dist;
			if(dist > max_dist)
				max_dist = dist;
		}
	}

	printf("min_dist = %lf\n", min_dist);
	printf("max_dist = %lf\n", max_dist);


	fp = fopen("fname.txt", "w");
	fprintf(fp, "bond%ld.dat", step);
	fclose(fp);

	fp = fopen("fname.txt", "r");
	fscanf(fp, "%s", fname);
	fclose(fp);

	fp = fopen(fname,"w");
	num = 0;
	sum = 0.0;
	for(i=0;i<Npoint;i++){
		for(j=0;j<i;j++){
			dist = sqrt(pow(pt[i].x - pt[j].x,2)+pow(pt[i].y - pt[j].y,2));
			if(dist <= 1.4*18.4752){//min_dist){
				printf("num exist\n");
				sum += dist;
				length[num] = dist;
				num++;
				fprintf(fp,"%lf %lf %lf %lf\n", pt[i].x, pt[i].y, pt[j].x, pt[j].y);
			}
		}
	}
	fclose(fp);

	mean = sum * 1.0 / num;
	
	theta = 0.0;
	for(i=0;i<num;i++){
		theta += pow(length[i] - mean, 2);
	}
	
	theta /= (num- 1);
	theta = sqrt(theta);
	theta /= mean;			// theta stands for the relative error

	fp = fopen("cal_theta.dat","a");
	if(num>100000)
		fprintf(fp,"there are too many num, more than 100000!!!!!\n");
	else{
		fprintf(fp,"step = %ld\n",step);
		fprintf(fp,"%lf\t%lf\n", mean, theta);
		fprintf(fp,"min_dist = %lf\n", min_dist);
		fprintf(fp,"max_dist = %lf\n", max_dist);
		fprintf(fp,"Npoint = %d\n", Npoint);
		fprintf(fp,"Num_bond= %d\n", num);
		fprintf(fp,"mean = %lf\n", mean);
		fprintf(fp,"theta = %lf\n\n\n\n", theta);
		fprintf(fp0,"%.8e\t%.8e\t%.8e\t%.8e\t",theta,mean,max_dist,min_dist);
	}
	fclose(fp);	

	free(pt);

	if (theta<thetaCritical)
	return 0;
	else return 1;

	
}





