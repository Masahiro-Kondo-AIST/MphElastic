//================================================================================================//
//------------------------------------------------------------------------------------------------//
//    MPH-Elastic : Moving Particle Hydrodynamics for Elastic body                                      //
//------------------------------------------------------------------------------------------------//
//    Developed by    : Masahiro Kondo                                                            //
//    Distributed in  : 2023                                                                      //
//    Lisence         : GPLv3                                                                     //
//    For instruction : see README                                                                //
//    For theory      : see the following references                                              //
//     [1] JSCES Paper No.20070031,  https://doi.org/10.11421/jsces.2007.20070031                 //
//     [2] Int.J.Numer.Meth.Engng.81 (2010) 1514-1528,  https://doi.org/10.1002/nme.2744          //
//    Copyright (c) 2010  Masahiro Kondo                                                          //
//    Copyright (c) 2023  National Institute of Advanced Industrial Science and Technology (AIST) //
//================================================================================================//


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <vector>

#include "vec3T.hpp"
#include "log.h"
#include "errorfunc.h"

using namespace std;

#define DIM 3

#define NORMAL    0
#define PERIODIC  1
#define SYMMETRIC 2

#define FORCE_VOLUME 1
#define DISPLACEMENT 2

#define DEFAULT_LOG  "sample.log"
#define DEFAULT_DATA "sample.data"
#define DEFAULT_GRID "sample.grid"
#define DEFAULT_PROF "sample%03d.prof"
#define DEFAULT_VTK  "sample%03d.vtk"
#define DEFAULT_ENE  "sample.ene"

// Simulation Time
static double EndTime;
static double Time;
static double Dt;

// Prof Output
static double OutputInterval=0.0;
static double OutputNext=0.0;
static double VtkOutputInterval=0.0;
static double VtkOutputNext=0.0;

// Simulation Domain
typedef struct{
    double pos;
    int prop;
}domain_t;
static domain_t DomainMin[DIM] = {{0.0,NORMAL},{0.0,NORMAL},{0.0,NORMAL}};
static domain_t DomainMax[DIM] = {{0.0,NORMAL},{0.0,NORMAL},{0.0,NORMAL}};


typedef struct{
    double density;
    double young;
    double poisson;
    double mu;
    double lambda;
    double artificial;
}property_t;

static int ParticlePropertyCount;
static vector<property_t> ParticleProperty;



// Particle
static int ParticleCount;
static double AverageParticleSpacing;
static vector<int> Property;
static vector< vec3<double> > Position;
static vector< vec3<double> > Velocity;
static vector< mat3<double> > Normalizer;
static vector< mat3<double> > DeformGradient;
static vector< mat3<double> > Strain;
static vector< mat3<double> > Stress;
static vector< vec3<double> > Acceleration;
static vector< vec3<double> > ExternalForce;
static vector< vec3<double> > InitialPosition;
static vector< vector<int> >  InitialNeighbor;
static vector< double >       Volume;
static vector< double >       Radius;

static int TotalParticleCount;
static int VirtualParticleCount;
static vector<int> Original;
static vector< mat3<double> > Conversion;

// Energy
static double KineticEnergy=0.0;
static double StrainPotential=0.0;
static double ArtificialPotential =0.0;
static double ExternalWork=0.0;
static double PreviousPotential=0.0;

// Equilibrium
static int Equilibrium=0;
static int ResetVelocityCount=0;

// BackGroundCells
#define CellId(iCX,iCY,iCZ)  ((iCX)*CellCount[1]*CellCount[2]+(iCY)*CellCount[2]+(iCZ))
static int CellCount[DIM];
static double MaxRadius = 0.0;
static double CellWidth = 0.0;
static vector<int> CellParticleCount;
static vector<int> CellParticleBegin;
static vector<int> CellParticle;


static void readDataFile(char *filename);
static void readGridFile(char *filename);
static void writeProfFile(char *filename);
static void writeVtkFile(char *filename);
static void writeEnergyFile(char *filename);

static void initializeExternalForce();
static void initializeDomain();
static void generateBoundary();
static void calculateNeighbor( vector< vector<int> > &neighbor, const vector< vec3<double> > &position );
static void calculateMuLambda();
static void calculateNormalizer();
static void resetAcceleration();
static void calculateStressForce();
static void calculateArtificialForce();
static void calculateConvection();
static double weight(const vec3<double> r, const double radius);


	clock_t cFrom, cTill, cStart, cEnd;
	clock_t cInitialize=0, cReset=0, cStress=0, cArtificial=0, cOutput=0, cEnergyOutput=0, cConvection=0;


int main(int argc, char* argv[])
{
	char logfilename[1024]  = DEFAULT_LOG;
	char datafilename[1024] = DEFAULT_DATA;
	char gridfilename[1024] = DEFAULT_GRID;
	char proffilename[1024] = DEFAULT_PROF;
	char vtkfilename[1024]  = DEFAULT_VTK;
	char energyfilename[1024]=DEFAULT_ENE;
	{
		if(argc>1)strcpy(datafilename,argv[1]);
		if(argc>2)strcpy(gridfilename,argv[2]);
		if(argc>3)strcpy(proffilename,argv[3]);
		if(argc>4)strcpy( vtkfilename,argv[4]);
		if(argc>5)strcpy(energyfilename,argv[5]);
		if(argc>6)strcpy( logfilename,argv[6]);
	}
	log_open(logfilename);
	{
		time_t t=time(NULL);
		log_printf("start reading files at %s\n",ctime(&t));
	}
	readDataFile(datafilename);
	readGridFile(gridfilename);
	
	
	
	{
		time_t t=time(NULL);
		log_printf("start initialization at %s\n",ctime(&t));
	}
	initializeDomain();
	generateBoundary();
	calculateNeighbor( InitialNeighbor, InitialPosition );
	calculateMuLambda();
	calculateNormalizer();

	
	{
		time_t t=time(NULL);
		log_printf("start main roop at %s\n",ctime(&t));
	}
	int iStep=(int)(Time/Dt);
	OutputNext = Time;
	VtkOutputNext = Time;
	cStart = clock();
	cFrom = cStart;
	
	cTill = clock(); cInitialize += (cTill-cFrom); cFrom = cTill;
	while (Time < EndTime + 1.0e-5*Dt){
		resetAcceleration();
		cTill = clock(); cReset += (cTill-cFrom); cFrom = cTill;
		calculateStressForce();
		cTill = clock(); cStress += (cTill-cFrom); cFrom = cTill;
		calculateArtificialForce();
		cTill = clock(); cArtificial += (cTill-cFrom); cFrom = cTill;
		if( Time + 1.0e-5*Dt >= OutputNext ){
			char filename[256];
			sprintf(filename,proffilename,iStep);
			writeProfFile(filename);
			log_printf("@ Prof Output Time:%lf\n",Time);
			OutputNext += OutputInterval;
		}
		if( Time + 1.0e-5*Dt >= VtkOutputNext ){
			char filename[256];
			sprintf(filename,vtkfilename,iStep);
			writeVtkFile(filename);
			log_printf("@ Vtk Output Time:%lf\n",Time);
			VtkOutputNext += VtkOutputInterval;
		}
		cTill = clock(); cOutput += (cTill-cFrom); cFrom = cTill;
		writeEnergyFile(energyfilename);
		cTill = clock(); cEnergyOutput += (cTill-cFrom); cFrom = cTill;
// 
		calculateConvection();
		cTill = clock(); cConvection += (cTill-cFrom); cFrom = cTill;
		
		Time += Dt;
		iStep++;
	}
	
	{
		time_t t=time(NULL);
		log_printf("end main roop at %s\n",ctime(&t));
	}
	
	cEnd = cTill;
    {
        time_t t=time(NULL);
    	log_printf("initialize:              %lf [CPU sec]\n", (double)cInitialize/CLOCKS_PER_SEC);
    	log_printf("reset acceleration:      %lf [CPU sec]\n", (double)cReset/CLOCKS_PER_SEC);
    	log_printf("stress force:            %lf [CPU sec]\n", (double)cStress/CLOCKS_PER_SEC);
    	log_printf("artificial force:        %lf [CPU sec]\n", (double)cArtificial/CLOCKS_PER_SEC);
    	log_printf("output:                  %lf [CPU sec]\n", (double)cOutput/CLOCKS_PER_SEC);
    	log_printf("energy output:           %lf [CPU sec]\n", (double)cEnergyOutput/CLOCKS_PER_SEC);
    	log_printf("convection:              %lf [CPU sec]\n", (double)cConvection/CLOCKS_PER_SEC);
    	log_printf("total:                   %lf [CPU sec]\n", (double)(cInitialize+cReset+cStress+cArtificial+cOutput+cEnergyOutput+cConvection)/CLOCKS_PER_SEC);
    	log_printf("total (check):           %lf [CPU sec]\n", (double)(cEnd-cStart)/CLOCKS_PER_SEC);
    }
	return 0;
	
}

static void readDataFile(char *filename)
{
	FILE * fp;
	char buf[1024];
	char command[1024];
	char sval[1024];
	// int ival;
	double dval;
	const int reading_global=0;
	const int reading_particle_property=1;
	int mode=reading_global;
	int iProperty;
	
	
	fp=fopen(filename,"r");
	mode=reading_global;
	while(fp!=NULL && !feof(fp) && !ferror(fp)){
		fgets(buf,sizeof(buf),fp);
		if(buf[0]=='#'){}
		else if(sscanf(buf," Dt %lf",&Dt)==1){mode=reading_global;}
		else if(sscanf(buf," OutputInterval %lf",&OutputInterval)==1){mode=reading_global;}
		else if(sscanf(buf," VtkOutputInterval %lf",&VtkOutputInterval)==1){mode=reading_global;}
		else if(sscanf(buf," EndTime %lf",&EndTime)==1){mode=reading_global;}
		else if(sscanf(buf," MinX %s",command)==1){
			if(command[0]=='#'){}
			else if(strcmp(command,"symmetric")==0){
				DomainMin[0].prop=SYMMETRIC;
			}
			else if(strcmp(command,"periodic")==0){
				DomainMin[0].prop=PERIODIC;
			}
			mode=reading_global;
		}
		else if(sscanf(buf," MinY %s",command)==1){
			if(command[0]=='#'){}
			else if(strcmp(command,"symmetric")==0){
				DomainMin[1].prop=SYMMETRIC;
			}
			else if(strcmp(command,"periodic")==0){
				DomainMin[1].prop=PERIODIC;
			}
			mode=reading_global;
		}
		else if(sscanf(buf," MinZ %s",command)==1){
			if(command[0]=='#'){}
			else if(strcmp(command,"symmetric")==0){
				DomainMin[2].prop=SYMMETRIC;
			}
			else if(strcmp(command,"periodic")==0){
				DomainMin[2].prop=PERIODIC;
			}
			mode=reading_global;
		}
		else if(sscanf(buf," MaxX %s",command)==1){
			if(command[0]=='#'){}
			else if(strcmp(command,"symmetric")==0){
				DomainMax[0].prop=SYMMETRIC;
			}
			else if(strcmp(command,"periodic")==0){
				DomainMax[0].prop=PERIODIC;
			}
			mode=reading_global;
		}
		else if(sscanf(buf," MaxY %s",command)==1){
			if(command[0]=='#'){}
			else if(strcmp(command,"symmetric")==0){
				DomainMax[1].prop=SYMMETRIC;
			}
			else if(strcmp(command,"periodic")==0){
				DomainMax[1].prop=PERIODIC;
			}
			mode=reading_global;
		}
		else if(sscanf(buf," MaxZ %s",command)==1){
			if(command[0]=='#'){}
			else if(strcmp(command,"symmetric")==0){
				DomainMax[2].prop=SYMMETRIC;
			}
			else if(strcmp(command,"periodic")==0){
				DomainMax[2].prop=PERIODIC;
			}
			mode=reading_global;
		}
		else if(sscanf(buf," ParticlePropertyCount %d", &ParticlePropertyCount)==1){
			ParticleProperty.resize(ParticlePropertyCount);
			for(int iProp=0;iProp<ParticlePropertyCount;++iProp){
				ParticleProperty[iProp].density=0.0;
				ParticleProperty[iProp].young=0.0;
				ParticleProperty[iProp].poisson=0.0;
				ParticleProperty[iProp].mu=0.0;
				ParticleProperty[iProp].lambda=0.0;
				ParticleProperty[iProp].artificial=0.0;
			}
			mode=reading_particle_property;
		}
		else if(sscanf(buf," ParticlePropertyId %d", &iProperty)==1 && mode==reading_particle_property){
			if(iProperty<0 || ParticlePropertyCount<=iProperty){
				iProperty=-1;
			}
		}
		else if(sscanf(buf," density %lf", &dval)==1 && mode==reading_particle_property){
			if(0<=iProperty && iProperty<ParticlePropertyCount ){
				ParticleProperty[iProperty].density = dval;
			}
		}
		else if(sscanf(buf," young %lf", &dval)==1 && mode==reading_particle_property){
			if(0<=iProperty && iProperty<ParticlePropertyCount ){
				ParticleProperty[iProperty].young = dval;
			}
		}
		else if(sscanf(buf," poisson %lf", &dval)==1 && mode==reading_particle_property){
			if(0<=iProperty && iProperty<ParticlePropertyCount ){
				ParticleProperty[iProperty].poisson = dval;
			}
		}
		else if(sscanf(buf," artificial %lf", &dval)==1 && mode==reading_particle_property){
			if(0<=iProperty && iProperty<ParticlePropertyCount ){
				ParticleProperty[iProperty].artificial = dval;
			}
		}
		else if(sscanf(buf," Equilibrium %s",command)==1){
			if(strcmp(command,"on")==0){
				Equilibrium=1;
			}
			mode=reading_global;
		}
		else{
			log_printf("Invalid line in data file \"%s\"\n", buf);
		}
	}
	fclose(fp);
	return;
}


static void readGridFile(char *filename)
{
    FILE *fp;
    char buf[1024];

 
    fp=fopen(filename,"r");
    
    fgets(buf,sizeof(buf),fp);
    sscanf(buf,"%lf",&Time);
    fgets(buf,sizeof(buf),fp);
    sscanf(buf,"%d %lf  %lf %lf %lf  %lf %lf %lf",
           &ParticleCount,
           &AverageParticleSpacing,
           &DomainMin[0].pos, &DomainMax[0].pos,
           &DomainMin[1].pos, &DomainMax[1].pos,
           &DomainMin[2].pos, &DomainMax[2].pos);

    Property.resize(ParticleCount);
    Position.resize(ParticleCount);
    Velocity.resize(ParticleCount);
    Normalizer.resize(ParticleCount);
    DeformGradient.resize(ParticleCount);
    Strain.resize(ParticleCount);
    Stress.resize(ParticleCount);
    Acceleration.resize(ParticleCount);
    ExternalForce.resize(ParticleCount);
    InitialPosition.resize(ParticleCount);
    InitialNeighbor.resize(ParticleCount);
    Volume.resize(ParticleCount);
    Radius.resize(ParticleCount);
    Original.resize(ParticleCount);
    Conversion.resize(ParticleCount);

    for(int iP=0;iP<ParticleCount;++iP){
        fgets(buf,sizeof(buf),fp);
        sscanf(buf,"%d  %lf %lf %lf  %lf %lf %lf  %lf %lf %lf  %lf %lf",
               &Property[iP],
               &Position[iP][0],&Position[iP][1],&Position[iP][2],
               &Velocity[iP][0],&Velocity[iP][1],&Velocity[iP][2],
               &InitialPosition[iP][0],&InitialPosition[iP][1],&InitialPosition[iP][2],
               &Volume[iP], &Radius[iP]
               );
    }
    fclose(fp);
    return;
    
}

static void writeProfFile(char *filename)
{
	double pressure=0.0;
	mat3<double> deviatoric;
	double vonMises;
	FILE *fp;
	fp=fopen(filename,"w");
	
	fprintf(fp,"%e\n",Time);
	fprintf(fp,"%d %e  %e %e %e %e %e %e\n",
		ParticleCount, 
		AverageParticleSpacing, 
		DomainMin[0].pos, DomainMax[0].pos,
		DomainMin[1].pos, DomainMax[1].pos,
		DomainMin[2].pos, DomainMax[2].pos
	);
	
	for(int iP=0;iP<ParticleCount;++iP){
		fprintf(fp,"%d  %e %e %e  %e %e %e  %e %e %e  %e %e  %e %e\n",
			Property[iP],
			Position[iP][0],Position[iP][1],Position[iP][2],
			Velocity[iP][0],Velocity[iP][1],Velocity[iP][2],
			InitialPosition[iP][0],InitialPosition[iP][1],InitialPosition[iP][2],
			Volume[iP], Radius[iP]);
	}
	fclose(fp);
}


static void writeVtkFile(char *filename)
{
		
	FILE *fp=fopen(filename, "w");
	
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "Unstructured Grid Example\n");
	fprintf(fp, "ASCII\n");
	
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp, "POINTS %d float\n", ParticleCount);
	for(int iP=0;iP<ParticleCount;++iP){
		fprintf(fp, "%e %e %e\n", (float)Position[iP][0], (float)Position[iP][1], (float)Position[iP][2]);
	}
	fprintf(fp, "CELLS %d %d\n", ParticleCount, 2*ParticleCount);
	for(int iP=0;iP<ParticleCount;++iP){
		fprintf(fp, "1 %d ",iP);
	}
	fprintf(fp, "\n");
	fprintf(fp, "CELL_TYPES %d\n", ParticleCount);
	for(int iP=0;iP<ParticleCount;++iP){
		fprintf(fp, "1 ");
	}
	fprintf(fp, "\n");
	
	fprintf(fp, "\n");
	
	fprintf(fp, "POINT_DATA %d\n", ParticleCount);
	fprintf(fp, "SCALARS label float 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for(int iP=0;iP<ParticleCount;++iP){
		fprintf(fp, "%d\n", Property[iP]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "SCALARS Volume float 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for(int iP=0;iP<ParticleCount;++iP){
		fprintf(fp, "%e\n",(float)Volume[iP]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "SCALARS Radius float 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for(int iP=0;iP<ParticleCount;++iP){
		fprintf(fp, "%e\n", (float)Radius[iP]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "SCALARS Pressure float 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for(int iP=0;iP<ParticleCount;++iP){
		const double pressure = -(Stress[iP][0][0]+Stress[iP][1][1]+Stress[iP][2][2])/DIM;
		fprintf(fp, "%e\n", (float)pressure);
	}
	fprintf(fp, "\n");
	for(int iD=0;iD<DIM;++iD){
		for(int jD=0;jD<DIM;++jD){
			fprintf(fp, "\n");    fprintf(fp, "SCALARS stress[%d][%d] float 1\n",iD,jD);
			fprintf(fp, "LOOKUP_TABLE default\n");
			for(int iP=0;iP<ParticleCount;++iP){
				fprintf(fp, "%e\n", (float)Stress[iP][iD][jD]);
			}
		}
	}
	fprintf(fp, "VECTORS velocity float\n");
	for(int iP=0;iP<ParticleCount;++iP){
		fprintf(fp, "%e %e %e\n", (float)Velocity[iP][0], (float)Velocity[iP][1], (float)Velocity[iP][2]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "VECTORS displacement float\n");
	for(int iP=0;iP<ParticleCount;++iP){
		const double displacement[DIM]={Position[iP][0]-InitialPosition[iP][0],Position[iP][1]-InitialPosition[iP][1],Position[iP][2]-InitialPosition[iP][2]};
		fprintf(fp, "%e %e %e\n", (float)displacement[0], (float)displacement[1], (float)displacement[2]);
	}
	fprintf(fp, "\n");
	
	fflush(fp);
	fclose(fp);
}


static void writeEnergyFile(char *filename)
{
    static FILE *fp;
    static int init_flag=0;
    if(init_flag==0){
        fp=fopen(filename,"w");
        init_flag=1;
    }
    KineticEnergy=0.0;
    for(int iP=0;iP<ParticleCount;++iP){
        KineticEnergy += 0.5*ParticleProperty[Property[iP]].density*Volume[iP] * (Velocity[iP]*Velocity[iP]);
    }
    fprintf(fp,"%e %e %e %e %e %e\n",
            Time,
            KineticEnergy,
            StrainPotential,
            ArtificialPotential,
            ExternalWork,
            KineticEnergy + StrainPotential + ArtificialPotential - ExternalWork
            );
    fflush(fp);
}


static void initializeDomain()
{
    
    CellWidth = AverageParticleSpacing;
    MaxRadius = 0.0;
    for(int iP=0;iP<ParticleCount;++iP){
        if(MaxRadius<Radius[iP]){
            MaxRadius=Radius[iP];
        }
    }
    
    int range = (int)(MaxRadius/CellWidth)+1;
    CellCount[0] = (int)((DomainMax[0].pos - DomainMin[0].pos)/CellWidth)+1 +2*range;
    CellCount[1] = (int)((DomainMax[1].pos - DomainMin[1].pos)/CellWidth)+1 +2*range;
    CellCount[2] = (int)((DomainMax[2].pos - DomainMin[2].pos)/CellWidth)+1 +2*range;

    CellParticleCount.resize( CellCount[0]*CellCount[1]*CellCount[2] );
    CellParticleBegin.resize( CellCount[0]*CellCount[1]*CellCount[2] );
}

void generateBoundary( void ){

    double maxradius;
    for(int iP=0;iP<ParticleCount;++iP){
        if(maxradius<Radius[iP]){
            maxradius=Radius[iP];
        }
    }
    const mat3<double> unit = mat3<double>(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
    for(int iP=0;iP<ParticleCount;++iP){
        Original[iP] = iP;
        Conversion[iP] = unit;
    }

    VirtualParticleCount = 0;
    TotalParticleCount = ParticleCount;
    for(int iDim=0;iDim<DIM;++iDim){
        if(DomainMin[iDim].prop==NORMAL){
        }
        else if(DomainMin[iDim].prop==SYMMETRIC){
            for(int iP=0;iP<TotalParticleCount;++iP){
                if(DomainMin[iDim].pos < InitialPosition[iP][iDim] && InitialPosition[iP][iDim] < maxradius+DomainMin[iDim].pos){
                    int original = Original[iP];
                    vec3<double> position = InitialPosition[iP];
                    position[iDim] = 2.0*DomainMin[iDim].pos-position[iDim];
                    mat3<double> conversion = Conversion[iP];
                    conversion[iDim][iDim] *= -1;
                    Original.push_back(original);
                    InitialPosition.push_back(position);
                    Conversion.push_back(conversion);
                    VirtualParticleCount++;
                }
            }
        }
        else if(DomainMin[iDim].prop==PERIODIC){
            for(int iP=0;iP<TotalParticleCount;++iP){
                if(DomainMin[iDim].pos < InitialPosition[iP][iDim] && InitialPosition[iP][iDim] < maxradius+DomainMin[iDim].pos){
                    int original = Original[iP];
                    vec3<double> position = InitialPosition[iP];
                    position[iDim] += DomainMax[iDim].pos - DomainMin[iDim].pos;
                    mat3<double> conversion = Conversion[iP];
                    Original.push_back(original);
                    InitialPosition.push_back(position);
                    Conversion.push_back(conversion);
                    VirtualParticleCount++;
                }
            }
        }

        if(DomainMax[iDim].prop==NORMAL){
        }
        else if(DomainMax[iDim].prop==SYMMETRIC){
            for(int iP=0;iP<TotalParticleCount;++iP){
                if(DomainMax[iDim].pos-maxradius < InitialPosition[iP][iDim] && InitialPosition[iP][iDim] < DomainMax[iDim].pos){
                    int original = Original[iP];
                    vec3<double> position = InitialPosition[iP];
                    position[iDim] = 2.0*DomainMax[iDim].pos-position[iDim];
                    mat3<double> conversion = Conversion[iP];
                    conversion[iDim][iDim] *= -1;
                    Original.push_back(original);
                    InitialPosition.push_back(position);
                    Conversion.push_back(conversion);
                    VirtualParticleCount++;
                }
            }
        }
        else if(DomainMax[iDim].prop==PERIODIC){
            for(int iP=0;iP<TotalParticleCount;++iP){
                if(DomainMax[iDim].pos-maxradius < InitialPosition[iP][iDim] && InitialPosition[iP][iDim] < DomainMax[iDim].pos){
                    int original = Original[iP];
                    vec3<double> position = InitialPosition[iP];
                    position[iDim] -= DomainMax[iDim].pos-DomainMin[iDim].pos;
                    mat3<double> conversion = Conversion[iP];
                    Original.push_back(original);
                    InitialPosition.push_back(position);
                    Conversion.push_back(conversion);
                    VirtualParticleCount++;
                }
            }
        }
        
        TotalParticleCount = ParticleCount + VirtualParticleCount;
    }
}



static void calculateNeighbor( vector< vector<int> > &neighbor, const vector< vec3<double> > &position )
{
 
    int range = (int)(MaxRadius/CellWidth) + 1;
    
    // store to cells
    CellParticle.resize( TotalParticleCount );
    for(int iC=0;iC<CellCount[0]*CellCount[1]*CellCount[2];++iC){
        CellParticleCount[iC]=0;
    }
    for(int iP=0; iP<TotalParticleCount; ++iP){
        const int iCX=(int)((position[iP][0]-DomainMin[0].pos)/CellWidth)+range;
        const int iCY=(int)((position[iP][1]-DomainMin[1].pos)/CellWidth)+range;
        const int iCZ=(int)((position[iP][2]-DomainMin[2].pos)/CellWidth)+range;
        const int iC=CellId(iCX,iCY,iCZ);
        CellParticleCount[iC]++;
    }
    int iCellParticle = 0;
    for(int iC=0;iC<CellCount[0]*CellCount[1]*CellCount[2];++iC){
        CellParticleBegin[iC] = iCellParticle;
        iCellParticle+=CellParticleCount[iC];
    }
    for(int iC=0;iC<CellCount[0]*CellCount[1]*CellCount[2];++iC){
        CellParticleCount[iC]=0;
    }
    for(int iP=0; iP<TotalParticleCount; ++iP){
        const int iCX=(int)((position[iP][0]-DomainMin[0].pos)/CellWidth)+range;
        const int iCY=(int)((position[iP][1]-DomainMin[1].pos)/CellWidth)+range;
        const int iCZ=(int)((position[iP][2]-DomainMin[2].pos)/CellWidth)+range;
        const int iC=CellId(iCX,iCY,iCZ);
        CellParticle[ CellParticleBegin[iC] + CellParticleCount[iC] ] = iP;
        CellParticleCount[iC]++;
    }
    
    // calculate neighbor
    for(int iP=0;iP<ParticleCount;++iP){
        neighbor[iP].resize(0);
    }  
    for(int iP=0;iP<ParticleCount;++iP){
        const int iCX=(int)((position[iP][0]-DomainMin[0].pos)/CellWidth)+range;
        const int iCY=(int)((position[iP][1]-DomainMin[1].pos)/CellWidth)+range;
        const int iCZ=(int)((position[iP][2]-DomainMin[2].pos)/CellWidth)+range;
        
        for(int jCX=iCX-range;jCX<=iCX+range;++jCX){
            for(int jCY=iCY-range;jCY<=iCY+range;++jCY){
                for(int jCZ=iCZ-range;jCZ<=iCZ+range;++jCZ){
                    const int jC=CellId(jCX,jCY,jCZ);
                    for(int jCP=CellParticleBegin[jC];jCP<CellParticleBegin[jC]+CellParticleCount[jC];++jCP){
                        int jP=CellParticle[jCP];
                        if(iP!=jP){
                            const vec3<double> rij = 2.0*Radius[iP]/(Radius[iP]+Radius[jP])*(position[jP]-position[iP]);
                            if(rij*rij<=Radius[iP]*Radius[iP]){
                                neighbor[iP].push_back(jP);
                            }
                        }
                    }
                }
            }
        }
    }
}


static double weight(const vec3<double> rij, const double radius)
{
    return (radius*radius)/(rij*rij) - 1.0;
}

static void calculateMuLambda()
{
    for(int iProp=0;iProp<ParticlePropertyCount;++iProp){
        const double E = ParticleProperty[iProp].young;
        const double v = ParticleProperty[iProp].poisson;
        ParticleProperty[iProp].lambda = E*v/((1.0+v)*(1-2.0*v));
        ParticleProperty[iProp].mu     = 0.5*E/(1+v);
    }
}

static void calculateNormalizer()
{
    
    for(int iP=0;iP<ParticleCount;++iP){
        Normalizer[iP]=mat3<double>(0);
        const int cN = (int)InitialNeighbor[iP].size();
        for(int iN=0;iN<cN;++iN){
            const int jP = InitialNeighbor[iP][iN];
            const double cij = 2.0*Radius[iP]/(Radius[iP]+Radius[jP]);
            const vec3<double> rij0 = (InitialPosition[jP] - InitialPosition[iP]);
            const double wij = weight(cij*rij0, Radius[iP]);
            Normalizer[iP] += wij * (rij0 % rij0);
        }
        Normalizer[iP] = Normalizer[iP].inverse();
    }
}

static void resetAcceleration()
{
    for(int iP=0;iP<ParticleCount;++iP){
        Acceleration[iP]=vec3<double>(0.0,0.0,0.0);
    }
}
    
static void calculateStressForce()
{
    StrainPotential = 0.0;
    
    for(int iP=0;iP<ParticleCount;++iP){
        const property_t &propI = ParticleProperty[Property[iP]];
        DeformGradient[iP]=mat3<double>(0);
        const int cN = (int)InitialNeighbor[iP].size();
        for(int iN=0;iN<cN;++iN){
            const int jP = InitialNeighbor[iP][iN];
            const int jOP = Original[jP];
            const double cij = 2.0*Radius[iP]/(Radius[iP]+Radius[jP]);
            const vec3<double> ui   = Position[iP] - InitialPosition[iP];
            const vec3<double> uj   = Conversion[jP]*(Position[jOP]-InitialPosition[jOP]);
            const vec3<double> rij0  = (InitialPosition[jP] - InitialPosition[iP]);
            const vec3<double> rij = (rij0 + (uj-ui));
            const double wij = weight(cij*rij0, Radius[iP]);
            DeformGradient[iP] += wij * (rij % rij0);
        }

        DeformGradient[iP] = DeformGradient[iP] * Normalizer[iP];
        Strain[iP] = 0.5*( (DeformGradient[iP].trans()*DeformGradient[iP]) - mat3<double>(1.0) );
        Stress[iP] = 2.0*propI.mu*Strain[iP] + propI.lambda * Strain[iP].trace() * mat3<double>(1.0);
        StrainPotential +=0.5*(
                               Stress[iP][0][0]*Strain[iP][0][0]
                               +Stress[iP][0][1]*Strain[iP][0][1]
                               +Stress[iP][0][2]*Strain[iP][0][2]
                               +Stress[iP][1][0]*Strain[iP][1][0]
                               +Stress[iP][1][1]*Strain[iP][1][1]
                               +Stress[iP][1][2]*Strain[iP][1][2]
                               +Stress[iP][2][0]*Strain[iP][2][0]
                               +Stress[iP][2][1]*Strain[iP][2][1]
                               +Stress[iP][2][2]*Strain[iP][2][2]
                               )*Volume[iP];
    }

    for(int iP=0;iP<ParticleCount;++iP){
        const property_t &propI = ParticleProperty[Property[iP]];
        const mat3<double> matrix = DeformGradient[iP] * Stress[iP] * Normalizer[iP];
        const int cN = (int)InitialNeighbor[iP].size();
        for(int iN=0;iN<cN;++iN){
            const int jP = InitialNeighbor[iP][iN];
            const int jOP = Original[jP];
            const double cij = 2.0*Radius[iP]/(Radius[iP]+Radius[jP]);
            const property_t &propJ = ParticleProperty[Property[jOP]];
            const vec3<double> ui   = Position[iP] - InitialPosition[iP];
            const vec3<double> uj   = Conversion[jP]*(Position[jOP]-InitialPosition[jOP]);
            const vec3<double> rij0  = (InitialPosition[jP] - InitialPosition[iP]);
            const vec3<double> rij = (rij0 + (uj-ui));
            const double wij = weight(cij*rij0, Radius[iP]);
            const vec3<double> stressforce = wij * (matrix * rij0);

            Acceleration[iP] += 1.0/propI.density * stressforce;
            Acceleration[jOP]-= 1.0/propJ.density * Conversion[jP].inverse()*stressforce * Volume[iP]/Volume[jOP];
        }
    }
    
}

static void calculateArtificialForce()
{
    ArtificialPotential = 0.0;
    for(int iP=0;iP<ParticleCount;++iP){
        const property_t &propI = ParticleProperty[Property[iP]];
        const int cN = InitialNeighbor[iP].size();
          for(int iN=0;iN<cN;++iN){
            const int jP = InitialNeighbor[iP][iN];
            const int jOP = Original[jP];
            const double cij = 2.0*Radius[iP]/(Radius[iP]+Radius[jP]);
            const property_t &propJ = ParticleProperty[Property[jOP]];
            const vec3<double> ui   = Position[iP] - InitialPosition[iP];
            const vec3<double> uj   = Conversion[jP]*(Position[jOP]-InitialPosition[jOP]);
            const vec3<double> rij0  = (InitialPosition[jP] - InitialPosition[iP]);
            const vec3<double> rij = (rij0 + (uj-ui));
            const double wij = weight(cij*rij0, Radius[iP]);
            const double young = propI.artificial;
            const double sum_rijrijwij = 8.0/15.0 * M_PI * Radius[iP]*Radius[iP]*Radius[iP]*Radius[iP]*Radius[iP]/Volume[iP];
            const vec3<double> artificialforce
                = DIM / sum_rijrijwij * young * (rij - DeformGradient[iP]*rij0) * wij;
            Acceleration[iP] += 1.0/propI.density * artificialforce;
            Acceleration[jOP] -= 1.0/propJ.density * Conversion[jP].inverse()*artificialforce * Volume[iP]/Volume[jOP];

            ArtificialPotential += 0.5 * artificialforce * (rij - DeformGradient[iP]*rij0) * Volume[iP];
        }
    }
}


static void calculateConvection()
{
    // Equilibrium
    if( Equilibrium!=0){
        static int init_flag=0;
        if(init_flag==0){
            PreviousPotential = StrainPotential + ArtificialPotential - ExternalWork;
            init_flag=1;
        }
        else{
            if( PreviousPotential > StrainPotential + ArtificialPotential - ExternalWork ){
                PreviousPotential = StrainPotential + ArtificialPotential - ExternalWork;
            }
            else{
                for(int iP=0;iP<ParticleCount;++iP){
                    Position[iP] -= Velocity[iP]*Dt;
                }
                Time -= Dt;
                for(int iP=0;iP<ParticleCount;++iP){
                    Velocity[iP] = vec3<double>(0.0,0.0,0.0);
                }
                ResetVelocityCount++;
                log_printf("Time:%e, %d th Reset Velocity\n",Time,ResetVelocityCount);
            }
        }
    }
    for(int iP=0;iP<ParticleCount;++iP){
        Velocity[iP] += Acceleration[iP]*Dt;
        Position[iP] += Velocity[iP]*Dt;
        ExternalWork += ExternalForce[iP]*Velocity[iP]*Dt;
    }
}



