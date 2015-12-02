#include<iostream>
#include<fstream>
#include "include/template.h"
using namespace std;

#ifdef USE_MPI
#define OMPI_IMPORTS
#include "mpi.h"
#endif


class FunctionReverse:public Function{
	Function *ff;
	public:
		FunctionReverse(Function *f):Function(f->getName(),f->getRange(0),f->getRange(1),-f->getFBest(),f->getNumDim()){
			ff=f;
			char buff[100];
			strcpy(buff,"-");
			strcat(buff,f->getName());
			setName(buff);
		}
		virtual double evaluate(double *xs){
			return -(ff->evaluate(xs));
		}
};
class EA{
	char name[100];
	protected:
	void setName(const char *n){
		strcpy(name,n);
	}
	public:
		//input:f,MaxFEs
		//output:x,fx
		virtual void getMin(Function *f,int MaxFEs,vector<double>&x,double &fx)=0;
		void getMax(Function *f,int MaxFEs,vector<double>&x,double &fx){
			Function rf=FunctionReverse(f);
			getMin(&rf,MaxFEs,x,fx);
			fx=-fx;
		}
		virtual const char *getName()const{
			return name;
		}
};
typedef SettingParser SearchParam ;
class Save{
	private:
			ofstream out;
			int ix;
	public:
		Save(const char *title,const char *x,const char *y){
			char buff[150];
			sprintf(buff,"%s.figdata",title,x,y);
			out.open(buff);
			out<<x<<"\t"<<y<<endl;
			ix=0;
		}
		void add(int x,double y){
			out<<x<<"\t"<<y<<endl;
		}
		void add(double y){
			out<<ix<<"\t"<<y<<endl;
			ix+=1;
		}
		~Save(){
		}
};
class DE:public EA{
	private:
		//about function:f
		Function *f;
		double xmin;
		double xmax;
		int numDim;
		//algorithm related parameters.
		int PopSize;
		double F,CR;
		//
		vector<vector<double> >x;//x,trail x.
		vector<vector<double> >tmpX;
		vector<double>fx;
		vector<double>tmpFx;
		vector<double>tx;
		//
		int bestI;
		//
	private:
		void updateX(){
			//main process
			ASSERT(PopSize>=3);
			RandomPermutation perm(PopSize);
			for(int i=0;i<PopSize;i++){
				perm.generate();
				int a=bestI; int b=perm.next(); int c=perm.next();
				int randDim=rand()%numDim;
				for(int j=0;j<numDim;j++){
					if(j==randDim||drand()<CR){
						tx[j]=x[a][j]+F*(x[b][j]-x[c][j]);
						if(tx[j]<xmin || tx[j]>xmax){
							tx[j]=drand(xmin,xmax);
						}
					}else{
						tx[j]=x[i][j];
					}
				}
				double ftx=f->evaluate(&tx[0]);
				if(ftx<fx[i]){
					x[i]=tx;
					fx[i]=ftx;
					if(ftx<fx[bestI]){
						bestI=i;
					}
				}
			}
		}
	private:
		double getBestFx()const{
			return fx[bestI];
		}
		void update(int maxGeneration){
			Save s(f->getName(),"Generation","F*/F");
			s.add(f->getFBest()/getBestFx());
			for(int g=1;g<=maxGeneration;g++){
				updateX();
				s.add(f->getFBest()/getBestFx());
			}
		}
	public:
		void initParam(SearchParam *param){
			PopSize=param->getInt("PopSize");
			F=param->getDouble("F");
			CR=param->getDouble("CR");
			setName(param->getString("Name"));
		}
		virtual void getMin(Function *f,int MaxFEs,vector<double>&out_x,double &out_fx){
			this->f=f;
			numDim=f->getNumDim();
			xmin=f->getRange(0);
			xmax=f->getRange(1);
			//population initializing....
			//allocate space.
			tx.resize(numDim);
			x.resize(PopSize);
			fx.resize(PopSize);
			for(int i=0;i<PopSize;i++){
				x[i].resize(numDim);
			}
			bestI=0;
			for(int i=0;i<PopSize;i++){
				for(int d=0;d<numDim;d++){
					x[i][d]=drand(xmin,xmax);
				} 
				fx[i]=f->evaluate(&x[i][0]);
				if(fx[i]<fx[bestI]){
					bestI=i;
				}
			}
			//update, main process.
			update(MaxFEs/PopSize-1);
			out_x=x[bestI];
			out_fx=fx[bestI];
		}
};

/*
#if ALGORITHM==4
int PDE(int processId,int PopSizerocess,Function*f,vector<double>&bestX,double &bestF){
	srand(time(NULL));//srand in every process.
	//MPI:
	const int TAG=99;
	const int TAG2=98;
	MPI_Status status;
	//DE:
	const int MaxFE=300000;
	const int NumAlgorithm=min(PopSizerocess-1,4);
	int PopSize=50;
	int numDim=f->getNumDim();
	//
	DE de;
	if(processId!=0){
		de.init(processId-1,PopSize);
		de.begin(f);
	}
	ASSERT(NumAlgorithm>0);
	const int MaxGeneration=MaxFE/(NumAlgorithm*PopSize);
	ASSERT(MaxGeneration>0);
	//master
	vector<double>meanFs;
	vector<vector<int> > migrateMap;
	//slave:
	vector<int>migrateVec;
	vector<double>XF;
	vector<double>bestXF;
	XF.resize(numDim+1);

	vector<int>popSizes;
	if(processId==0){
		//for master process
		meanFs.resize(NumAlgorithm);
		migrateMap.resize(NumAlgorithm);
		popSizes.resize(NumAlgorithm);
		for(int i=0;i<NumAlgorithm;i++){
			migrateMap[i].resize(NumAlgorithm);
			popSizes[i]=PopSize;
		}
	}else{
		migrateVec.resize(NumAlgorithm);
	}
	const int MinPopSize=5;
	for(int g=1;g<=MaxGeneration;g++){
		if(processId==0){
			//master process:
#ifdef DEBUG
			if(g==1)cout<<"Master:Generation("<<g<<") begin"<<endl;
#endif
			double Pmigrate=0.01+0.99*(exp((double)10.0*g/(double)MaxGeneration)-1.0)/(exp(10.0)-1.0);
			for(int i=1;i<=NumAlgorithm;i++){
				MPI_Recv(&meanFs[i-1],1,MPI_DOUBLE,i,TAG,MPI_COMM_WORLD,&status);
			}
			//build migrateMap.
			for(int i=0;i<NumAlgorithm;i++){
				migrateMap[i][i]=0;
				for(int j=i+1;j<NumAlgorithm;j++){
					if(drand()<Pmigrate){
						if(meanFs[i]<meanFs[j]){
							if(popSizes[j]<=MinPopSize){
								migrateMap[i][j]=0;
							}else{
							//i is better than j.
								migrateMap[i][j]=-1;//recve from j, j->i
								popSizes[i]++;
								popSizes[j]--;
							}
						}else{
							if(popSizes[i]<=MinPopSize){
								migrateMap[i][j]=0;
							}else{
							//j is better than i
								migrateMap[i][j]=1;//send to j. i->i
								popSizes[i]--;
								popSizes[j]++;
							}
						}
					}else{
							migrateMap[i][j]=0;//do nothing.
					}
					migrateMap[j][i]=-migrateMap[i][j];
				}
			}
#ifdef DEBUG
			cout<<"MigrateMap:"<<endl;
			for(int i=0;i<NumAlgorithm;i++){
				for(int j=0;j<NumAlgorithm;j++){
					cout<<migrateMap[i][j]<<",";
				}
				cout<<endl;
			}
			cout<<endl;
			if(g!=MaxGeneration)cout<<"Master:Generation("<<g+1<<") begin"<<endl;
#endif
			for(int i=1;i<=NumAlgorithm;i++){
				MPI_Send(&migrateMap[i-1][0],migrateMap[0].size(),MPI_INT,i,TAG,MPI_COMM_WORLD);
			}
		}else{
			//slave processes:
			de.update(1);
			//cout<<"Process("<<processId<<"):"<<"end de.update()"<<endl;
			double meanF=de.getMeanF();
#ifdef DEBUG
			cout<<"Process("<<processId<<"):"<<"meanF:"<<meanF<<endl;
#endif
			MPI_Send(&meanF,1,MPI_DOUBLE,0,TAG,MPI_COMM_WORLD);
		
			MPI_Recv(&migrateVec[0],migrateVec.size(),MPI_INT,0,TAG,MPI_COMM_WORLD,&status);
			//cout<<"Process("<<processId<<"):"<<"end recv migrateMap"<<endl;
			//Attention: without explicit synchronization, there might be a deadlock or error.
			for(int i=0;i<NumAlgorithm;i++){
				if(migrateVec[i]==1){
					//1:send
					ASSERT(de.getPopSize()>1);
					int randI=rand()%de.getPopSize();
					XF=de.del(randI);
#ifdef DEBUG
			cout<<"Process("<<processId<<"):"<<"sending to P("<<i+1<<")"<<endl;
#endif
					MPI_Send(&XF[0],XF.size(),MPI_DOUBLE,i+1,TAG,MPI_COMM_WORLD);
				}else if(migrateVec[i]==-1){
					//-1:recv
					MPI_Recv(&XF[0],XF.size(),MPI_DOUBLE,i+1,TAG,MPI_COMM_WORLD,&status);
					de.add(XF);
				}
			}
		}
	}
	//return only in  process0.
	if(processId==0){
		for(int i=1;i<=NumAlgorithm;i++){
			MPI_Recv(&XF[0],XF.size(),MPI_DOUBLE,i,TAG2,MPI_COMM_WORLD,&status);
			if(i==1){
				bestXF=XF;
			}else{
				if(XF.back()<bestXF.back()){
					bestXF=XF;
				}
			}
		}
		bestF=bestXF.back();
		bestX=bestXF;
		bestX.pop_back();
		return 0;
	}else{
		de.getOutput(XF,bestF);
		XF.push_back(bestF);
		MPI_Send(&XF[0],XF.size(),MPI_DOUBLE,0,TAG2,MPI_COMM_WORLD);
	}
	return -1;
}
vector<double> runPDE(int id,int idSize,Function*f,int maxRun){
	vector<double>results;
	results.resize(maxRun);
	for(int run=0;run<maxRun;run++){
		vector<double>bestX;
		double bestF;
		PDE(id,idSize,f,bestX,bestF);
		results[run]=fabs(bestF-(f->getBest()));
	}
	return results;
}
#endif
vector<double> runSerialDE(DE &de,Function*f,int maxRun){
	vector<double>results;
	const int MaxFE=300000;
	int numDim=f->getNumDim();
	vector<double>bestX;
	double bestF;
	results.resize(maxRun);
	for(int run=0;run<maxRun;run++){
		de.solve(f,MaxFE/de.getPopSize(),bestX,bestF);
		results[run]=fabs(bestF-(f->getBest()));
	}
	return results;
}

*/
int unused_main(int argc,char *argv[]){
//int main(int argc,char *argv[]){
	SearchParam param("DE1.json");
	Trace(param.getDouble("CR"));
	Trace(param.getInt("MaxFEs"));
}
//int unused_main(int argc,char *argv[]){
int main(int argc,char *argv[]){
	srand(time(NULL));
	//parse the setting.
	SearchParam param("DE1.json");
	FunctionFactoryMy &funGenerator=FunctionFactoryMy::Instance(param.getInt("NumDim"));
	//const int numTestFunction=funGenerator.getNumFunction();
	const int numTestFunction=1;
	DE de;
	de.initParam(&param);

	//cout<<"Runing DE("<<de.getName()<<") "<<maxRun<<"times."<<endl;
	//printf("F\tmean\tstd\n");
	Tic::tic("begin");
	vector<double>x;
	double fx;
	for(int i=0;i<numTestFunction;i++){
		Function*f=funGenerator.getFunction(i);
		de.getMin(f,param.getInt("MaxFEs"),x,fx);
		printf("%s(%g,%g)",f->getName(),fx,f->getFBest());
		//vector<double>results=runSerialDE(de,f,maxRun);
		//double min,max,mean,std;
		//calStatistics(results,min,max,mean,std);
		//printf("%s\t%g\t%g\n",f->getShortName(),mean,std);
	}
	Tic::tic("end");
	return 0;
}
