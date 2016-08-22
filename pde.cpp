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
			sprintf(buff,"%s.figdata",title);
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
#define SaveData s.add(f->getFBest()/getBestFx());
			Save s(f->getName(),"Generation","F*/F");
				SaveData;
			for(int g=1;g<=maxGeneration;g++){
				updateX();
				SaveData;
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
int main(int argc,char *argv[]){
	srand(time(NULL));
	//parse the setting.
	SearchParam param("DE1.json");
	FunctionFactoryMy *funGenerator=FunctionFactoryMy::Instance(param.getInt("NumDim"));
	const int numTestFunction=funGenerator->getNumFunction();
	DE de;
	de.initParam(&param);
	cout<<"Runing "<<de.getName()<<" "<<endl;
	Tic::tic("begin");
	vector<double>x;
	double fx;
	cout<<"FunName(MyBestF,Optima)"<<endl;
	for(int i=0;i<numTestFunction;i++){
		Function*f=funGenerator->getFunction(i);
		de.getMin(f,param.getInt("MaxFEs"),x,fx);
		printf("%s(%g,%g)",f->getName(),fx,f->getFBest());
		Tic::tic("end");
	}
	Tic::tic("end");
	return 0;
}
