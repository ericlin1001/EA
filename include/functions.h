#ifdef FUNCTION_H
#else
#define FUNCTION_H
#include "utils.h"
#include "cec14_test_func.h"
static double u(double x,double a,double k,double m){
	if(x>a)return k*pow(x-a,m);
	if(x<-a)return k*pow(-x-a,m);
	return 0;
}
class Function{
	private:
		char name[50];
		double xlow,xup;
		double fbest;
		bool isFindMin;
		int numDim;
		//
		int feCounter;
	private:
	public:
		virtual double evaluate(double *xs){
			feCounter++;
			return 0;
		}
		inline double evaluate(vector<double>&xs){
			return evaluate(&xs[0]);
		}
	protected:
		void setName(const char *s){
			strcpy(name,s);
		}
	public:
		Function(const char *s,double xlow,double xup,double fbest,int numDim){
			this->xlow=xlow;
			this->xup=xup;
			this->fbest=fbest;
			this->isFindMin=isFindMin;
			this->numDim=numDim;
			setName(s);
			feCounter=0;
		}
		double getFBest()const{return fbest;}
		int getFECounter()const{return feCounter;}
		int getNumDim()const{return numDim;}
		double getRange(int botOrUp){
			if(botOrUp==0)return xlow;
			return xup;
		}
		const char *getName()const{return name;}
};

#define DefFunction(name,xlow,xup,fbest,isFindMin) class name : public Function{\
	public: name(int numDim):Function(#name,xlow,xup,fbest,numDim){}\
			virtual double evaluate(double *xs){\
				Function::evaluate(xs);\
				int size=getNumDim();
#define EndDef }};	
DefFunction(PDEF3,-10,10,0,true)
	double res=0.0;
	for(int i=0;i<size-1;i++){
		double t=xs[i+1]-xs[i];
		double t1=xs[i]-1;
		res+=t*t*100.0+t1*t1;
	}
	return res;
EndDef

DefFunction(PDEF4,-500,500,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		res+=-xs[i]*sin(sqrt(fabs(xs[i])));
	}
	res+=(double)418.9829*(double)size;
	return res;
EndDef

DefFunction(F1,-100,100,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		double x=xs[i];
		res+=x*x;
	}
return res;
	EndDef

DefFunction(F2,-10,10,0,true)
	double res=0.0;
	double sum=0.0;
	double mul=1.0;
	for(int i=0;i<size;i++){
		double fabsx=fabs(xs[i]);
		sum+=fabsx;
		mul*=fabsx;
	}
res=sum+mul;
return res;
EndDef

DefFunction(F3,-100,100,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		double insum=0.0;
		for(int j=0;j<=i;j++){
			insum+=xs[j];
		}
		res+=insum*insum;
	}
return res;
EndDef

DefFunction(F4,-100,100,0,true)
	double res=fabs(xs[0]);
	for(int i=1;i<size;i++){
		double tmp=fabs(xs[i]);
		if(tmp<res)res=tmp;
	}
return res;
EndDef

//untest:
DefFunction(F5,-30,30,0,true)
	double res=0.0;
	for(int i=0;i<size-1;i++){
		double tmp=pow(xs[i+1]-xs[i]*xs[i],2)*100.0+pow(xs[i]-1.0,2);
		res+=tmp;
	}
return res;
EndDef

DefFunction(F6,-100,100,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		int tmp=floor(xs[i]+0.5);
		res+=tmp*tmp;
	}
return res;
EndDef

DefFunction(F7,-1.28,1.28,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		double tmp=pow(xs[i],4)*(double)(i+1);
		res+=tmp;
	}
	res+=drand();
return res;
EndDef

DefFunction(F8,-500,500,-12569.5,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		double tmp=-xs[i]*sin(sqrt(fabs(xs[i])));
		res+=tmp;
	}
return res;
EndDef

DefFunction(F9,-5.12,5.12,0,true)
	double res=0.0;
	for(int i=0;i<size;i++){
		double tmp=pow(xs[i],2)-(double)10.0*cos(xs[i]*2.0*M_PI)+10.0;
		res+=tmp;
	}
return res;
EndDef

DefFunction(F10,-32,32,0,true)
	double res=0.0;
	double sumx2=0.0;
	double sumcosx=0.0;
	for(int i=0;i<size;i++){
		sumx2+=pow(xs[i],2);
		sumcosx+=cos(xs[i]*M_PI*2.0);
	}
res=-20.0*exp(-0.2*sqrt(sumx2/(double)size))-exp(sumcosx/(double)size)+20.0+M_E;
return res;
EndDef

DefFunction(F11,-600.0,600.0,0,true)
	double res=0.0;
	double sumx2=0.0;
	double mulcos=1.0;
	for(int i=0;i<size;i++){
		sumx2+=pow(xs[i],2);
		mulcos*=cos(xs[i]/sqrt((double)i+1));
	}
res=sumx2/4000.0-mulcos+1.0;
return res;
EndDef

DefFunction(F12,-50,50,0,true)
	double res=0.0;
	double y1=1.0+(xs[0]+1.0)/4.0;
	double yd=1.0+(xs[size-1]+1.0)/4.0;
	double sumy=0.0;
	double sumu=0.0;
	//
	double yi,yi1;
	yi=y1;
	for(int i=0;i<size-1;i++){
		yi1=1.0+(xs[i+1]+1.0)/4.0;
		sumy+=pow(yi-1.0,2)*(1.0+10.0*pow(sin(M_PI*yi1),2));
		yi=yi1;
	}
for(int i=0;i<size;i++){
	sumu+=u(xs[i],10,100,4);
}
res=M_PI/(double)size*(10.0*pow(sin(M_PI*y1),2)+sumy+pow(yd-1,2))
	+sumu;
	return res;
	EndDef

DefFunction(F13,-50,50,0,true)
	double res=0.0;
	double sumx=0.0;
	double sumu=0.0;
	for(int i=0;i<size-1;i++){
		sumx+=pow(xs[i]-1,2)*(1+pow(sin(3.0*M_PI*xs[i+1]),2));
	}
for(int i=0;i<size;i++){
	sumu+=u(xs[i],5,100,4);
}
double xd=xs[size-1];
res=0.1*(pow(sin(3.0*M_PI*xs[0]),2)+sumx+
		pow(xd-1.0,2)*(1+pow(sin(2.0*M_PI*xd),2)))+sumu;
return res;
EndDef

#define DefCEC14Func(name,no) DefFunction(name,-100,100,no*100,true) \
	double res;\
	cec14_test_func(xs,&res,size,1,no);\
	return res;\
EndDef

DefCEC14Func(CEC14_F1,1);
DefCEC14Func(CEC14_F2,2);
DefCEC14Func(CEC14_F3,3);
DefCEC14Func(CEC14_F4,4);
DefCEC14Func(CEC14_F5,5);
DefCEC14Func(CEC14_F6,6);
DefCEC14Func(CEC14_F7,7);
DefCEC14Func(CEC14_F8,8);
DefCEC14Func(CEC14_F9,9);
DefCEC14Func(CEC14_F10,10);
DefCEC14Func(CEC14_F11,11);
DefCEC14Func(CEC14_F12,12);
DefCEC14Func(CEC14_F13,13);
DefCEC14Func(CEC14_F14,14);
DefCEC14Func(CEC14_F15,15);
DefCEC14Func(CEC14_F16,16);
DefCEC14Func(CEC14_F17,17);
DefCEC14Func(CEC14_F18,18);
DefCEC14Func(CEC14_F19,19);
DefCEC14Func(CEC14_F20,20);
DefCEC14Func(CEC14_F21,21);
DefCEC14Func(CEC14_F22,22);
DefCEC14Func(CEC14_F23,23);
DefCEC14Func(CEC14_F24,24);
DefCEC14Func(CEC14_F25,25);
DefCEC14Func(CEC14_F26,26);
DefCEC14Func(CEC14_F27,27);
DefCEC14Func(CEC14_F28,28);
DefCEC14Func(CEC14_F29,29);
DefCEC14Func(CEC14_F30,30);

class FunctionFactory{
	private:
	protected:
		vector<Function*>fs;
		FunctionFactory(){
		}
		void addFunction(Function*f){
			fs.push_back(f);
		}
	public:
		static FunctionFactory*instance;
		virtual void init(int numDim){
			cout<<"FunctionFactory init"<<endl;
			fs.resize(4);
			fs[0]=new F1(numDim);
			fs[1]=new F3(numDim);
			fs[2]=new PDEF3(numDim);
			fs[3]=new PDEF4(numDim);
		}
		static FunctionFactory *Instance(int numDim){
			if(instance==0){
				instance=new FunctionFactory();
				instance->init(numDim);
			}
			return instance;
		}
		Function*getFunction(int index)const{
			return fs[index];
		}
		int getNumFunction()const{
			return fs.size();
		}
		~FunctionFactory(){
			for(int i=0;i<getNumFunction();i++){
				delete fs[i];
			}
		}
};

FunctionFactory*FunctionFactory::instance=0;
class FunctionFactoryMy:public FunctionFactory{
	private:
		FunctionFactoryMy():FunctionFactory(){
		}
	public:
		static FunctionFactoryMy *Instance(int numDim){
			if(instance==0){
				instance=new FunctionFactoryMy();
				instance->init(numDim);
			}
			return (FunctionFactoryMy*)instance;
		}
		virtual void init(int numDim){
			addFunction(new CEC14_F1(numDim));
			addFunction(new CEC14_F2(numDim));
			addFunction(new CEC14_F3(numDim));
			addFunction(new CEC14_F4(numDim));
			addFunction(new CEC14_F5(numDim));
			addFunction(new CEC14_F6(numDim));
			addFunction(new CEC14_F7(numDim));
			addFunction(new CEC14_F8(numDim));
			addFunction(new CEC14_F9(numDim));
			addFunction(new CEC14_F10(numDim));
			addFunction(new CEC14_F11(numDim));
			addFunction(new CEC14_F12(numDim));
			addFunction(new CEC14_F13(numDim));
			addFunction(new CEC14_F14(numDim));
			addFunction(new CEC14_F15(numDim));
			addFunction(new CEC14_F16(numDim));
			addFunction(new CEC14_F17(numDim));
			addFunction(new CEC14_F18(numDim));
			addFunction(new CEC14_F19(numDim));
			addFunction(new CEC14_F20(numDim));
			addFunction(new CEC14_F21(numDim));
			addFunction(new CEC14_F22(numDim));
			addFunction(new CEC14_F23(numDim));
			addFunction(new CEC14_F24(numDim));
			addFunction(new CEC14_F25(numDim));
			addFunction(new CEC14_F26(numDim));
			addFunction(new CEC14_F27(numDim));
			addFunction(new CEC14_F28(numDim));
			addFunction(new CEC14_F29(numDim));
			addFunction(new CEC14_F30(numDim));
		}
};
#endif
