#ifdef SETTING_PARSER_H
#else
#define SETTING_PARSER_H 
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"
#include<vector>
#include<iostream>
using namespace std;
using namespace rapidjson;

class SettingParser{
	private:
		Document doc;
		char file[256];
	private:
		void parseFile(const char *s){
			FILE* fp;
			char readBuffer[256];
			fp=fopen(file,"r");
			if(fp==NULL){
				cerr<<"Error:File("<<file<<") doesn't exit."<<endl;
				return ;
			}
			FileReadStream is(fp,readBuffer,sizeof(readBuffer));
			doc.ParseStream(is);
			fclose(fp);
			assert(doc.IsObject());
		}
	public:
		SettingParser(){}
		SettingParser(const char *file){
			strcpy(this->file,file);
			parseFile(file);
		}
		SettingParser(SettingParser&p){
			strcpy(file,p.file);
			parseFile(file);
		}
		SettingParser&operator=(SettingParser&p){
			strcpy(file,p.file);
			parseFile(file);
		}
		void getAuthorInfo(const char *name,char *info){
			const Value&authorInfo=doc["AuthorInfo"];
			strcpy(info,authorInfo[name].GetString());
		}
		int getInt(const char *name)const{
			return doc[name].GetInt();
		}
		double getDouble(const char *name)const{
			return doc[name].GetDouble();
		}

		void getVector(const char *name,vector<double>&arr)const{
			const Value&a=doc[name];
			assert(a.IsArray());
			arr.resize(a.Size());
			for(SizeType i=0;i<a.Size();i++){
				arr[i]=a[i].GetDouble();
			}
		}
		void getVector(const char *name,vector<int>&arr)const{
			const Value&a=doc[name];
			assert(a.IsArray());
			arr.resize(a.Size());
			for(SizeType i=0;i<a.Size();i++){
				arr[i]=a[i].GetInt();
			}
		}
		void getBiVector(const char *name,vector<vector<double> >&arr)const{
			const Value&a=doc[name];
			assert(a.IsArray());
			arr.resize(a.Size());
			for(SizeType i=0;i<a.Size();i++){
				const Value&sa=a[i];
				assert(sa.IsArray());
				vector<double>&sarr=arr[i];
				sarr.resize(sa.Size());
				for(SizeType j=0;j<sa.Size();j++){
					sarr[j]=sa[j].GetDouble();
				}
			}
		}
		void getBiVector(const char *name,vector<vector<int> >&arr)const{
			const Value&a=doc[name];
			assert(a.IsArray());
			arr.resize(a.Size());
			for(SizeType i=0;i<a.Size();i++){
				const Value&sa=a[i];
				assert(sa.IsArray());
				vector<int>&sarr=arr[i];
				sarr.resize(sa.Size());
				for(SizeType j=0;j<sa.Size();j++){
					sarr[j]=sa[j].GetInt();
				}
			}
		}
		const char *getString(const char *name)const{
			return doc[name].GetString();
		}
		void getString(const char *name,char *dest)const {
			strcpy(dest,doc[name].GetString());
		}
		~SettingParser(){
		}
};
#endif
