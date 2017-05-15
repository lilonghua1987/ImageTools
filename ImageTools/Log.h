#pragma once

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "Matrix.h"

#ifdef _WIN32
#include <direct.h>
#include <io.h>
#elif _LINUX
#include <stdarg.h>
#include <sys/stat.h>
#endif

#ifdef _WIN32
#define ACCESS _access
#define MKDIR(a) _mkdir((a))
#elif _LINUX
#define ACCESS access
#define MKDIR(a) mkdir((a),0755)
#endif

#pragma warning(disable: 4250)

//using FileModel = ios_base::openmode;

typedef ios_base::openmode FileModel;

enum FileModelType
{
	FileMOdelIn = ios_base::in,
	FileModelOut = ios_base::out,
	FileModelAte = ios_base::ate,
	FileModelApp = ios_base::app,
	FileModelTrunc = ios_base::trunc,
	FileModelBinary = ios_base::binary,
	FileModelCreate = ios_base::_Nocreate,
	FileModelRePlace = ios_base::_Noreplace
};

class FileStorage:public fstream
{
public:
	FileStorage():fstream(),xmlFlag(false){};

	FileStorage(const char fileName[],FileModel model = FileModelOut)
		:fstream(fileName,model)
		,xmlFlag(false)
	{
		assert(fileName != nullptr);
		fileExt(fileName);
		if (xmlFlag)
		{
			*this << xmlHead << "\n";
			*this << xmlRootNode << "\n";
		}
	};

	void open(const char fileName[],FileModel model = FileModelOut)
	{
		if (!this->is_open() && !this->good())
		{
			this->fstream::open(fileName,model);
			fileExt(fileName);
			if (xmlFlag)
			{
				*this << xmlHead << "\n";
				*this << xmlRootNode << "\n";
			}
		}
	}

	/*template<class T>
	friend FileStorage& operator << (FileStorage& fs, const T& t)
	{
		fs.fstream::operator<<(t);
		return fs;
	}*/

	~FileStorage();

private:
	void fileExt(const char fileName[])
	{
		if (fileName != nullptr)
		{
			std::string name = std::string(fileName);
			name = name.substr(name.find_last_of('.'));
			if(!name.empty() && name == ".xml") xmlFlag = true;
		}
	}

private:
	bool xmlFlag;
	const static char xmlHead[];
	const static char xmlRootNode[];
	const static char xmlRootNodeEnd[];

};

class Log
{
public:
	Log(void);
	~Log(void);

	Log(const char *filename);

	void systemLog(const char format[],double out);
	void systemLog(const char format[],float out);
	//void systemLog(const char format[],LONG64 out);
	void systemLog(const char format[],int out);
	void systemLog(const char format[],char out);
	void systemLog(const char format[]);

	template<typename T>
	void systemLog(const Matrix<T>& m);

	template<typename T>
	void systemLog(T& out);

private:
	int CreateDirectory(const char* pDir);

private:
	const char* PATH;
	const char* fileName;
	const char* logPath;
};


template<typename T>
void Log::systemLog(const Matrix<T>& m)
{
	string path(PATH);path.append(fileName);
	fstream file(path,ios::app);
	if(!file){
		return;
	}
	file.seekg(0,ios::end);   //让文件指针定位到文件末尾 

	for (int i = 0; i < m.row; i++)
	{
		for (int j = 0; j < m.column; j++)
		{
			for (int c = 0; c < m.channel; c++)
			{
				file<<m.at(i,j,c)<<" ";
			}
		}
		file<<"\n";
	}
	file.close();
}


template<typename T>
void Log::systemLog(T& out)
{
	string path(PATH);path.append(fileName);
	fstream file(path,ios::app);
	if(!file){
		return;
	}
	file.seekg(0,ios::end);   //让文件指针定位到文件末尾 
	file<<out;
	file.close();
}
