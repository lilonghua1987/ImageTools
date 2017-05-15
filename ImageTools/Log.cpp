#include "Log.h"


//FileStorage
const char FileStorage::xmlHead[] = "<?xml version=\"1.0\"?>";
const char FileStorage::xmlRootNode[] = "<llh_storage>";
const char FileStorage::xmlRootNodeEnd[] = "</llh_storage>";

FileStorage::~FileStorage()
{
	if(xmlFlag) *this << xmlRootNodeEnd;
	//this->fstream::~basic_fstream();
	if (this->is_open()) this->close();
}


Log::Log(void)
	:PATH("log\\")
	,fileName("log.txt")
	,logPath(NULL)
{
	CreateDirectory(PATH);
	string path(PATH);path.append(fileName);
	logPath = path.c_str();

	//打开一个文件时，如果文件不存在，创建该文件
	//file.open("log.txt",ios::app);
}


Log::~Log(void)
{
}


Log::Log(const char* fileName)
	:PATH("log\\")
	,fileName(NULL)
	,logPath(NULL)
{
	CreateDirectory(PATH);
	Log::fileName = fileName;
	string path(PATH);path.append(Log::fileName);
	logPath = path.c_str();	
}


int Log::CreateDirectory(const char* pDir)
{
	int i = 0;
	int iRet;
	int iLen;
	char* pszDir;

	if(NULL == pDir)
	{
		return 0;
	}
	
	pszDir = _strdup(pDir);
	iLen = strlen(pszDir);

	// 创建中间目录
	for (i = 0;i < iLen;i ++)
	{
		if (pszDir[i] == '\\' || pszDir[i] == '/')
		{ 
			pszDir[i] = '\0';

			//如果不存在,创建
			iRet = ACCESS(pszDir,0);
			if (iRet != 0)
			{
				iRet = MKDIR(pszDir);
				if (iRet != 0)
				{
					return -1;
				} 
			}
			//支持linux,将所有\换成/
			pszDir[i] = '/';
		} 
	}

	iRet = MKDIR(pszDir);
	free(pszDir);
	return iRet;
}


void Log::systemLog(const char format[],double out)
{
	string path(PATH);path.append(fileName);
	fstream file(path,ios::app);
	if(!file)
	{
		return;
	}
	file.seekg(0,ios::end);   //让文件指针定位到文件末尾 
	file<<format<<out;
	file.close();
}


void Log::systemLog(const char format[],float out)
{
	string path(PATH);path.append(fileName);
	fstream file(path,ios::app);
	if(!file)
	{
		return;
	}
	file.seekg(0,ios::end);   //让文件指针定位到文件末尾 
	file<<format<<out;
	file.close();
}


//void Log::systemLog(const char format[],LONG64 out){
//	fstream file("temp\\log.txt",ios::app);
//	if(!file){
//		return;
//	}
//	file.seekg(0,ios::end);   //让文件指针定位到文件末尾 
//	file<<format<<","<<out;
//	file.close();
//}


void Log::systemLog(const char format[],int out)
{
	string path(PATH);path.append(fileName);
	fstream file(path,ios::app);
	if(!file)
	{
		return;
	}
	file.seekg(0,ios::end);   //让文件指针定位到文件末尾 
	file<<format<<out;
	file.close();
}

void Log::systemLog(const char format[],char out)
{
	string path(PATH);path.append(fileName);
	fstream file(path,ios::app);
	if(!file)
	{
		return;
	}
	file.seekg(0,ios::end);   //让文件指针定位到文件末尾 
	file<<format<<out;
	file.close();
}

void Log::systemLog(const char format[])
{
	string path(PATH);path.append(fileName);
	fstream file(path,ios::app);
	if(!file)
	{
		return;
	}
	file.seekg(0,ios::end);   //让文件指针定位到文件末尾 
	file<<format;
	file.close();
}

//void Log::systemLog(const Image& img)
//{
//	string path(PATH);path.append(fileName);
//	fstream file(path,ios::app);
//	if(!file){
//		return;
//	}
//	for (int i = 0; i < img.height; i++)
//	{
//		for (int j = 0; j < img.width; j++)
//		{
//			Pixel<uchar> p = img.getPixel(i,j);
//			if (img.channel == 1)
//			{
//				file<<p.red;
//			}else
//			{
//				file<<"["<<p.red<<","<<p.green<<","<<p.blue<<"]";
//			}
//			file<<" ";
//		}
//		file<<"\n";
//	}
//	file.close();
//}