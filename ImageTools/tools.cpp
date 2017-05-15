#include "tools.h"

using namespace std;


void RunTimer::start()
{	
	m_begin=clock();
}


float RunTimer::stop()
{ 
	m_end=clock();
	return ( float(m_end-m_begin)/CLOCKS_PER_SEC );
}


void RunTimer::timeDisplay(std::string disp)
{ 
	std::cout << " Running time <" << disp << ">: " << stop() << " seconds." << std::endl;
}



void RunTimer::fpsDisplay(std::string disp)
{ 
	std::cout << " Running time <" << disp << ">: " << 1.0f/stop() << " frame per second." << std::endl;
}


//
tools::tools(void)
{
}


tools::~tools(void)
{
}


string tools::fileNameFromTime(string prefix,string suffix)
{
	string fileName;
	fileName.append(prefix);

	time_t	nowtime = time(NULL);
	struct	tm	*p;
	p = gmtime(&nowtime);
	char	filename[256] = {0};
	sprintf(filename,"%d_%d_%d_%d_%d_%d",1900+p->tm_year,1+p->tm_mon,p->tm_mday,p->tm_hour+8,p->tm_min,p->tm_sec);
	fileName.append(filename);
	fileName.append(suffix);
	return fileName;
}


string tools::fileNameFromTime(const char* prefix,const char* suffix)
{
	string fileName;
	fileName.append(prefix);

	time_t	nowtime = time(NULL);
	struct	tm	*p;
	p = gmtime(&nowtime);
	char	filename[256] = {0};
	sprintf(filename,"%d_%d_%d_%d_%d_%d",1900+p->tm_year,1+p->tm_mon,p->tm_mday,p->tm_hour+8,p->tm_min,p->tm_sec);
	fileName.append(filename);
	fileName.append(suffix);
	return fileName;
}


string tools::fileNameFromTime(string path,string name,string suffix)
{
	string fileName;
	if (!path.empty())
	{
		if ((path.find_last_of("\\") == path.length()-1) || (path.find_last_of("//") == path.length()-1)|| (path.find_last_of("/") == path.length()-1))
		{
			fileName.append(path);
		}else
		{
			fileName.append(path);
			fileName.append("\\");
		}		
	}

	if (!name.empty())
	{
		fileName.append(name);
	}

	if(suffix.empty()) return string();

	time_t	nowtime = time(NULL);
	struct	tm	*p;
	p = gmtime(&nowtime);
	char	filename[256] = {0};
	sprintf(filename,"%d_%d_%d_%d_%d_%d",1900+p->tm_year,1+p->tm_mon,p->tm_mday,p->tm_hour+8,p->tm_min,p->tm_sec);
	fileName.append(filename);
	fileName.append(suffix);
	return fileName;
}


string tools::fileNameFromTime(const char* path,const char* name,const char* suffix)
{
	string fileName;
	if (path)
	{
		int iLen = strlen(path);		
		if (path[iLen-1] == '\\' || path[iLen-1] == '/')
		{ 
			fileName.append(path);
		}else
		{
			fileName.append(path);
			fileName.append("\\");
		}
	}

	if (name)
	{
		fileName.append(name);
	}

	time_t	nowtime = time(NULL);
	struct	tm	*p;
	p = gmtime(&nowtime);
	char	filename[256] = {0};
	sprintf(filename,"%d_%d_%d_%d_%d_%d",1900+p->tm_year,1+p->tm_mon,p->tm_mday,p->tm_hour+8,p->tm_min,p->tm_sec);
	fileName.append(filename);
	fileName.append(suffix);
	return fileName;
}


string tools::fileName(string prefix,long w,long h,string suffix)
{
	string fileName;
	fileName.append(prefix);
	fileName.append(to_string((long long )w));
	fileName.append("_");
	fileName.append(to_string((long long )h));
	fileName.append(suffix);
	return fileName;
}


string tools::fileName(const char* prefix,long w,long h,const char* suffix)
{
	string fileName;
	fileName.append(prefix);
	fileName.append(to_string((long long )w));
	fileName.append("_");
	fileName.append(to_string((long long )h));
	fileName.append(suffix);
	return fileName;
}


string tools::fileName(string path,string name,long w,long h,string suffix)
{
	string fileName;
	if (!path.empty())
	{
		if ((path.find_last_of("\\") == path.length()-1) || (path.find_last_of("//") == path.length()-1)|| (path.find_last_of("/") == path.length()-1))
		{
			fileName.append(path);
		}else
		{
			fileName.append(path);
			fileName.append("\\");
		}		
	}

	if (!name.empty())
	{
		fileName.append(name);
	}
	fileName.append(to_string((long long )w));
	fileName.append("_");
	fileName.append(to_string((long long )h));
	fileName.append(suffix);
	return fileName;
}


string tools::fileName(const char* path,const char* name,long w,long h,const char* suffix)
{
	string fileName;
	if (path)
	{
		int iLen = strlen(path);		
		if (path[iLen-1] == '\\' || path[iLen-1] == '/')
		{ 
			fileName.append(path);
		}else
		{
			fileName.append(path);
			fileName.append("\\");
		}
	}

	if (name)
	{
		fileName.append(name);
	}

	fileName.append(to_string((long long )w));
	fileName.append("_");
	fileName.append(to_string((long long )h));
	fileName.append(suffix);
	return fileName;
}


void tools::getWords(std::vector<std::string>& words, const std::string& str, const char splite)
{
	string word;
	for (int i = 0; i < str.length(); i++)
	{
		char temp = str.at(i);
		if ((i == (str.length() - 1)) && (temp != splite))
		{
			word += temp;
			temp = splite;
		}
		if (temp != splite)
		{
			word += temp;
			continue;
		}
		else
		{
			words.push_back(word);
			word.swap(string());
		}
	}
}


string tools::getFileName(const string name)
{
	if(name.empty()) return string();
	int pos = name.find_last_of("\\");
	if((pos+1) == name.length()) return string();

	return name.substr(pos+1,name.length());
}


int tools::CreatDir(const char * const pDir)
{
	int i = 0;
	int iRet = 0;
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


void  tools::listFileByDir(string dir, string fileExtension, vector<string>* fileList)
{
	string fileFolder = dir,searchDir = dir; 

	if (!dir.empty())
	{
		/*int i = dir.find_last_of("\\"),j = dir.find_last_of("//");*/
	
		if ((dir.find_last_of("\\") == dir.length()-1) || (dir.find_last_of("//") == dir.length()-1)|| (dir.find_last_of("/") == dir.length()-1))
		{
			searchDir.append("*");
		}else
		{
			searchDir.append("\\*");
			fileFolder.append("\\");
		}
	}else
	{
		searchDir.append("\\*");
	}
  
  
    // 遍历文件夹  
  
    struct _finddata_t fileInfo;    // 文件信息结构体  
  
    // 1. 第一次查找  
    long findResult = _findfirst(searchDir.c_str(), &fileInfo);            
    if (findResult == -1)  
    {  
        _findclose(findResult);   
		return;  
    }  
      
    // 2. 循环查找  
    do   
    {  
		//判断是否有子目录
        if (fileInfo.attrib & _A_SUBDIR)    
        {
            //这个语句很重要
            if( (strcmp(fileInfo.name,".") != 0 ) &&(strcmp(fileInfo.name,"..") != 0))   
            {
                string newPath = fileFolder + fileInfo.name;
                listFileByDir(newPath, fileExtension,fileList);
            }
        }
        else  if ( fileInfo.attrib == _A_ARCH)  // 是存档类型文件  
        {
			string fileName = fileFolder + fileInfo.name;
			if(fileExtension.empty())
			{
				fileList->push_back(fileName);
			}else if (!fileName.empty() && fileName.find_last_of(fileExtension))
			{
				fileList->push_back(fileName);
			}
        }
  
    } while (!_findnext(findResult, &fileInfo));    
  
  
    _findclose(findResult);  
}