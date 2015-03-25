from pathlib import *
from filetools import *

def getTumorType(projectPath):
	return projectPath.strip("/").split("/")[-4]

def getLabName(projectPath):
	return projectPath.strip("/").split("/")[-2]

def getInstitutionName(projectPath):
	return projectPath.strip("/").split("/")[-3]

def getProjectNumber(projectPath):
	return projectPath.strip("/").split("/")[-1]

def getStudyId(projectPath):
	# Using new pathlib Path object (which overloads "/")
	# for path concatentation
    data=parseMetaData( Path(projectPath) / "meta_study.txt" )
    return data["cancer_study_identifier"]

def parseMetaData(fname):
    fp=smartOpen(fname)
    data=dict()
    for line in fp:
        line=line.strip()
        pos=line.find(": ")
        if pos>-1:
            data[line[:(pos)]]=line[(pos+2):]
        else:
            raise ValueError("Invalid meta data line %s" % (line))

    return data

def smartOpen(pathType,mode="r"):
    if isinstance(pathType,PosixPath):
        fp=pathType.open(mode=mode)
    elif isinstance(pathType,str):
        fp=open(pathType,mode=mode)
    elif isinstance(pathType,file):
        fp=pathType
    else:
        raise ValueError("Invalid filepath type <%s>" % (type(pathType)))
    return fp

