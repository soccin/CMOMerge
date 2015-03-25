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