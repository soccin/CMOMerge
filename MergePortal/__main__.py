import sys
from string import Template
from pathlib import *
import re
import datetime
import subprocess
import os.path

from filetools import *
from lib import *
from templates import *


PDIR=os.path.dirname(os.path.realpath(__file__))
SDIR=os.path.split(PDIR)[0]

#--git-dir=$REPO/.git --work-tree=$REPO
GIT_VERSION=subprocess.check_output(
	["git",
	 "--git-dir=%s/.git" % SDIR,
	 "--work-tree=%s" % SDIR,
	"describe","--always"
	]).strip()
print "Version (%s)" % GIT_VERSION

import argparse
import shutil
parser = argparse.ArgumentParser()
parser.add_argument("--tumorType","-t",help="Set tumor type")
parser.add_argument("--labName","-l",help="Set lab name")
parser.add_argument("--projectNumber","-n",help="Set project number")
parser.add_argument("--institutionName","-i",help="Set institution Name (mskcc)")
parser.add_argument("--mergeBatches","-m",help="Batches in merge")
parser.add_argument("--cnaGeneList",help="Set explicit gene list of CNA merge")
parser.add_argument("--project", action='append', help="project root directory for merge, and optional updated clinical file, seperated by :")
parser.add_argument("--force","-f", action="store_true", default=False, help="Force overwrite")
parser.add_argument("--root",default="",help="Set location of hg repository")
args=parser.parse_args()

if not args.project:
	parser.print_help()
	sys.exit()

projectList=[]
updatedClinicalFile=dict()

for pi in args.project:
	piParse = pi.split(":")
	if len(piParse) == 1:
		projectList.append(piParse[0])
	elif len(piParse) == 2:
		projectList.append(piParse[0])
		updatedClinicalFile[piParse[0]] = piParse[1]
		print "Project: ", piParse[0], "will use updated clinical file: ", piParse[1]
	else:
		print >>sys.stderr
		print >>sys.stderr, "Incorrect --project parameter: ", pi
		print >>sys.stderr
		sys.exit()
baseProject = projectList[0]


if not args.cnaGeneList:
	geneList=None
else:
	geneList=args.cnaGeneList

if not args.tumorType:
	tumorTypeSet = set()
	for project in projectList:
		tumorTypeSet.add(getTumorType(project))
	if len(tumorTypeSet) == 1:
		tumorType=getTumorType(baseProject)
	else:
		print >>sys.stderr
		print >>sys.stderr, "Inconsistant tumor types"
		for project in projectList:
			print >>sys.stderr, project, " = ", getTumorType(project)
		print >>sys.stderr
		print >>sys.stderr, "Must set tumor type explicitly with --tumorType (-t)"
		print >>sys.stderr
		sys.exit()
else:
	tumorType=args.tumorType

tumorType=tumorType.lower()

if not args.labName:
	labNameSet = set()
	for project in projectList:
		labNameSet.add(getLabName(project))
	if len(labNameSet) == 1:
		labName=getLabName(baseProject)
	else:
		print >>sys.stderr
		print >>sys.stderr, "Inconsistant lab names"
		for project in projectList:
			print >>sys.stderr, project, " = ", getLabName(project)
		print >>sys.stderr
		print >>sys.stderr, "Must set lab name explicitly with --labName (-l)"
		print >>sys.stderr
		sys.exit()
else:
	labName=args.labName

if not args.institutionName:
	institutionName="cbe"
else:
	institutionName=args.institutionName

if args.projectNumber:
	projectNumber=args.projectNumber
else:
	print >>sys.stderr
	print >>sys.stderr, "Merge project number must be specified with -n option"
	print >>sys.stderr
	for project in projectList:
		print >>sys.stderr, "   ",project, " = ", getProjectNumber(project)
	print >>sys.stderr
	print >>sys.stderr
	sys.exit()

studyId="_".join([tumorType,institutionName,labName,projectNumber])
print "-" * 80
print "Merged study id =", studyId

if args.mergeBatches:
	mergeBatches=args.mergeBatches
else:
	projectNumSet = set()
	projectTagSet = set()
	projectBatchSet = []
	for project in projectList:
		projectNum=getProjectNumber(project)
		projectTag=projectNum.split("_")[0]
		projectBatch="_".join(projectNum.split("_")[1:]).upper()
		projectBatch="A" if not projectBatch else projectBatch
		projectNumSet.add(projectNum)
		projectTagSet.add(projectTag)
		projectBatchSet.append(projectBatch)
	if len(projectTagSet) > 1:
		baseProjectNum=getProjectNumber(baseProject)
		projectNumSet.remove(baseProjectNum)
		mergeBatches=",".join([baseProjectNum] + list(projectNumSet))
	else:
		mergeBatches=",".join(sorted(projectBatchSet))

print "mergeBatches =", mergeBatches
print "projectList =", projectList
print "projectTag =", projectTag

if args.root=="":
	outPath=Path("/".join([tumorType,institutionName,labName,projectNumber]))
else:
	outPath=Path("/".join([args.root,tumorType,institutionName,labName,projectNumber]))

caseListDir=Path("case_lists")


if outPath.exists():
	if not args.force:
		print >>sys.stderr, "\nOutput folder", outPath, "exists"
		print >>sys.stderr, "Will not overwrite\n"
		sys.exit()
	else:
		print "\n Overwriting", outPath
		import shutil
		shutil.rmtree(str(outPath))

outPath.mkdir(parents=True)
(outPath / caseListDir).mkdir()

basePath=Path(baseProject)


caseFiles=["cases_all.txt","cases_cna.txt","cases_cnaseq.txt","cases_sequenced.txt"]
for caseFile in caseFiles:
	samples = set()
	for project in projectList:
		projectPath=Path(project)
		samples |= getCaseList(projectPath / caseListDir / caseFile)
	writeCaseLists(outPath / caseListDir, caseFile, samples, studyId)


rbindFiles=getFileTemplates("""
	data_clinical.txt
	data_mutations_extended.txt
	_data_cna_hg19.seg
""")

for fTuple in rbindFiles:
	print "-" * 80
	print "fileSuffix to rbind =", fTuple
	if fTuple[0]=="data_clinical.txt":
		unionFieldNames=True
		(mergeList,mergedFile)=getPathsForMerge(projectList,studyId,outPath,fTuple,updatedClinicalFile)
	elif fTuple[0]=="data_mutations_extended.txt":
		unionFieldNames=True
		(mergeList,mergedFile)=getPathsForMerge(projectList,studyId,outPath,fTuple)
	else:
		unionFieldNames=False
		(mergeList,mergedFile)=getPathsForMerge(projectList,studyId,outPath,fTuple)

	for mf in mergeList:
		print "inputFile =", mf
	print "mergedFile =", mergedFile
	print
	mergedTable=rbind(mergeList,unionFieldNames)
	writeTable(mergedTable,mergedFile)


cnaTuple=("data_CNA.txt",None)
(mergeList,mergedFile)=getPathsForMerge(projectList,studyId,outPath,cnaTuple)
mergedTable=mergeCNAData(mergeList,geneList)
writeTable(mergedTable,mergedFile)

today=str(datetime.date.today())
newData=dict(
	studyId=studyId,tumorType=tumorType,upperTumorType=tumorType.upper(),
	projectNumber="%s [%s]" % (projectNumber,projectTag),
	lastUpdateDate=today
	)

for metaFile in metaFiles:
	print "Merging", metaFile
	fTuple=getFileTemplates(metaFile)[0]
	baseFile=resolvePathToFile(basePath,fTuple,dict(studyId=getStudyId(baseProject)))
	print baseFile
	baseData=parseMetaData(baseFile)
	baseData.update(newData)
	if "name" in baseData:
		baseData["name"]=baseData["name"].replace(getProjectNumber(baseProject),newData["projectNumber"])
	if "description" in baseData:
		baseData["description"]=re.sub(r"2\d\d\d-\d\d-\d\d",today,baseData["description"])
		pos=baseData["description"].find(" (BATCHES:")
		if pos>-1:
			baseData["description"]=baseData["description"][:pos]+" (BATCHES: %s)" % mergeBatches
		else:
			baseData["description"]+=" (ver %s; BATCHES: %s)" % (GIT_VERSION, mergeBatches)
	outFile=resolvePathToFile(outPath,fTuple,dict(studyId=studyId))
	print outFile
	writeMetaFile(outFile,metaFiles[metaFile].substitute(baseData))
	print


clinSuppTuple=("data_clinical_supp.txt",None)
(mergeList,mergedFile)=getPathsForMerge(projectList,studyId,outPath,clinSuppTuple)
#
# Only do supp clinical merge if the base file
# exists, mergeSuppData knows now how to deal with
# missing file on cdr side
#
baseFile=mergeList[0]
if baseFile.exists():
	mergedTable=mergeSuppData(mergeList)
	if mergedTable:
		writeTable(mergedTable,mergedFile)
print




clinSuppExtraPattern="data_clinical_supp_*.txt"
mergeList=getPathsForMergeRegEx(projectList, clinSuppExtraPattern)
for clinSuppFile in mergeList:
	print "Copying", clinSuppFile
	outputFile=os.path.join(str(outPath), os.path.basename(clinSuppFile))
	if os.path.isfile(outputFile):
		print >>sys.stderr, "\nFile ", outputFile, "already exists"
		print >>sys.stderr, "There might be a file name conflict in multiple projects\n"
		sys.exit()
	shutil.copyfile(clinSuppFile, outputFile)
print


timelinePattern="data_timeline_*.txt"
mergeList=getPathsForMergeRegEx(projectList, timelinePattern)
for timelineFile in mergeList:
	print "Copying", timelineFile
	outputFile=os.path.join(str(outPath), os.path.basename(timelineFile))
	if os.path.isfile(outputFile):
		print >>sys.stderr, "\nFile ", outputFile, "already exists"
		print >>sys.stderr, "There might be a file name conflict in multiple projects\n"
		sys.exit()
	shutil.copyfile(timelineFile, outputFile)
if(len(mergeList)>0):
	metaFile="meta_timeline.txt"
	outputFile=os.path.join(str(outPath), metaFile)
	print "Creating", outputFile
	writeMetaFile(outputFile, metaFilesOptional[metaFile].substitute(newData))
print





















