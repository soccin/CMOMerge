import sys
from string import Template
from pathlib import *
import re
import datetime
import subprocess
import os.path
from collections import defaultdict

from filetools import *
from lib import *
from templates import *
from gene_panel import *


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
parser.add_argument("--virtualGroupFile","-v",default="",help="Virtual group file")
parser.add_argument("--metaStudyFile","-s",default="",help="Replacement meta study file")
parser.add_argument("--skipCaseLists", action="store_true", default=False, help="Do not add any case lists")
args=parser.parse_args()

if not args.project:
	print >>sys.stderr, "--project is required"
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

# if virtual project
sample_to_group = defaultdict(str)

# read clusters file
if args.virtualGroupFile:
	print "Reading:", args.virtualGroupFile
	with open(args.virtualGroupFile, "r") as virtual_group_file:
		if virtual_group_file.readline(): # skip first line
			for line in virtual_group_file:
				line = line.strip("\n")
				cells = line.split("\t")
				sample_to_group[cells[1].strip()] = cells[0].strip()

include_samples = set(sample_to_group.keys())

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
		print >>sys.stderr, "	 ",project, " = ", getProjectNumber(project)
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
batches = mergeBatches.split(",")
shortMergeBatches = mergeBatches
if len(batches) > 3:
	shortMergeBatches = "%s,,%d_total" % (batches[0], len(batches))
print "mergeBatches =", mergeBatches
print "shortMergeBatches=", shortMergeBatches
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

if not args.skipCaseLists:
	caseFiles=["cases_all.txt","cases_cna.txt","cases_cnaseq.txt","cases_sequenced.txt"]
	samples_to_projects=defaultdict(set)
	for caseFile in caseFiles:
		samples = set()
		for project in projectList:
			projectPath=Path(project)
			# filter this if virtual project by passing set(sample_to_group.keys())
			project_samples = getCaseList(projectPath / caseListDir / caseFile, include_samples)
			samples |= project_samples
			for sample in project_samples:
				samples_to_projects[sample].add(project)
		writeCaseLists(outPath / caseListDir, caseFile, samples, studyId)

	# we can only check this if we read the case lists (at least for now)
	replicate_samples = [ (s, p) for (s, p) in samples_to_projects.iteritems() if len(p) > 1]
	if replicate_samples:
		print >>sys.stderr, "\nReplicate samples found:"
		for (sample, projects) in replicate_samples:
			print >>sys.stderr, "Sample", sample, "is found in projects:", ", ".join(projects) , "\n"
		sys.exit(64) # this will be the exit code for this specific issue

rbindFiles=getFileTemplates("""
	data_clinical.txt
	data_mutations_extended.txt
	_data_cna_hg19.seg
""")

for fTuple in rbindFiles:
	print "-" * 80
	print "fileSuffix to rbind =", fTuple
	# TODO virtual project
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
	mergedTable=rbind(mergeList,unionFieldNames,sample_to_group)
	writeTable(mergedTable,mergedFile,replace_cancer_type=(fTuple[0]=="data_clinical.txt"))

write_gene_panel_files(studyId,outPath,projectList,include_samples)

cnaTuple=("data_CNA.txt",None)
(mergeList,mergedFile)=getPathsForMerge(projectList,studyId,outPath,cnaTuple)
mergedTable=mergeCNAData(mergeList,geneList,include_samples)
writeTable(mergedTable,mergedFile)

# if this is not a virtual study and we have a meta study file
# use the fields in it, rather than just using that file as the final file
newMetaStudyData=dict()
if not args.virtualGroupFile and args.metaStudyFile:
	newMetaStudyData=parseMetaData(args.metaStudyFile) 

today=str(datetime.date.today())
newData=dict()
newData.update(dict(
	studyId=studyId,tumorType=tumorType,upperTumorType=tumorType.upper(),
	projectNumber="%s [%s]" % (projectNumber,projectTag),
	lastUpdateDate=today
	))

for metaFile in metaFiles:
	fTuple=getFileTemplates(metaFile)[0]
	outFile=resolvePathToFile(outPath,fTuple,dict(studyId=studyId))
	print outFile
	baseFile=resolvePathToFile(basePath,fTuple,dict(studyId=getStudyId(baseProject)))
	# if virtual study just copy this file
	if args.virtualGroupFile and args.metaStudyFile and baseFile.name == "meta_study.txt":
		print "Copying", metaFile
		shutil.copyfile(args.metaStudyFile, str(outFile))
	else:
		print "Merging", metaFile
		print baseFile
		baseData=parseMetaData(baseFile)
		baseData.update(newData)
		if baseFile.name == "meta_study.txt":
			baseData.update(newMetaStudyData)
		if "name" in baseData:
			baseData["name"]=baseData["name"].replace(getProjectNumber(baseProject),newData["projectNumber"])
			#baseData["name"]=re.sub(r"^(.+?) - ([^ ]+) (.+)$", r"\1 - %s \3" % (labName.capitalize()[:-1]), baseData["name"])
		if "description" in baseData:
			baseData["description"]=re.sub(r"2\d\d\d-\d\d-\d\d",today,baseData["description"])
			pos=baseData["description"].find(" (BATCHES:")
			displayBatches = shortMergeBatches
			if fTuple[0] == "meta_study.txt":
				displayBatches = mergeBatches
			if pos>-1:
				baseData["description"]=baseData["description"][:pos]+" (BATCHES: %s)" % displayBatches
			elif args.virtualGroupFile:
				baseData["description"]+=" (ver %s;)" % (GIT_VERSION)
			else:
				baseData["description"]+=" (ver %s; BATCHES: %s)" % (GIT_VERSION, displayBatches)
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
	print "Copying (and filtering if virtual study)", clinSuppFile
	outputFile=os.path.join(str(outPath), os.path.basename(clinSuppFile))
	if os.path.isfile(outputFile):
		print >>sys.stderr, "\nFile ", outputFile, "already exists"
		print >>sys.stderr, "There might be a file name conflict in multiple projects\n"
		sys.exit()
	#shutil.copyfile(clinSuppFile, outputFile)
	table=rbind([Path(clinSuppFile)],False,sample_to_group)
	writeTable(table,Path(outputFile))
print

# TODO add back once we can filter on patient id (REAL patient id not virtual one)
if not args.skipCaseLists:
	timelinePattern="data_timeline_*.txt"
	mergeList=getPathsForMergeRegEx(projectList, timelinePattern)
	for timelineFile in mergeList:
		print "Copying", timelineFile
		if sample_to_group and len(sample_to_group) > 0:
			print >>sys.stderr, "Do not currently replace patient ids with group ids for timeline file in virtual study\n"
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

fusionDataFile="data_fusions.txt"
dataMergeList=getPathsForMergeRegEx(projectList,fusionDataFile)
if len(dataMergeList)>0:
	unionFieldNames=False
	dataOutputFile=os.path.join(str(outPath),fusionDataFile)
	fusionMetaFile="meta_fusions.txt"
	for mf in dataMergeList:
		print "inputFile =", mf
	print "mergedFile =", dataOutputFile
	print
	mergedTable=rbind([Path(mf) for mf in dataMergeList],unionFieldNames,sample_to_group)
	writeTable(mergedTable,dataOutputFile)

	metaMergeList=getPathsForMergeRegEx(projectList,fusionMetaFile)
	# write meta file here too if we find any
	if len(metaMergeList)>0:
		metaOutputFile=os.path.join(str(outPath),fusionMetaFile)
		print "Creating", metaOutputFile
		print "Using", metaMergeList[0]
		baseData=parseMetaData(metaMergeList[0])
		baseData.update(newData)
		writeMetaFile(metaOutputFile, metaFilesOptional[fusionMetaFile].substitute(baseData))

