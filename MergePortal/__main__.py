import sys
from string import Template
from pathlib import *
import re
import datetime

from filetools import *
from lib import *

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--tumorType","-t",help="Set tumor type")
parser.add_argument("--labName","-l",help="Set lab name")
parser.add_argument("--projectNumber","-p",help="Set project number")
parser.add_argument("baseProject", help="Base project root directory for merge")
parser.add_argument("rightProject", help="other project root directory for merge")
args=parser.parse_args()

baseProject=args.baseProject
cdrProject=args.rightProject

if not args.tumorType:
	if getTumorType(baseProject) == getTumorType(cdrProject):
		tumorType=getTumorType(baseProject)
	else:
		print >>sys.stderr
		print >>sys.stderr, "Inconsistant tumor types"
		print >>sys.stderr, "   base =", getTumorType(baseProject)
		print >>sys.stderr, "   merge =", getTumorType(cdrProject)
		print >>sys.stderr
		print >>sys.stderr, "Must set tumor type explicitly with --tumorType (-t)"
		print >>sys.stderr
		sys.exit()
else:
	tumorType=args.tumorType

if not args.labName:
	if getLabName(baseProject) == getLabName(cdrProject):
		labName=getLabName(baseProject)
	else:
		print >>sys.stderr
		print >>sys.stderr, "Inconsistant lab names"
		print >>sys.stderr, "   base =", getLabName(baseProject)
		print >>sys.stderr, "   merge =", getLabName(cdrProject)
		print >>sys.stderr
		print >>sys.stderr, "Must set lab name explicitly with --labName (-l)"
		print >>sys.stderr
		sys.exit()
else:
	labName=args.labName

if getInstitutionName(baseProject)!=getInstitutionName(cdrProject):
		print >>sys.stderr
		print >>sys.stderr, "Inconsistant lab names"
		print >>sys.stderr, "   base =", getInstitutionName(baseProject)
		print >>sys.stderr, "   merge =", getInstitutionName(cdrProject)
		print >>sys.stderr
		print >>sys.stderr, "Must set lab name explicitly with --labName (-l)"
		print >>sys.stderr
		sys.exit()

institutionName=getInstitutionName(baseProject)
if args.projectNumber:
	projectNumber=args.projectNumber
else:
	projectNumber=getProjectNumber(baseProject)

studyId="_".join([tumorType,institutionName,labName,projectNumber])
print "Merged study id =", studyId

outPath=Path("/".join(["_merge",tumorType,institutionName,labName,projectNumber]))
caseListDir=Path("case_lists")

if outPath.exists():
	print >>sys.stderr, "\nOutput folder", outPath, "exists"
	print >>sys.stderr, "Will not overwrite\n"
	sys.exit()
else:
	outPath.mkdir(parents=True)
	(outPath / caseListDir).mkdir()

basePath=Path(baseProject)
cdrPath=Path(cdrProject)

caseFiles=["cases_all.txt","cases_cna.txt","cases_cnaseq.txt","cases_sequenced.txt"]
for caseFile in caseFiles:
	samples=getCaseList(basePath / caseListDir / caseFile) \
			| getCaseList(cdrPath / caseListDir / caseFile)
	writeCaseLists(outPath / caseListDir, caseFile, samples, studyId)


rbindFiles=getFileTemplates("""
	data_clinical.txt
	data_mutations_extended.txt
	_data_cna_hg19.seg
""")

for fTuple in rbindFiles:
	print "-" * 80
	print "fileSuffix to rbind =", fTuple
	(baseFile,cdrFile,mergedFile)=get3PathsForMerge(baseProject,cdrProject,
													studyId,outPath,fTuple)
	print "baseFile =", baseFile
	print "cdrFile =", cdrFile
	print "mergedFile =", mergedFile
	print
	mergedTable=rbind(baseFile,cdrFile)
	writeTable(mergedTable,mergedFile)

cnaTuple=("data_CNA.txt",None)
(baseFile,cdrFile,mergedFile)=get3PathsForMerge(baseProject,cdrProject,
												studyId,outPath,cnaTuple)
mergedTable=mergeCNAData(baseFile,cdrFile)
writeTable(mergedTable,mergedFile)

from templates import *

today=str(datetime.date.today())
newData=dict(
	studyId=studyId,tumorType=tumorType,upperTumorType=tumorType.upper(),
	projectNumber=projectNumber,lastUpdateDate=today
	)

for metaFile in metaFiles:
	print "Merging", metaFile
	fTuple=getFileTemplates(metaFile)[0]
	baseFile=resolvePathToFile(basePath,fTuple,dict(studyId=getStudyId(baseProject)))
	print baseFile
	baseData=parseMetaData(baseFile)
	baseData.update(newData)
	if "name" in baseData:
		baseData["name"]=baseData["name"].replace(getProjectNumber(baseProject),projectNumber)
	if "description" in baseData:
		baseData["description"]=re.sub(r"2\d\d\d-\d\d-\d\d",today,baseData["description"])
	outFile=resolvePathToFile(outPath,fTuple,dict(studyId=studyId))
	print outFile
	writeMetaFile(outFile,metaFiles[metaFile].substitute(baseData))
	print

#######################################################################
#######################################################################
#######################################################################
#######################################################################


"""
meta_mutations_extended.txt
_meta_cna_hg19_seg.txt
"""

