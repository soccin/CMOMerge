import sys
from string import Template
from pathlib import *

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
	baseFile=resolvePathToFile(
				basePath,
				fTuple,
				dict(studyId=getStudyId(baseProject)))
	cdrFile=resolvePathToFile(
				cdrPath,
				fTuple,
				dict(studyId=getStudyId(cdrProject)))
	mergedFile=resolvePathToFile(
				outPath,
				fTuple,
				dict(studyId=studyId))
	print "baseFile =", baseFile
	print "cdrFile =", cdrFile
	print "mergedFile =", mergedFile
	print
	mergedTable=rbind(baseFile,cdrFile)
	writeTable(mergedTable,mergedFile)



#######################################################################
#######################################################################
#######################################################################
#######################################################################


"""
cbind:
	data_CNA.txt

meta_CNA.txt
meta_mutations_extended.txt
meta_study.txt
_meta_cna_hg19_seg.txt
"""

