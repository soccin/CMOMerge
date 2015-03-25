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
parser.add_argument("mergeProject", help="other project root directory for merge")
args=parser.parse_args()

baseProject=args.baseProject
mergeProject=args.mergeProject

if not args.tumorType:
	if getTumorType(baseProject) == getTumorType(mergeProject):
		tumorType=getTumorType(baseProject)
	else:
		print >>sys.stderr
		print >>sys.stderr, "Inconsistant tumor types"
		print >>sys.stderr, "   base =", getTumorType(baseProject)
		print >>sys.stderr, "   merge =", getTumorType(mergeProject)
		print >>sys.stderr
		print >>sys.stderr, "Must set tumor type explicitly with --tumorType (-t)"
		print >>sys.stderr
		sys.exit()
else:
	tumorType=args.tumorType

if not args.labName:
	if getLabName(baseProject) == getLabName(mergeProject):
		labName=getLabName(baseProject)
	else:
		print >>sys.stderr
		print >>sys.stderr, "Inconsistant lab names"
		print >>sys.stderr, "   base =", getLabName(baseProject)
		print >>sys.stderr, "   merge =", getLabName(mergeProject)
		print >>sys.stderr
		print >>sys.stderr, "Must set lab name explicitly with --labName (-l)"
		print >>sys.stderr
		sys.exit()
else:
	labName=args.labName

if getInstitutionName(baseProject)!=getInstitutionName(mergeProject):
		print >>sys.stderr
		print >>sys.stderr, "Inconsistant lab names"
		print >>sys.stderr, "   base =", getInstitutionName(baseProject)
		print >>sys.stderr, "   merge =", getInstitutionName(mergeProject)
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
mergePath=Path(mergeProject)

caseFiles=["cases_all.txt","cases_cna.txt","cases_cnaseq.txt","cases_sequenced.txt"]
for caseFile in caseFiles:
	samples=getCaseList(basePath / caseListDir / caseFile) \
			| getCaseList(mergePath / caseListDir / caseFile)
	writeCaseLists(outPath / caseListDir, caseFile, samples, studyId)


rbindFiles=getFileTemplates("""data_clinical.txt
data_mutations_extended.txt
_data_cna_hg19.seg
""")

print "files to rbind", rbindFiles
print "base=", basePath, getStudyId(baseProject)
print "merge=", mergePath, getStudyId(mergeProject)
print "output=", outPath, studyId


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

