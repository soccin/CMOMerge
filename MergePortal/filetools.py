import sys
import csv
from itertools import izip_longest
from string import Template
from pathlib import *
from lib import *

from globals import *

def get3PathsForMerge(baseProject,cdrProject,studyId,outPath,fTuple):
	baseFile=resolvePathToFile(
				Path(baseProject),
				fTuple,
				dict(studyId=getStudyId(baseProject)))
	cdrFile=resolvePathToFile(
				Path(cdrProject),
				fTuple,
				dict(studyId=getStudyId(cdrProject)))
	mergedFile=resolvePathToFile(
				outPath,
				fTuple,
				dict(studyId=studyId))
	return (baseFile,cdrFile,mergedFile)

def rbind(fname1,fname2):
	cin1=csv.DictReader(smartOpen(fname1),delimiter=CSVDELIM)
	cin2=csv.DictReader(smartOpen(fname2),delimiter=CSVDELIM)
	if cin1.fieldnames != cin2.fieldnames:
		print >>sys.stderr
		print >>sys.stderr, "Colnames do not match"
		print >>sys.stderr, fname1, cin1.fieldnames
		print >>sys.stderr, fname2, cin2.fieldnames
		print >>sys.stderr
		sys.exit()

	data=[]
	for rec in cin1:
		data.append(dict(rec))
	for rec in cin2:
		data.append(dict(rec))

	return (cin1.fieldnames,data)

def writeTable(table,outfile):
	fp=smartOpen(str(outfile),mode="w")
	cout=csv.DictWriter(fp,table[0],delimiter=CSVDELIM,lineterminator="\n")
	cout.writeheader()
	for r in table[1]:
		cout.writerow(r)

fixGeneNames={"FAM123B":"AMER1"}

def getCNADataTable(cin):
	table=dict()
	for r in cin:
		gene=r["Hugo_Symbol"]
		if gene in fixGeneNames:
			gene=fixGeneNames[gene]
		rr=dict(r)
		rr.pop("Hugo_Symbol")
		table[gene]=rr
	return table

def mergeCNAData(fname1,fname2):
	cin1=csv.DictReader(smartOpen(fname1),delimiter=CSVDELIM)
	cin2=csv.DictReader(smartOpen(fname2),delimiter=CSVDELIM)

	if cin1.fieldnames[0]!=cin2.fieldnames[0] or cin1.fieldnames[0] != "Hugo_Symbol" :
		print >>sys.stderr
		print >>sys.stderr, "Col1 do not match or not Gene Symbols"
		print >>sys.stderr, fname1, cin1.fieldnames[0]
		print >>sys.stderr, fname2, cin2.fieldnames[0]
		print >>sys.stderr
		sys.exit()

	header=cin1.fieldnames+cin2.fieldnames[1:]

	tbl1=getCNADataTable(cin1)
	tbl2=getCNADataTable(cin2)

	if set(tbl1.keys())!=set(tbl2.keys()):
		print >>sys.stderr
		print >>sys.stderr, "Inconsistent gene sets"
		print >>sys.stderr, fname1, set(tbl1.keys)
		print >>sys.stderr, fname2, set(tbl2.keys)
		print >>sys.stderr
		sys.exit()

	data=[]
	genes=sorted(tbl1.keys())
	for gi in genes:
		rr=dict(tbl1[gi])
		rr.update(tbl2[gi])
		rr['Hugo_Symbol']=gi
		data.append(rr)

	return (header,data)

def getCaseList(path):
	caseList=set()
	with path.open() as fp:
		for line in fp:
			if line.startswith("case_list_ids:"):
				return set(line.strip().split()[1:])
	raise ValueError("Missing case_list_ids line in %s" % path)

caseFileData={
	"cases_all.txt":dict(
		stable_id="_all",
		case_list_category="all_cases_in_study",
		case_list_name="All Tumors",
		case_list_description="All tumor samples"
	),
	"cases_cna.txt":dict(
		stable_id="_cna",
		case_list_category="all_cases_with_cna_data",
		case_list_name="Tumors CNA",
		case_list_description="All tumors with CNA data"
	),
	"cases_cnaseq.txt":dict(
		stable_id="_cnaseq",
		case_list_category="all_cases_with_mutation_and_cna_data",
		case_list_name="Tumors with sequencing and CNA data",
		case_list_description="All tumor samples that have CNA and sequencing data"
	),
	"cases_sequenced.txt":dict(
		stable_id="_sequenced",
		case_list_category="all_cases_with_mutation_data",
		case_list_name="Sequenced Tumors",
		case_list_description="All sequenced tumors"
	)
}

def writeCaseLists(outDir, caseFile, samples, studyId):
	data=caseFileData[caseFile]
	with (outDir / caseFile).open("wb") as fp:
		print >>fp, "cancer_study_identifier:", studyId
		print >>fp, "stable_id:", studyId+data["stable_id"]
		print >>fp, "case_list_category:", data["case_list_category"]
		print >>fp, "case_list_name:", data["case_list_name"]
		print >>fp, "case_list_description:", data["case_list_description"]
		print >>fp, "case_list_ids:", "\t".join(map(str,sorted(samples)))

def writeMetaFile(outfile,outdata):
	fp=smartOpen(str(outfile),"w")
	print >>fp, outdata
	fp.close()

def getFileTemplates(fileLists):
	files=[]
	for fi in fileLists.strip().split():
		if fi.startswith("_"):
			template=Template("${studyId}"+fi)
		else:
			template=None
		files.append((fi,template))
	return files

def resolvePathToFile(path,fnameTuple,templateData=dict()):
	print "RPTF:1:", fnameTuple
	if fnameTuple[1]:
		print "RPTF:2:", templateData, fnameTuple[1].template
		fname=fnameTuple[1].substitute(templateData)
		print "RPTF:3:", fname
	else:
		fname=fnameTuple[0]

	fullPath = path / fname
	return fullPath
