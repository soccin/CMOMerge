import sys
import csv
from itertools import izip_longest
from string import Template
from pathlib import *
from lib import *
import os.path
from collections import defaultdict

from globals import *

def get3PathsForMerge(projectList,studyId,outPath,fTuple,updatedFile=None):
	mergeList = []
	for project in projectList:
		if updatedFile and project in updatedFile:
			filename=updatedFile[project]
			print "RPTF:0:", filename
		else:
			filename=resolvePathToFile(
						Path(project),
						fTuple,
						dict(studyId=getStudyId(project)))
		mergeList.append(filename)

	mergedFile=resolvePathToFile(
				outPath,
				fTuple,
				dict(studyId=studyId))
	return (mergeList,mergedFile)

def rbind(fnameList, unionFields=False):
	fieldnames=[]
	data=[]
	for fname in fnameList:
		cin=csv.DictReader(smartOpen(fname),delimiter=CSVDELIM)
		if not fieldnames:
			fieldnames=list(cin.fieldnames)
		elif fieldnames != cin.fieldnames:
			if unionFields:
				for fi in cin.fieldnames:
					if fi not in fieldnames:
						fieldnames.append(fi)
			else:
				print >>sys.stderr, "\n\nfiletools::rbind"
				print >>sys.stderr, "Colnames do not match"
				print >>sys.stderr, fieldnames
				print >>sys.stderr, fname, cin.fieldnames
				print >>sys.stderr
				raise ValueError("Inconsistent colnames")
				sys.exit()
		for rec in cin:
			data.append(dict(rec))
	return (fieldnames,data)


def writeTable(table,outfile):
	fp=smartOpen(str(outfile),mode="w")
	cout=csv.DictWriter(fp,table[0],delimiter=CSVDELIM,lineterminator="\n")
	cout.writeheader()
	for r in table[1]:
		try:
			cout.writerow(r)
		except TypeError:
			print r
			raise
			sys.exit()

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

def mergeCNAData(fnameList,geneList=None):

	print
	print "geneList =",geneList
	print
	genes=set()
	if geneList:
		SDIR=os.path.dirname(os.path.realpath(__file__))
		with open(os.path.join(SDIR, geneList+".genes")) as fp:
			for line in fp:
				for gi in line.strip().split():
					genes.add(gi)
		genes=sorted(genes)

	header=[]
	mergedTable=defaultdict(dict)
	for fname in fnameList:
		cin=csv.DictReader(smartOpen(fname),delimiter=CSVDELIM)
		if cin.fieldnames[0] != "Hugo_Symbol" :
			print >>sys.stderr
			print >>sys.stderr, "Col1 is not Gene Symbols"
			print >>sys.stderr, fname, cin.fieldnames[0]
			print >>sys.stderr
			sys.exit()
		if not header:
			header = list(cin.fieldnames);
		else:
			header += cin.fieldnames[1:]
		tbl=getCNADataTable(cin)
		if not genes:
			genes=sorted(tbl.keys())
		elif not geneList and set(genes)!=set(tbl.keys()):
			print >>sys.stderr
			print >>sys.stderr, "Inconsistent gene sets"
			print >>sys.stderr
			print >>sys.stderr, fname
			print >>sys.stderr
			sys.exit()
		for gi in genes:
			mergedTable[gi].update(tbl[gi])

	data=[]
	for gi in genes:
		rr=dict(mergedTable[gi])
		rr['Hugo_Symbol']=gi
		data.append(rr)
	return (header,data)

def mergeSuppData(fnameList):
	commonSuppFields=set();
	baseSuppFields=[]
	dataRaw=[]
	for fname in fnameList:
		if fname.exists():
			cin=csv.DictReader(smartOpen(fname),delimiter=CSVDELIM)
			if not baseSuppFields:
				baseSuppFields=cin.fieldnames
			if not commonSuppFields:
				commonSuppFields=set(cin.fieldnames)
			else:
				commonSuppFields &= cin.fieldnames
			for r in cin:
				dataRaw.append(r)
		else:
			print "Missing supplemental clinical file: ", fname
			clinicalFile=str(fname).replace("_supp","")
			cin=csv.DictReader(smartOpen(clinicalFile),delimiter=CSVDELIM)
			for r in cin:
				rr=dict()
				rr["SAMPLE_ID"]=r["SAMPLE_ID"]
				rr["PATIENT_ID"]=r["PATIENT_ID"]
				dataRaw.append(rr)

	header=[]
	for fi in baseSuppFields:
		if fi in commonSuppFields:
			header.append(fi)

	data=[]
	for dataRec in dataRaw:
		rr=dict()
		for fi in header:
			if fi in dataRec:
				rr[fi]=dataRec[fi]
			else:
				rr[fi]="na"
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
	print >>fp, outdata.encode('utf-8')
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
