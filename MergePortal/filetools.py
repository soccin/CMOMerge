import sys
import csv
from itertools import izip_longest
from string import Template
from pathlib import *
from lib import *
import os
from collections import defaultdict
import fnmatch

from globals import *

def getPathsForMerge(projectList,studyId,outPath,fTuple,updatedFile=None):
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


def getPathsForMergeRegEx(projectList, pattern, recursive=False):
	mergeList = []
	for project in projectList:
		projectPath = str(Path(project))
		for root, dirnames, filenames in os.walk(projectPath):
			for filename in fnmatch.filter(filenames, pattern):
				mergeList.append(os.path.join(root, filename))
			if(not recursive):
				break # only search for the top level of the directory
	return mergeList


def rbind(fnameList, unionFields, sample_to_group):
	fieldnames=[]
	data=[]
	include_samples = set()
	sample_id_field = None
	if sample_to_group and len(sample_to_group) > 0:
		include_samples = set(sample_to_group.keys())
		if sample_to_group and len(sample_to_group) > 0: # filter samples
			if fnameList and (fnameList[0].name == "data_clinical.txt" or fnameList[0].name == "data_gene_matrix.txt" or fnameList[0].name.startswith("data_clinical_supp_")):
				sample_id_field = "SAMPLE_ID"
			elif fnameList and (fnameList[0].name == "data_mutations_extended.txt" or fnameList[0].name == "data_fusions.txt"):
				sample_id_field = "Tumor_Sample_Barcode"
			elif fnameList and fnameList[0].name.endswith("_data_cna_hg19.seg"):
				sample_id_field = "ID"
			else:
				raise ValueError("Need to filter samples for virtual project and do not know sample column name for file " + fnameList[0].name if fnameList else "?")
				sys.exit()
	for fname in fnameList:
		print fname
		fh=smartOpen(fname)
		skipComments(fh)
		cin=csv.DictReader(fh,delimiter=CSVDELIM)
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
			if not sample_id_field or (sample_id_field in rec and rec[sample_id_field] in include_samples):
				if sample_id_field and sample_id_field in rec and rec[sample_id_field] in include_samples and (fname.name == "data_clinical.txt" or fnameList[0].name.startswith("data_clinical_supp_")):
					rec["PATIENT_ID"] = sample_to_group[rec[sample_id_field]]
				data.append(dict(rec))
	return (fieldnames,data)


def writeTable(table,outfile,replace_cancer_type=None):
	if replace_cancer_type and "CANCER_TYPE" in table[0] and "ONCOTREE_CODE" in table[0]:
		table[0].remove("CANCER_TYPE")
		for r in table[1]:
			if "CANCER_TYPE" in r and ("ONCOTREE_CODE" not in r or not r["ONCOTREE_CODE"]):
				r["ONCOTREE_CODE"] = r["CANCER_TYPE"]
			if "CANCER_TYPE" in r:
				del r["CANCER_TYPE"]

	if (len(table[1]) > 0) or (outfile.name == "data_mutations_extended.txt"):
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
	else:
		print "WARN: %s will have no data, not writing file" % (outfile)

fixGeneNames={"FAM123B":"AMER1"}

def getCNADataTable(cin, include_samples):
	table=dict()
	for r in cin:
		gene=r["Hugo_Symbol"]
		if gene in fixGeneNames:
			gene=fixGeneNames[gene]
		rr=dict(r)
		rr.pop("Hugo_Symbol")
		if len(include_samples) > 0:
			for field in r.keys():
				if field != "Hugo_Symbol" and field not in include_samples:
					del rr[field]
		table[gene]=rr
	return table

def mergeCNAData(fnameList,geneList=None,include_samples=set()):

	print
	print "geneList =",geneList
	print

	# First read in all tables
	cnaTables=[]
	header=[]
	for fname in fnameList:
		cin=csv.DictReader(smartOpen(fname),delimiter=CSVDELIM)
		if cin.fieldnames[0] != "Hugo_Symbol" :
			print >>sys.stderr
			print >>sys.stderr, "Col1 is not Gene Symbols"
			print >>sys.stderr, fname, cin.fieldnames[0]
			print >>sys.stderr
			sys.exit()
		if not header:
			if len(include_samples) > 0:
				header = list(set(["Hugo_Symbol"]) | (set(cin.fieldnames) & include_samples))
			else:
				header = list(cin.fieldnames);
		else:
			if len(include_samples) > 0:
				header += list(set(cin.fieldnames[1:]) & include_samples) 
			else:
				header += cin.fieldnames[1:]
		# remove data not in header
		#if len(include_samples) > 0:
		#	for row in cin:
		#		for field in cin.fieldnames[1:]:
		#			if field not in header:
		#				print >>sys.stderr, "deleting", field, "from", header
		#				del row[field]
		#				del cin.fieldnames[field]
		cnaTables.append(getCNADataTable(cin, include_samples))

	genes=set()
	if geneList:
		SDIR=os.path.dirname(os.path.realpath(__file__))
		with open(os.path.join(SDIR, geneList+".genes")) as fp:
			for line in fp:
				for gi in line.strip().split():
					genes.add(gi)
	else:
		# Get union of all genes from all batches
		for tbl in cnaTables:
			genes.update(tbl.keys())
	genes=sorted(genes)

	mergedTable=defaultdict(dict)
	for tbl in cnaTables:
		samples=tbl[tbl.keys()[0]].keys()
		nullRow=dict([(x,"NA") for x in samples])
		for gi in genes:
			mergedTable[gi].update(tbl.get(gi,nullRow))

	data=[]
	for gi in genes:
		rr=dict(mergedTable[gi])
		rr['Hugo_Symbol']=gi
		data.append(rr)
	return (header,data)

def mergeSuppData(fnameList, sample_to_group):
	commonSuppFields=set();
	baseSuppFields=[]
	dataRaw=[]
	include_samples = set()
	if sample_to_group:
		include_samples = set(sample_to_group.keys())
		sample_id_field = None
		if sample_to_group: # filter samples
			if fname.name == "data_clinical.txt":
				sample_id_field = "SAMPLE_ID"
			elif fname.name == "data_mutations_extended.txt":
				sample_id_field = "Tumor_Sample_Barcode"
			elif fname.name.endswith("_data_cna_hg19.seg"):
				sample_id_field = "ID"
			else:
				raise ValueError("Need to filter samples for virtual project and do not know sample column name for file " + fname)
				sys.exit()
	for fname in fnameList:
		if fname.exists():
			cin=csv.DictReader(smartOpen(fname),delimiter=CSVDELIM)
			if not baseSuppFields:
				baseSuppFields=list(cin.fieldnames)
			if not commonSuppFields:
				commonSuppFields=set(cin.fieldnames)
			else:
				commonSuppFields &= set(cin.fieldnames)
			for r in cin:
				if not sample_id_field or r[sample_id_field] in include_samples:
					if r[sample_id_field] in include_samples:
						r["PATIENT_ID"] = sample_to_group[r[sample_id_field]]
					dataRaw.append(r)
		else:
			print "Missing supplemental clinical file: ", fname
			clinicalFile=str(fname).replace("_supp","")
			cin=csv.DictReader(smartOpen(clinicalFile),delimiter=CSVDELIM)
			for r in cin:
				if not sample_id_field or r[sample_id_field] in include_samples:
					rr=dict()
					if r[sample_id_field] in include_samples:
						r["PATIENT_ID"] = sample_to_group[r[sample_id_field]]
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



def getCaseList(path, include_samples):
	caseList=set()
	with path.open() as fp:
		for line in fp:
			if line.startswith("case_list_ids:"):
				project_samples = set(line.strip().split()[1:])
				if len(include_samples) > 0:
					return include_samples & project_samples # return intersection of sets
				else:
					return project_samples
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

def skipComments(fileHandle):
	startReadingHere = fileHandle.tell()
	line = fileHandle.readline()
	while line.startswith("#"):
		startReadingHere = fileHandle.tell()
		line = fileHandle.readline() 

	# not a comment, reset to start of line
	fileHandle.seek(startReadingHere)
	return
