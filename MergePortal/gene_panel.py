import os.path
from pathlib import *
from lib import *
import csv

from globals import *

def write_gene_panel_files(studyId, studyPath, projectList, include_samples):
	# some files might be missing, but we can get it from data_clinical instead
	samples_to_gene_panels = {}
	for project in projectList:
		data_gene_panel = os.path.join(str(Path(project)), "data_gene_matrix.txt")
		if not os.path.isfile(data_gene_panel):
			data_clinical_file = os.path.join(str(Path(project)), "data_clinical.txt")
			print "WARN: missing %s get data from %s instead" % (data_gene_panel, data_clinical_file)

			# read data clinical file
			fh=smartOpen(data_clinical_file)
			samples = csv.DictReader(fh,delimiter=CSVDELIM)
			for sample in samples:
				if len(include_samples) == 0 or sample["SAMPLE_ID"] in include_samples:
					samples_to_gene_panels[sample["SAMPLE_ID"]] = sample["GENE_PANEL"]
		else:
			# read data from gene panel file instead 
			fh=smartOpen(data_gene_panel)
			samples = csv.DictReader(fh,delimiter=CSVDELIM)
			for sample in samples:
				if len(include_samples) == 0 or sample["SAMPLE_ID"] in include_samples:
					samples_to_gene_panels[sample["SAMPLE_ID"]] = sample["mutations"]

	# write gene panel meta file
	meta_gene_panel = os.path.join(str(studyPath), "meta_gene_matrix.txt")
	write_meta_gene_matrix(meta_gene_panel, studyId)
	# write gene panel data file
	data_gene_panel = os.path.join(str(studyPath), "data_gene_matrix.txt")
	write_data_gene_matrix(data_gene_panel, samples_to_gene_panels)

def write_data_gene_matrix(fileName, samples_to_gene_panels):
	with open(fileName, 'w') as fh:
		fh.write("SAMPLE_ID\tmutations\n")
		for sample, panel in samples_to_gene_panels.iteritems():
			fh.write("%s\t%s\n" % (sample, panel))

def write_meta_gene_matrix(fileName, studyId):
	with open(fileName, 'w') as fh:
		fh.write("cancer_study_identifier: %s\n" % (studyId))
		fh.write("data_filename: data_gene_matrix.txt\n")
