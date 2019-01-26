from string import Template

metaFiles={
    "meta_study.txt":Template(
"""type_of_cancer: $tumorType
cancer_study_identifier: $studyId
name: $name
date_of_last_update: $lastUpdateDate
description: $description
groups: $groups
short_name: $upperTumorType (CMO $projectNumber)"""
        ),
    "meta_CNA.txt":Template(
"""cancer_study_identifier: $studyId
stable_id: ${studyId}_cna
profile_name: $profile_name
profile_description: $profile_description
genetic_alteration_type: $genetic_alteration_type
datatype: $datatype"""
        ),
    "meta_mutations_extended.txt":Template(
"""cancer_study_identifier: $studyId
stable_id: ${studyId}_mutations
profile_name: $profile_name
profile_description: $profile_description
genetic_alteration_type: $genetic_alteration_type
datatype: $datatype"""
        ),
    "_meta_cna_hg19_seg.txt":Template(
"""cancer_study_identifier: $studyId
reference_genome_id: $reference_genome_id
description: $description
data_filename: ${studyId}_data_cna_hg19.seg"""
        )
}

metaFilesOptional={
   "meta_timeline.txt":Template(
"""cancer_study_identifier: $studyId
"""
),
  "meta_fusions.txt":Template(
"""cancer_study_identifier: $studyId
stable_id: ${studyId}_mutations
datatype: $datatype
genetic_alteration_type: $genetic_alteration_type
show_profile_in_analysis_tab: $show_profile_in_analysis_tab
profile_description: $profile_description
profile_name: $profile_name"""
)
}
