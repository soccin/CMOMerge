# CMOMerge
Merge two CMO projects from their portal repository folders

Usage::
CMOMerge/App.sh 
usage: MergePortal [-h] [--tumorType TUMORTYPE] [--labName LABNAME]
                   [--projectNumber PROJECTNUMBER]
                   [--institutionName INSTITUTIONNAME]
                   [--mergeBatches MERGEBATCHES]
                   [--cdrClinicalFile CDRCLINICALFILE]
                   [--cnaGeneList CNAGENELIST]
                   baseProject rightProject

For projects with no merging issues you should simply be able to do

./CMOMerge/App.sh /path/to/portal/repo/batch1 /path/to/portal/repo/batch2


