# CMOMerge
Merge multiple CMO projects from their portal repository folders

Usage::

```
CMOMerge/App.sh
usage: MergePortal [-h] [--tumorType TUMORTYPE] [--labName LABNAME]
                   [--projectNumber PROJECTNUMBER]
                   [--institutionName INSTITUTIONNAME]
                   [--mergeBatches MERGEBATCHES]
                   [--cnaGeneList CNAGENELIST]
                   [--project PROJECTPATH:CLINICALFILE] (:CLINICALFILE is optional)
```

For projects with no merging issues you should simply be able to do

```bash
./CMOMerge/App.sh \
    --project /path/to/portal/repo/batch1 \
    --project /path/to/portal/repo/batch2
```

and this will create a merged project in the _merge subfolder. E.g.:

```bash
./CMOMerge/App.sh \
    --project bic-mskcc/cecc/cmo/levine2/5529_b \
    --project bic-mskcc/cecc/mskcc/levine2/5529
```

will create:

> _merge/cecc/cbe/levine2/5529_b


