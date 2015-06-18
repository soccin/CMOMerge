# CMOMerge
Merge multiple CMO projects from their portal repository folders

## Usage

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

## Copy Number Gene Lists

If there is an incompatibility in the genes in the copy number data you will get an:

```
Inconsistent gene sets
```

You will then you need to specify a common gene set to use. This is done with the

```bash
--cnaGeneList
```

option and currently the following gene lists are available:

* ```impact341.genes```: The genes in IMPACT341
* ```hemepactA.genes```: The intersection of hemepactV2 and hemepactV3


