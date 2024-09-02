# In vivo Macrophages Annotation

The analysis process would take the scripts as the following order:

1. `Data_processing.Rmd`
2. `Myeloids_selection.Rmd`
3. `Macrophages_annotation.Rmd`
4. `Data_processing.Rmd`

Remember to change the value of `rawdata.dir`, `saveddata.dir` etc.

`functions.R` contains all the basic functions that would be used in the scripts. But only one function `load_seurat()` is explictly needed.

Use the function `load_seurat(accession, study.code, rawdata.dir, saveddata.dir)` to load dataset individually. either one of `accession` or `study.code` would be enough to identify datasets.

If there is no locally saved Seurat object corresponding to required study, the function would automatically prompt question asking whether you want to download and process it now. The function will then automatically create folder named after the GEO accession in `rawdata.dir` to save raw data, and finish all the procedure till it Seurat object and will automatically save the object in `saveddata.dir`

The funtcion `load_seurat()` could also be used to load Seurat object after UMAP as long as the object was properly stored following the scripts. Use the parameter `after_UMAP = T` (`F` by default) to activate it.