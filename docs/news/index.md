# Changelog

## DEGgo 1.0.0

### New features

- Initial stable release.
- Added support for differential expression analysis with DESeq2, edgeR,
  and limma.
- Added single-contrast and pairwise differential expression workflows.
- Added automated quality-control and input-validation workflows.
- Added publication-ready PCA plots, volcano plots, heatmaps, and
  summary tables.
- Added Gene Ontology enrichment analysis and visualization.
- Added automated HTML report generation.
- Added PowerPoint export.
- Added reproducibility reports containing analysis parameters and
  session information.
- Added example datasets, tutorials, and complete analysis workflows.

### Experimental design assessment

- Added design_qc() for assessing experimental structure before
  differential expression analysis.
- Added validation of count matrices, metadata, sample identifiers,
  groups, batches, and experimental-unit columns.
- Added automated transformation of count matrices using VST, rlog, or
  log2-CPM.
- Added principal-component analysis for identifying factors associated
  with global expression variation.
- Added optional PERMANOVA analysis using the vegan package.
- Added assessment of treatment, batch, tissue, sex, donor, patient,
  cage, animal, and other design variables.
- Added detection of repeated experimental units when unit_col is
  supplied or conservatively detected.
- Added identification of paired, repeated-measures, nested, and
  potentially pseudoreplicated designs.
- Added automated recommendations for DESeq2 or dream based on the
  detected experimental structure.
- Added suggested model formulas and alternative analysis strategies.
- Added the deggo_design_qc result class and a dedicated print method
  for concise design summaries.

### Rhythmicity analysis

- Added run_deggo_rhythmicity() for circadian and time-course
  rhythmicity analysis.
- Added MetaCycle-based rhythmicity detection.
- Added single-component cosinor regression with automatic or manual
  fitting.
- Added optional differential-rhythmicity testing between two groups.
- Added rhythmicity diagnostic plots and gene-symbol annotation support.
- Added structured rhythmicity result tables and reproducibility
  parameters.

### Improvements

- Improved plot_volcano() documentation.
- Improved count and metadata preprocessing.
- Added automatic detection of common gene identifier columns.
- Added optional removal of Ensembl version suffixes.
- Added optional aggregation of duplicated gene identifiers.
- Improved sample matching and metadata validation.
- Improved handling of repeated-measures and experimental-unit designs.
- Improved recommendations for avoiding pseudoreplication.
- Improved output organization and reproducibility.
- Enhanced differential-expression summary tables and visualizations.
- Improved Gene Ontology enrichment summaries.
- Improved automated HTML report generation.
- Removed unintended global-variable dependencies.
- Improved CRAN portability by replacing non-ASCII characters in R
  source files.
- Resolved package warnings and code/documentation issues detected by R
  CMD check.
- General code cleanup, refactoring, testing, and documentation updates.

## DEGgo 0.1.2

### Improvements

- Refactored the internal analysis workflow for improved
  maintainability.
- Added analysis.R and cleaning.R modules.
- Renamed enrichement.R to enrichment.R.
- Simplified and optimized run_deggo() internals.
- Improved Gene Ontology enrichment workflow and visualization.
- Enhanced annotation and preprocessing modules.
- Updated function documentation and examples.
- Expanded unit tests and improved package robustness.
- Removed obsolete internal helper documentation.

### Bug fixes

- Fixed documentation inconsistencies.
- Resolved R CMD check notes related to documentation and plotting.
- Performed general code cleanup and internal refactoring.

## DEGgo 0.1.0

### Initial public release

#### New features

- Added an automated bulk RNA-seq quality-control workflow.
- Added differential expression analysis using DESeq2, edgeR, or limma.
- Added single-contrast and pairwise analysis modes.
- Added automated PCA, volcano plot, and heatmap generation.
- Added Gene Ontology enrichment analysis and visualization.
- Added interactive HTML reporting with reproducibility tracking.
- Added support for custom OrgDb annotation databases through the orgdb
  argument, enabling analyses of non-model organisms.

#### Improvements

- Improved output organization and reproducibility files.
- Enhanced differential-expression summary tables and visualization.
- Added automated report generation with integrated figures and
  enrichment results.

## DEGgo 0.0.0.9000

### Development

- Initial development version.
- Added bulk RNA-seq quality control, differential expression, Gene
  Ontology enrichment, and plotting helpers.
- Added the initial package structure.
