# DEGgo 0.1.2

## Improvements

* Refactored the internal analysis workflow for improved maintainability.
* Added `analysis.R` and `cleaning.R` modules.
* Renamed `enrichement.R` to `enrichment.R`.
* Simplified and optimized `run_deggo()` internals.
* Improved Gene Ontology enrichment workflow and visualization.
* Enhanced annotation and preprocessing modules.
* Updated function documentation and examples.
* Expanded unit tests and improved package robustness.
* Removed obsolete internal helper documentation.

## Bug fixes

* Fixed documentation inconsistencies.
* Resolved `R CMD check` notes related to documentation and plotting.
* General code cleanup and internal refactoring.

---

# DEGgo 0.1.1

## Improvements

* Improved output organization and reproducibility.
* Enhanced DEG summary tables and visualization.
* Improved automated HTML report generation.
* General code cleanup and documentation updates.

---

# DEGgo 0.1.0

## Initial public release

### New features

* Automated bulk RNA-seq quality control workflow.
* Differential expression analysis using DESeq2, edgeR, or limma.
* Single-contrast and pairwise analysis modes.
* Automated PCA, volcano plot, and heatmap generation.
* Gene Ontology enrichment analysis and visualization.
* Interactive HTML reporting with reproducibility tracking.
* Support for custom OrgDb annotation databases through the `orgdb` argument, enabling analyses of non-model organisms.

### Improvements

* Improved output organization and reproducibility files.
* Enhanced DEG summary tables and visualization.
* Added automated report generation with integrated figures and enrichment results.

---

# DEGgo 0.0.0.9000

## Development

* Initial development version.
* Added bulk RNA-seq QC, differential expression, GO enrichment, and plotting helpers.
* Initial package structure.
