# A data-driven testing module for gramtools

## Data production

From a fasta reference and a vcf against it, we can:
	* Slice both
	* Randomly pick derived alleles in the sliced vcf --> makes a mosaic sample
	* Randomly mutate the mosaic sample --> makes a mutant sample
	* Generate reads from the sample, using art_ilmn.

This can be done for any number of samples, specified at command line.

## Validation

* The mosaic sample is stored as vcf, and can serve to bench how well quasimap + infer are able to retrieve the correct mosaic.

* The mutant sample is stored as a vcf, rebased against the original (sliced) reference, and can serve to bench how well discover is able to retrieve the variation on top of the mosaic.

## Pipelines

* fullRun: runs a full build/quasimap/infer/discover on a sample
* Multisample : runs multi sample pipeline

## TODOs

[] Combining the mutant vcfs + the original (sliced) vcf into genotyped multi-sample vcf. This can serve as a reference of what multi sample pipeline should ideally output.
[] Automated assessment using gramtools outputs versus the truths.
