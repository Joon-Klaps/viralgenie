# nf-core/viralmetagenome: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.1dev - YYYY-MM-DD

### `Added`

- ([#313](https://github.com/nf-core/viralmetagenome/pull/313)) - Add param `--normalise_reads` to digitally normalise reads with BBNorm before de novo assembly (by @Joon-Klaps)
- ([#317](https://github.com/nf-core/viralmetagenome/pull/317)) - Add `tests/lib/UTILS.groovy` and rewrite the pipeline-level nf-tests as scenario lists, so assertions are defined once instead of copied across seven files (by @Joon-Klaps)
- ([#318](https://github.com/nf-core/viralmetagenome/pull/318)) - Update all nf-core modules and subworkflows (by @Joon-Klaps)
- ([#320](https://github.com/nf-core/viralmetagenome/pull/320)) - Fix [#281](https://github.com/nf-core/viralmetagenome/issues/281) - add opt-in param `--use_host_filtered_reads` to route host-filtered reads into the iterative consensus refinement and final variant-calling mapping steps. Defaults to `false` to preserve existing behaviour (by @Joon-Klaps)
- ([#318](https://github.com/nf-core/viralmetagenome/pull/318)) - Document how to build a custom annotation database from NCBI Virus, and fix stale `customisation/databases.md` cross-references left over from a docs restructure (by @Joon-Klaps)
- ([#318](https://github.com/nf-core/viralmetagenome/pull/318)) - Add param `--annotation_metadata` to read the consensus annotation fields from a csv/tsv table (optionally gzipped) instead of parsing them out of the annotation database fasta headers (by @Joon-Klaps)
- ([#312](https://github.com/nf-core/viralmetagenome/pull/312)) - Fix [#282](https://github.com/nf-core/viralmetagenome/issues/282) - drop unmapped reads from the contig-coverage alignment and unmapped read pairs from the alignments used in polishing and consensus refinement, and add param `--keep_unmapped` to carry them through instead (by @Joon-Klaps)

### `Fixed`

- ([#302](https://github.com/nf-core/viralmetagenome/pull/302)) - Fix [#285](https://github.com/nf-core/viralmetagenome/issues/285) & update modules of the subworkflow `fastq_kraken_kaiju` (by @Joon-Klaps)
- ([#315](https://github.com/nf-core/viralmetagenome/pull/315)) - Fix [#314](https://github.com/nf-core/viralmetagenome/issues/314) - make `cluster.tsv` not an 'intermediate' file, as it contains a lot of info and is useful for downstream analyses (by @Joon-Klaps)
- ([#318](https://github.com/nf-core/viralmetagenome/pull/318)) - Replace the deprecated `cat/cat` module with `find/concatenate`, and `tabix/tabix` module with `htslib/bgziptabix` (by @Joon-Klaps)

### `Dependencies`

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `bcftools` | 1.22        | 1.23.1      |
| `fastp`    | 1.0.1       | 1.3.6       |
| `htslib`   | 1.22.1      | 1.24        |
| `minimap2` | 2.29        | 2.30        |
| `mosdepth` | 0.3.11      | 0.3.14      |
| `multiqc`  | 1.34        | 1.35        |
| `picard`   | 3.4.0       | 3.5.0       |
| `prokka`   | 1.14.6      | 1.15.6      |
| `samtools` | 1.22.1      | 1.24        |
| `snpeff`   | 5.4.0a      | 5.4.0c      |
| `vsearch`  | 2.21.1      | 2.31.0      |

### `Parameters`

| Old parameter | New parameter               |
| ------------- | --------------------------- |
|               | `--keep_unmapped`           |
|               | `--normalise_reads`         |
|               | `--arguments_bbnorm`        |
|               | `--annotation_metadata`     |
|               | `--use_host_filtered_reads` |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.

## v1.1.3 - 2026-07-02

### `Added`

### `Fixed`

- default value of `--reference_pool` from `current` to `C-RVDBv31.0` as RVDB changed their default link from `fasta.gz` to `fasta.tar.gz` which can no longer be directly integrated within the pipeline. ([#298](https://github.com/nf-core/viralmetagenome/pull/298)) (by @Joon-Klaps)

### `Dependencies`

### `Parameters`

| Old parameter | New parameter               |
| ------------- | --------------------------- |
|               | `--use_host_filtered_reads` |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.

## v1.1.2 - 2026-06-03

### `Added`

- Add param `--cluster_with_reference_pool` to optionally include the reference pool in the clustering step ([#278](https://github.com/nf-core/viralmetagenome/pull/278)) (by @nrminor, @Joon-Klaps)
- Template update to nf-core/tools v4.0.2 ([#289](https://github.com/nf-core/viralmetagenome/pull/289)) (by @Joon-Klaps)
- Bump version from dev to `1.1.2` for release ([#292](https://github.com/nf-core/viralmetagenome/pull/292)) (by @Joon-Klaps)

### `Fixed`

- Cast Lazymaps to Hashmaps to become thread safe with various groovy operators ([#291](https://github.com/nf-core/viralmetagenome/pull/291)) (by @aersoares81, @Joon-Klaps)
- Extract Python package versions via `python -c "import X; print(X.__version__)"` instead of `pip show`, so versions are recorded reliably across conda/pip layouts ([#292](https://github.com/nf-core/viralmetagenome/pull/292)) (by @Joon-Klaps)

### `Dependencies`

| Dependency  | Old version | New version |
| ----------- | ----------- | ----------- |
| `nf-schema` | 2.5.1       | 2.7.2       |
| `nft-utils` | 0.0.3       | 1.0.0       |

### `Parameters`

| Old parameter | New parameter                   |
| ------------- | ------------------------------- |
|               | `--cluster_with_reference_pool` |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.

## v1.1.1 - 2026-03-31

### `Added`

- Template update for nf-core/tools v3.5.2 ([#270](https://github.com/nf-core/viralmetagenome/pull/270)) (by @Joon-Klaps)
- Modules update to latest nf-core versions ([#271](https://github.com/nf-core/viralmetagenome/pull/271)) (by @Joon-Klaps)
- Update local modules to use topics ([#272](https://github.com/nf-core/viralmetagenome/pull/272)) (by @Joon-Klaps)
- Add better output file descriptions for prokka [#274](https://github.com/nf-core/viralmetagenome/pull/274) (by @Joon-Klaps)
- Better item descriptors (`it -> [it[0], it[1]` to `meta, fasta , _foo -> [meta,fasta]`) for better readability ([#265](https://github.com/nf-core/viralmetagenome/pull/265)) (by @Joon-Klaps)
- Include the new parameter `--transpose_overview_tables` , fixes [#220](https://github.com/nf-core/viralmetagenome/issues/220) ([#275](https://github.com/nf-core/viralmetagenome/pull/275)) (by @Joon-Klaps)
- Bump version from dev to `1.1.1` for release ([#277](https://github.com/nf-core/viralmetagenome/pull/276)) (by @Joon-Klaps)

### `Fixed`

- Follow the rules of strict syntax health [#265](https://github.com/nf-core/viralmetagenome/pull/265) (by @Joon-Klaps)
- Moved the module specifics from base.config to modules.config ([#265](https://github.com/nf-core/viralmetagenome/pull/265))
- Fix bug from regex in input schema still allowing underscores in sample and group names ([#269](https://github.com/nf-core/viralmetagenome/pull/269)) - (by @Joon-Klaps)

### `Dependencies`

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `snpeff`   | 5.3.0a      | 5.4.0a      |

### `Parameters`

| Old parameter | New parameter                 |
| ------------- | ----------------------------- |
|               | `--transpose_overview_tables` |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.

## v1.1.0 - 2026-01-23

### `Added`

- Add subworkflow tests of `fasta_fastq_clust` and `fasta_contig_clust` [#253](https://github.com/nf-core/viralmetagenome/pull/253) (by @Joon-Klaps)
- Add subworkflow tests for `preprocessing_illumina` and `fastq_fastqc_umitools_trimmomatic` [#255](https://github.com/nf-core/viralmetagenome/pull/255) (by @Joon-Klaps)

### `Fixed`

- Fix memory issues for `network_cluster` [#105](https://github.com/nf-core/viralmetagenome/issues/105) by switching to Clusty [#251](https://github.com/nf-core/viralmetagenome/pull/251) (by @Joon-Klaps)
- Improve memory efficiency and speed of FASTA sequence extraction [#252](https://github.com/nf-core/viralmetagenome/pull/252) (by @Joon-Klaps)
- Fix bug of `clusty` not handling single genome distance files well [#253](https://github.com/nf-core/viralmetagenome/pull/253) (by @Joon-Klaps)
- Fix reported trimmomatic bug [#254](https://github.com/nf-core/viralmetagenome/issues/254) 'no such variable "trim_read_count"' [#255](https://github.com/nf-core/viralmetagenome/pull/255) (by @Joon-Klaps)
- Fix descrepancy of documentation and actual arguments [#256](https://github.com/nf-core/viralmetagenome/pull/256) (by @Joon-Klaps)
- Fix performance issues with `extract_cluster.py` for large datasets [#259](https://github.com/nf-core/viralmetagenome/pull/259) (by @Joon-Klaps)

### `Dependencies`

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `clusty`   | N/A         | 1.2.2       |

### `Deprecated`

## v1.0.1 - 2025-12-03

### `Added`

- Update modules ([#237](https://github.com/nf-core/viralmetagenome/pull/237)) (by @Joon-Klaps)

### `Fixed`

- Template merge from nf-core/tools v3.5.1 ([#236](https://github.com/nf-core/viralmetagenome/pull/236)) (by @Joon-Klaps)
- Fix bug with numeric sample ids in custom_multiqc.py module ([#243](https://github.com/nf-core/viralmetagenome/pull/243)) (by @Joon-Klaps)
- Fix runtime environment migration of `network_cluster` from custom docker repo to seqera containers ([#247](https://github.com/nf-core/viralmetagenome/pull/247)) (by @Joon-Klaps)

### `Dependencies`

| Dependency         | Old version | New version |
| ------------------ | ----------- | ----------- |
| `bcftools`         | 1.21        | 1.22        |
| `blast`            | 2.16.0      | 2.17.0      |
| `bwa-mem2`         | 2.2.1       | 2.3         |
| `coreutils`        | 9.4         | 9.5         |
| `fastp`            | 0.24.0      | 1.0.1       |
| `htslib`           | 1.21        | 1.22.1      |
| `kraken2`          | 2.1.5       | 2.1.6       |
| `krakentools`      | 1.2         | 1.2.1       |
| `mmseqs2`          | 17.b804f    | 18.8cc5c    |
| `mosdepth`         | 0.3.10      | 0.3.11      |
| `picard`           | 3.3.0       | 3.4.0       |
| `samtools`         | 1.21        | 1.22.1      |
| `snpeff`           | 5.1         | 5.3.0a      |
| `umitools extract` | 1.1.5       | 1.1.6       |

### `Deprecated`

## v1.0.0 - 2025-10-04

Initial release of nf-core/viralmetagenome, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Add argument to control the maximum mpileup depth in `custom_mpileup.py` script ([#176](https://github.com/nf-core/viralmetagenome/pull/176)) (by @Joon-Klaps)
- Removing redundant `samtools_sort` after `BAM_DEDUPLICATE` ([#177](https://github.com/nf-core/viralmetagenome/pull/177)) (by @Joon-Klaps)
- Removing `seqkit replace` and move logic to `blast_filter.py` ([#178](https://github.com/nf-core/viralmetagenome/pull/178)) (by @Joon-Klaps)
- Giving Kaiju more memory ([#179](https://github.com/nf-core/viralmetagenome/pull/179)) (by @Joon-Klaps)
- Add option for sample merging based on group column in samplesheet ([#180](https://github.com/nf-core/viralmetagenome/pull/180)) (by @Joon-Klaps)
- Creating testdatasets using nf-core repo ([#183](https://github.com/nf-core/viralmetagenome/pull/183)) (by @Joon-Klaps)
- Add option to annotate snps with Snpeff for mapping constraint route ([#186](https://github.com/Joon-Klaps/viralgenie/pull/186)) (by @Joon-Klaps)
- Add `nf-tests` for `test` profile ([#189](https://github.com/nf-core/viralmetagenome/pull/189)) (by @Joon-Klaps)
- Update docs ([#200](https://github.com/nf-core/viralmetagenome/pull/200)) (by @Joon-Klaps)
- Template update for nf-core/tools v3.3.2 ([#202](https://github.com/nf-core/viralmetagenome/pull/202)) (by @Joon-Klaps)
- Add option for local desktop configuration profile ([#227](https://github.com/nf-core/viralmetagenome/pull/227)) (by @Joon-Klaps)
- Added option to blacklist certain NCBI accessions (`--blacklist`) from reference pool (`--reference_pool`) ([#228](https://github.com/nf-core/viralmetagenome/pull/228)) (by @Joon-Klaps)
- Add support for both string and integer inputs in samplesheets through `validation.lentientMode = true` ([#230](https://github.com/nf-core/viralmetagenome/pull/230)) (by @Joon-Klaps)

### `Fixed`

- Fix local module `ivar_variants_to_vcf` to handle empty tsv files ([#197](https://github.com/Joon-Klaps/viralmetagenome/pull/197/)) (by @Joon-Klaps)
- Migrate lib dir functions to utils_nfcore_viralmetagenome_pipeline ([#194](https://github.com/nf-core/viralmetagenome/pull/194)) (by @Joon-Klaps)
- Fix inconsistent dependency versions across modules ([#208](https://github.com/nf-core/viralmetagenome/pull/208)) (by @Joon-Klaps)
- Fix conda issue unrecognized arguments: --mkdir ([#210](https://github.com/nf-core/viralmetagenome/pull/210)) (by @Joon-Klaps)
- Fix writing no sequence for select_reference.py to the first reference of the multifasta ([#214](https://github.com/nf-core/viralmetagenome/pull/214)) (by @Joon-Klaps)
- Fix main language detection to ignore generated files ([#224](https://github.com/nf-core/viralmetagenome/pull/224)) (by @Joon-Klaps)

### `Dependencies`

### `Deprecated`

- Refactor `params.skip_annotation` to `params.skip_consensus_annotation`([#181](https://github.com/nf-core/viralmetagenome/pull/181)) (by @Joon-Klaps)
- Deprecate `params.skip_nocov_to_reference` ([#212](https://github.com/nf-core/viralmetagenome/pull/212)) (by @Joon-Klaps)
- Deprecate `BWAMEM` as mapping tool ([#212](https://github.com/nf-core/viralmetagenome/pull/212)) (by @Joon-Klaps)

## v0.1.2 - 2025-02-28

Second release of thenf-core/viralmetagenome pipeline. Focusing on user experience and bug fixes.

### `Added`

- Set default umitools dedup strategy to cluster ([#126](https://github.com/nf-core/viralmetagenome/pull/126)) (by @Joon-Klaps)
- Include both krakenreport &nodes.dmp in taxonomy filtering ([#128](https://github.com/nf-core/viralmetagenome/pull/128)) (by @Joon-Klaps)
- Add Sspace indiv to each assembler seperatly ([#132](https://github.com/nf-core/viralmetagenome/pull/132)) (by @Joon-Klaps)
- Add read & contig decomplexification using prinseq++ ([#133](https://github.com/nf-core/viralmetagenome/pull/133)) (by @Joon-Klaps)
- Add option to filter contig clusters based on cumulative read coverage ([#138](https://github.com/nf-core/viralmetagenome/pull/138)) (by @Joon-Klaps)
- Reffurbish mqc implementation ([#139](https://github.com/nf-core/viralmetagenome/pull/139)) (by @Joon-Klaps)
- Adding mash-screen output to result table ([#140](https://github.com/nf-core/viralmetagenome/pull/140)) (by @Joon-Klaps)
- Add logic to allow samples with no reference hits to be analysed ([#141](https://github.com/nf-core/viralmetagenome/pull/141)) (by @Joon-Klaps)
- Add visualisation for hybrid scaffold ([#143](https://github.com/nf-core/viralmetagenome/pull/143)) (by @Joon-Klaps)
- Add new module to inculde custom mpileup-vcf file for intra-host analyses ([#151](https://github.com/nf-core/viralmetagenome/pull/151)) (by @Joon-Klaps)
- Update docs ([#150](https://github.com/nf-core/viralmetagenome/pull/150)) (by @Joon-Klaps)
- Make custom-mpileup.py postion 1 index based and not 0 index to follow bcftools ([#153](https://github.com/nf-core/viralmetagenome/pull/153)) (by @Joon-Klaps)
- Update docs for more streamlined docs & figures ([#154](https://github.com/nf-core/viralmetagenome/pull/154)) (by @Joon-Klaps)
- Add column in custom mpileup - Shannon entropy ([#156](https://github.com/nf-core/viralmetagenome/pull/156)) (by @Joon-Klaps)
- Constrain -> Constraint & further python script debugging ([#161](https://github.com/nf-core/viralmetagenome/pull/161)) (by @Joon-Klaps)
- include empty samples in multiqc sample overview ([#162](https://github.com/nf-core/viralmetagenome/pull/162)) (by @Joon-Klaps)
- Include samtools stats pre dedup & post dedup in overview tables ([#163](https://github.com/nf-core/viralmetagenome/pull/163)) (by @Joon-Klaps)
- adding prokka for gene detection & annotation ([#165](https://github.com/nf-core/viralmetagenome/pull/165)) (by @Joon-Klaps)

### `Fixed`

- OOM with longer contigs for nocov_to_reference, uses more RAM now ([#125](https://github.com/nf-core/viralmetagenome/pull/125)) (by @Joon-Klaps)
- fixing null output from global prefix ([#147](https://github.com/nf-core/viralmetagenome/pull/147)) (by @Joon-Klaps)
- Fix empty filtered clusters ([#148](https://github.com/nf-core/viralmetagenome/pull/148)) (by @Joon-Klaps)
- Fixing missing columns from general stats & add general stats sample filtering ([#149](https://github.com/nf-core/viralmetagenome/pull/149)) (by @Joon-Klaps)
- process.shell template fix ([#157](https://github.com/nf-core/viralmetagenome/pull/157)) (by @Joon-Klaps) - see also [nf-core/tools #3416](https://github.com/nf-core/tools/pull/3416)

### `Parameters`

- New parameter mmseqs_cluster_mode default to 0 ([#130](https://github.com/nf-core/viralmetagenome/pull/130)) (by @Joon-Klaps) **DEPRECATED**
- Refactor module arguments to pipeline arguments ([#166](https://github.com/nf-core/viralmetagenome/pull/166)) (by @Joon-Klaps)

## v0.1.1 - 2024-05-08

Initial release of nf-core/viralmetagenome, created with the [nf-core](https://nf-co.re/) template.
