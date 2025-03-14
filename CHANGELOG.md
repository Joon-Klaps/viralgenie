# Viralgenie: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.3dev - [date]

## v0.1.2 - 2025-02-28

Second release of the viralgenie pipeline. Focusing on user experience and bug fixes.

### `Enhancement`

- Set default umitools dedup strategy to cluster ([#126](https://github.com/Joon-Klaps/viralgenie/pull/126))
- Include both krakenreport &nodes.dmp in taxonomy filtering ([#128](https://github.com/Joon-Klaps/viralgenie/pull/128))
- Add Sspace indiv to each assembler seperatly ([#132](https://github.com/Joon-Klaps/viralgenie/pull/132))
- Add read & contig decomplexification using prinseq++  ([#133](https://github.com/Joon-Klaps/viralgenie/pull/133))
- Add option to filter contig clusters based on cumulative read coverage ([#138](https://github.com/Joon-Klaps/viralgenie/pull/138))
- Reffurbish mqc implementation ([#139](https://github.com/Joon-Klaps/viralgenie/pull/139))
- Adding mash-screen output to result table ([#140](https://github.com/Joon-Klaps/viralgenie/pull/140))
- Add logic to allow samples with no reference hits to be analysed ([#141](https://github.com/Joon-Klaps/viralgenie/pull/141))
- Add visualisation for hybrid scaffold ([#143](https://github.com/Joon-Klaps/viralgenie/pull/143))
- Add new module to inculde custom mpileup-vcf file for intra-host analyses ([#151](https://github.com/Joon-Klaps/viralgenie/pull/151))
- Update docs ([#150](https://github.com/Joon-Klaps/viralgenie/pull/150))
- Make custom-mpileup.py postion 1 index based and not 0 index to follow bcftools ([#153](https://github.com/Joon-Klaps/viralgenie/pull/153))
- Update docs for more streamlined docs & figures ([#154](https://github.com/Joon-Klaps/viralgenie/pull/154))
- Add column in custom mpileup - Shannon entropy ([#156](https://github.com/Joon-Klaps/viralgenie/pull/156))
- Constrain -> Constraint & further python script debugging ([#161](https://github.com/Joon-Klaps/viralgenie/pull/161))
- include empty samples in multiqc sample overview ([#162](https://github.com/Joon-Klaps/viralgenie/pull/162))
- Include samtools stats pre dedup & post dedup in overview tables ([#163](https://github.com/Joon-Klaps/viralgenie/pull/163))
- adding prokka for gene detection & annotation ([#165](https://github.com/Joon-Klaps/viralgenie/pull/165))


### `Fixed`

- OOM with longer contigs for nocov_to_reference, uses more RAM now ([#125](https://github.com/Joon-Klaps/viralgenie/pull/125))
- fixing null output from global prefix ([#147](https://github.com/Joon-Klaps/viralgenie/pull/147))
- Fix empty filtered clusters ([#148](https://github.com/Joon-Klaps/viralgenie/pull/148))
- Fixing missing columns from general stats & add general stats sample filtering ([#149](https://github.com/Joon-Klaps/viralgenie/pull/149))
- process.shell template fix ([#157](https://github.com/Joon-Klaps/viralgenie/pull/157)) - see also [nf-core/tools #3416](https://github.com/nf-core/tools/pull/3416)

### `Parameters`
- New parameter mmseqs_cluster_mode default to 0 ([#130](https://github.com/Joon-Klaps/viralgenie/pull/130)) __DEPRECATED__
- Refactor module arguments to pipeline arguments ([#166](https://github.com/Joon-Klaps/viralgenie/pull/166))


## v0.1.1 - 2024-05-08

Initial release of Joon-Klaps/viralgenie, created with the [nf-core](https://nf-co.re/) template.
