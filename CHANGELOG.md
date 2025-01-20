# Viralgenie: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1dev - [date]

Initial release of Joon-Klaps/viralgenie, created with the [nf-core](https://nf-co.re/) template.

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

### `Fixed`

- OOM with longer contigs for lowcov_to_reference, uses more RAM now ([#125](https://github.com/Joon-Klaps/viralgenie/pull/125))
- fixing null output from global prefix ([#147](https://github.com/Joon-Klaps/viralgenie/pull/147))
- Fix empty filtered clusters ([#148](https://github.com/Joon-Klaps/viralgenie/pull/148))
- Fixing missing columns from general stats & add general stats sample filtering ([#149](https://github.com/Joon-Klaps/viralgenie/pull/149))

### `Parameters`
- New parameter mmseqs_cluster_mode default to 0 ([#130](https://github.com/Joon-Klaps/viralgenie/pull/130))
