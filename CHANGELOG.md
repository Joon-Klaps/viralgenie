# Viralgenie: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1dev - [date]

Initial release of Joon-Klaps/viralgenie, created with the [nf-core](https://nf-co.re/) template.

### `Enhancement`

- Set default umitools dedup strategy to cluster ([#126](https://github.com/Joon-Klaps/viralgenie/pull/126))
- Include both krakenreport &nodes.dmp in taxonomy filtering ([#128](https://github.com/Joon-Klaps/viralgenie/pull/128))
- Include coverage plot & subset contig results in mqc report ([#129](https://github.com/Joon-Klaps/viralgenie/pull/129))
- Add Sspace indiv to each assembler seperatly ([#132](https://github.com/Joon-Klaps/viralgenie/pull/132))
- Add read & contig decomplexification using prinseq++  ([#133](https://github.com/Joon-Klaps/viralgenie/pull/133))
- Add option to filter contig clusters based on cumulative read coverage ([#138](https://github.com/Joon-Klaps/viralgenie/pull/138))
- Reffurbish mqc implementation ([#139](https://github.com/Joon-Klaps/viralgenie/pull/139))

### `Fixed`

- OOM with longer contigs for lowcov_to_reference, uses more RAM now ([#125](https://github.com/Joon-Klaps/viralgenie/pull/125))

### `Parameters`
- New parameter mmseqs_cluster_mode default to 0 ([#130](https://github.com/Joon-Klaps/viralgenie/pull/130))
