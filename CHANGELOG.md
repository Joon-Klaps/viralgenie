# Viralgenie: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1dev - [date]

Initial release of Joon-Klaps/viralgenie, created with the [nf-core](https://nf-co.re/) template.

### `Enhancement`

- Set default umitools dedup strategy to cluster ([#126](https://github.com/Joon-Klaps/viralgenie/pull/126))
- Include sspace for contig extension ([#123](https://github.com/Joon-Klaps/viralgenie/pull/123))
- Include both krakenreport &nodes.dmp in taxonomy ([#128](https://github.com/Joon-Klaps/viralgenie/pull/128))

### `Fixed`

- OOM with longer contigs for lowcov_to_reference, uses more RAM now ([#125](https://github.com/Joon-Klaps/viralgenie/pull/125))

### `Parameters`
