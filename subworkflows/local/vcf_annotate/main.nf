// Based on https://github.com/nf-core/viralrecon/blob/master/subworkflows/local/variants_ivar.nf

include { SNPEFF_BUILD          } from '../../../modules/local/snpeff_build/main'
include { SNPEFF_SNPEFF         } from '../../../modules/nf-core/snpeff/snpeff/main'
include { SNPSIFT_EXTRACTFIELDS } from '../../../modules/local/snpsift_extractfields/main'

workflow VCF_ANNOTATE {
    take:
    ch_vcf_ref // channel: [ val(meta), [ vcf], [ref]]

    main:

    // A gff must be supplied & the genome shouldn't have undergone selection
    ch_vcf_ref = ch_vcf_ref.filter { meta, _vcf, _ref -> meta.gff != null && meta.gff != [] && !meta.selection }

    ch_gff = ch_vcf_ref
        .map { meta, _vcf, ref -> [[id: meta.cluster_id], ref, meta.gff] }
        .unique()

    SNPEFF_BUILD(
        ch_gff
    )

    ch_snpeff_in = ch_vcf_ref
        .map { meta, vcf, ref -> [[id: meta.cluster_id], meta, ref, vcf] }
        .combine(SNPEFF_BUILD.out.db, by: [0])
        .multiMap { db_id, meta, _ref, vcf, db, config ->
            vcf: [meta, vcf]
            db: db_id.id
            cache: [meta, db, config]
        }

    SNPEFF_SNPEFF(
        ch_snpeff_in.vcf,
        ch_snpeff_in.db,
        ch_snpeff_in.cache,
    )
    ch_vcf_ann = SNPEFF_SNPEFF.out.vcf

    SNPSIFT_EXTRACTFIELDS(
        ch_vcf_ann
    )

    emit:
    vcf      = ch_vcf_ann                    // channel: [ val(meta), [ tsv ] ]
    txt      = SNPSIFT_EXTRACTFIELDS.out.txt // channel: [ val(meta), [ txt ] ]
}
