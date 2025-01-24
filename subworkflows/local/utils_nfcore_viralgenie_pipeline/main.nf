//
<<<<<<< HEAD
// Subworkflow with functionality specific to the nf-core/viralgenie pipeline
=======
// Subworkflow with functionality specific to the Joon-Klaps/viralgenie pipeline
>>>>>>> TEMPLATE
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schemas/input.json"))
        .map{
            meta, read1, read2 ->
            def single_end = read1 && !read2
            if (single_end) {
                return [meta + [sample: meta.id, single_end: single_end] , [read1]]
            }
            else {
                return [meta + [sample: meta.id, single_end: single_end] , [read1, read2]]
            }
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "Viralgenie (Klaps et al.)",
            "nf-core (Ewels et al. 2020)",
            "Nextflow (Di Tommaso et al. 2017)",
            "Bbduk (Bushnell 2022)",
            "BCFtools (Danecek et al. 2021)",
            "BLAST+ (Camacho et al. 2009)",
            "Bowtie2 (Langmead and Salzberg 2012)",
            "BWA-MEM (Li 2013)",
            "BWA-MEM2 (Vasimuddin et al. 2019)",
            "CD-HIT (Fu et al. 2012)",
            "CheckV (Nayfach et al. 2021)",
            "FastQC (Andrews 2010)",
            "fastp (Chen et al. 2018)",
            "HUMID (Laros and van den Berg)",
            "iVar (Grubaugh et al. 2019)",
            "Kaiju (Menzel et al. 2016)",
            "Kraken2 (Wood et al. 2019)",
            "Leiden Algorithm (Traag et al. 2019)",
            "Mash (Ondov et al. 2016)",
            "MEGAHIT (Li et al. 2016)",
            "Minimap2 (Li 2018)",
            "MMseqs2 (Steinegger and Söding 2017)",
            "Mosdepth (Pedersen and Quinlan 2018)",
            "MultiQC (Ewels et al. 2016)",
            "Picard (Broad Institute)",
            "QUAST (Gurevich et al. 2013)",
            "SAMtools (Li 2011)",
            "SPAdes (Bankevich et al. 2012)",
            "SSPACE Basic (Boetzer et al. 2011)",
            "Trimmomatic (Bolger et al. 2014)",
            "Trinity (Haas et al. 2013)",
            "UMI-tools (Smith et al. 2017)",
            "vRhyme (Kieft et al. 2022)",
            "VSEARCH (Rognes et al. 2016)",
            "Anaconda (Anaconda 2016)",
            "Bioconda (Grüning et al. 2018)",
            "BioContainers (da Veiga Leprevost et al. 2017)",
            "Docker (Merkel 2014)",
            "Singularity (Kurtzer et al. 2017)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Klaps J, Lemey P, Kafetzopoulou L. Viralgenie: A metagenomics analysis pipeline for eukaryotic viruses. __Github__ https://github.com/Joon-Klaps/viralgenie.</li>",
            "<li>Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.</li>",
            "<li>Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.</li>",
            "<li>Bushnell B. (2022) BBMap, URL: http://sourceforge.net/projects/bbmap/</li>",
            "<li>Danecek, Petr et al. “Twelve years of SAMtools and BCFtools.” GigaScience vol. 10,2 (2021): giab008. doi:10.1093/gigascience/giab008</li>",
            "<li>Camacho, Christiam et al. “BLAST+: architecture and applications.” BMC bioinformatics vol. 10 421. 15 Dec. 2009, doi:10.1186/1471-2105-10-421</li>",
            "<li>Langmead, Ben, and Steven L Salzberg. “Fast gapped-read alignment with Bowtie 2.” Nature methods vol. 9,4 357-9. 4 Mar. 2012, doi:10.1038/nmeth.1923</li>",
            "<li>Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2.</li>",
            "<li>M. Vasimuddin, S. Misra, H. Li and S. Aluru, 'Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems,' 2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS), Rio de Janeiro, Brazil, 2019, pp. 314-324, doi: 10.1109/IPDPS.2019.00041.</li>",
            "<li>Fu, Limin et al. “CD-HIT: accelerated for clustering the next-generation sequencing data.” Bioinformatics (Oxford, England) vol. 28,23 (2012): 3150-2. doi:10.1093/bioinformatics/bts565</li>",
            "<li>Nayfach, Stephen et al. “CheckV assesses the quality and completeness of metagenome-assembled viral genomes.” Nature biotechnology vol. 39,5 (2021): 578-585. doi:10.1038/s41587-020-00774-7</li>",
            "<li>Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online].</li>",
            "<li>Chen, Shifu et al. “fastp: an ultra-fast all-in-one FASTQ preprocessor.” Bioinformatics (Oxford, England) vol. 34,17 (2018): i884-i890. doi:10.1093/bioinformatics/bty560</li>",
            "<li>Laros J, van den Berg R, __Github__ https://github.com/jfjlaros/HUMID</li>",
            "<li>Grubaugh, Nathan D et al. “An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar.” Genome biology vol. 20,1 8. 8 Jan. 2019, doi:10.1186/s13059-018-1618-7</li>",
            "<li>Menzel, Peter et al. “Fast and sensitive taxonomic classification for metagenomics with Kaiju.” Nature communications vol. 7 11257. 13 Apr. 2016, doi:10.1038/ncomms11257</li>",
            "<li>Wood, Derrick E., Jennifer Lu, and Ben Langmead. 2019. Improved Metagenomic Analysis with Kraken 2. Genome Biology 20 (1): 257. doi: 10.1186/s13059-019-1891-0.</li>",
            "<li>Traag, V A et al. “From Louvain to Leiden: guaranteeing well-connected communities.” Scientific reports vol. 9,1 5233. 26 Mar. 2019, doi:10.1038/s41598-019-41695-z</li>",
            "<li>Ondov, Brian D et al. “Mash: fast genome and metagenome distance estimation using MinHash.” Genome biology vol. 17,1 132. 20 Jun. 2016, doi:10.1186/s13059-016-0997-x</li>",
            "<li>Li, Dinghua et al. “MEGAHIT v1.0: A fast and scalable metagenome assembler driven by advanced methodologies and community practices.” Methods (San Diego, Calif.) vol. 102 (2016): 3-11. doi:10.1016/j.ymeth.2016.02.020</li>",
            "<li>Li, Heng. “Minimap2: pairwise alignment for nucleotide sequences.” Bioinformatics (Oxford, England) vol. 34,18 (2018): 3094-3100. doi:10.1093/bioinformatics/bty191</li>",
            "<li>Steinegger, Martin, and Johannes Söding. “MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.” Nature biotechnology vol. 35,11 (2017): 1026-1028. doi:10.1038/nbt.3988</li>",
            "<li>Pedersen, Brent S, and Aaron R Quinlan. “Mosdepth: quick coverage calculation for genomes and exomes.” Bioinformatics (Oxford, England) vol. 34,5 (2018): 867-868. doi:10.1093/bioinformatics/btx699</li>",
            "<li>Ewels, Philip et al. “MultiQC: summarize analysis results for multiple tools and samples in a single report.” Bioinformatics (Oxford, England) vol. 32,19 (2016): 3047-8. doi:10.1093/bioinformatics/btw354</li>",
            "<li>Gurevich, Alexey et al. “QUAST: quality assessment tool for genome assemblies.” Bioinformatics (Oxford, England) vol. 29,8 (2013): 1072-5. doi:10.1093/bioinformatics/btt086</li>",
            "<li>Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93. doi: 10.1093/bioinformatics/btr509. Epub 2011 Sep 8. PMID: 21903627; PMCID: PMC3198575.</li>",
            "<li>Bankevich, Anton et al. “SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing.” Journal of computational biology : a journal of computational molecular cell biology vol. 19,5 (2012): 455-77. doi:10.1089/cmb.2012.0021</li>",
            "<li>Boetzer, Marten et al. “Scaffolding pre-assembled contigs using SSPACE.” Bioinformatics (Oxford, England) vol. 27,4 (2011): 578-9. doi:10.1093/bioinformatics/btq683</li>",
            "<li>Bolger, Anthony M et al. “Trimmomatic: a flexible trimmer for Illumina sequence data.” Bioinformatics (Oxford, England) vol. 30,15 (2014): 2114-20. doi:10.1093/bioinformatics/btu170</li>",
            "<li>Haas, Brian J et al. “De novo transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis.” Nature protocols vol. 8,8 (2013): 1494-512. doi:10.1038/nprot.2013.084</li>",
            "<li>Smith, Tom et al. “UMI-tools: modeling sequencing errors in Unique Molecular Identifiers to improve quantification accuracy.” Genome research vol. 27,3 (2017): 491-499. doi:10.1101/gr.209601.116</li>",
            "<li>Kieft, Kristopher et al. “vRhyme enables binning of viral genomes from metagenomes.” Nucleic acids research vol. 50,14 (2022): e83. doi:10.1093/nar/gkac341</li>",
            "<li>Rognes, Torbjørn et al. “VSEARCH: a versatile open source tool for metagenomics.” PeerJ vol. 4 e2584. 18 Oct. 2016, doi:10.7717/peerj.2584</li>",
            "<li>Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.</li>",
            "<li>Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.</li>",
            "<li>da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.</li>",
            "<li>Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.</li>",
            "<li>Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.</li>",
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""


    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

def createFileChannel(param) {
    return param ? Channel.fromPath(param, checkIfExists: true).collect() : []
}

def createChannel(dbPath, dbName, skipFlag) {
    return dbPath && skipFlag ? Channel.fromPath(dbPath, checkIfExists: true).map { db -> [[id: dbName], db] } : Channel.empty()
}

def filterContigs(contig, min_len, n_100) {
    contig
        .map { meta, fasta -> [ meta, fasta, WorkflowCommons.getLengthAndAmbigous( fasta ) ] }
        .branch { meta, fasta, stats ->
            pass: stats.contig_size >= min_len.toInteger() && stats.n_100 <= n_100.toInteger()
                return [ meta, fasta ]
            fail: stats.contig_size < min_len.toInteger() || stats.n_100 > n_100.toInteger()
                return [ meta, fasta, stats ]}
}

def failedContigsToMultiQC(tsv_data, min_len, n_100) {
    tsv_data
        .map { meta, fasta, stats -> ["$meta.id\t$meta.sample\t$meta.cluster_id\t$meta.previous_step\t$stats.contig_size\t$stats.n_100"] }
        .collect()
        .map { tsv ->
            WorkflowCommons.multiqcTsvFromList(
                tsv,
                ['Id','sample name', 'cluster','step','contig size', 'N\'s %'],
                [
                    "id: 'failed_contig_quality'",
                    "anchor: 'WARNING: Filtered contigs'",
                    "section_name: 'Failed contig quality'",
                    "format: 'tsv'",
                    "description: 'Contigs that are not of minimum size ${min_len} or have more then ${n_100} ambigous bases per 100 kbp were filtered out'",
                    "plot_type: 'table'"
                ]
            )
        }
}

def failedMappedReadsToMultiQC(tsv_data, min_mapped_reads) {
    tsv_data
        .map { meta, bam, mapped_reads ->
            ["$meta.id\t$meta.sample\t$meta.cluster_id\t$meta.previous_step\t$mapped_reads"]
            }
        .collect()
        .map { tsv ->
            WorkflowCommons.multiqcTsvFromList(tsv,
                ['id','sample name', 'cluster','step','mapped reads'],
                [
                    "id: 'failed_mapped'",
                    "anchor: 'WARNING: Filtered contigs'",
                    "section_name: 'Minimum mapped reads'",
                    "format: 'tsv'",
                    "description: 'Contigs that did not have more then ${min_mapped_reads} mapped reads were filtered out'",
                    "plot_type: 'table'"
                ]
            )
        }
}

def noBlastHitsToMultiQC(tsv_data, assemblers) {
    tsv_data
        .map { meta, txt, fasta ->
            def n_fasta = fasta.countFasta()
            ["$meta.sample\t$n_fasta"]}
        .collect()
        .map { tsv ->
            WorkflowCommons.multiqcTsvFromList(tsv,
                ['sample name', "number of contigs"],
                [
                    "id: 'samples_without_blast_hits'",
                    "anchor: 'WARNING: Filtered samples'",
                    "section_name: 'Samples without blast hits'",
                    "format: 'tsv'",
                    "description: 'Samples that did not have any blast hits for their contigs (using ${assemblers}) were not included in further assembly polishing'",
                    "plot_type: 'table'"
                ]
            )
        }
}

def lowReadSamplesToMultiQC(tsv_data, min_trimmed_reads) {
    tsv_data
        .map { meta, read_count -> ["$meta.sample\t$read_count"] }
        .collect()
        .map { tsv ->
            WorkflowCommons.multiqcTsvFromList(
                tsv,
                ['Sample', "Number of reads"],
                [
                    "id: 'samples_low_reads'",
                    "anchor: 'WARNING: Filtered samples'",
                    "section_name: 'Samples with to few reads'",
                    "format: 'tsv'",
                    "description: 'Samples that did not have the minimum number of reads (<${min_trimmed_reads}) after trimming, complexity filtering & host removal'",
                    "plot_type: 'table'"
                ]
            )
        }
}

def noContigSamplesToMultiQC(tsv_data, assemblers) {
    tsv_data
        .map { meta, fasta ->
            def n_fasta = fasta.countFasta()
            ["$meta.sample\t$n_fasta"]
        }
        .collect()
        .map { tsv ->
            WorkflowCommons.multiqcTsvFromList(
                tsv,
                ['sample name', "number of contigs"],
                [
                    "id: 'samples_without_contigs'",
                    "anchor: 'WARNING: Filtered samples'",
                    "section_name: 'Samples without contigs'",
                    "format: 'tsv'",
                    "description: 'Samples that did not have any contigs (using ${assemblers}) were not included in further assembly polishing'",
                    "plot_type: 'table'"
                ]
            )
        }
}

