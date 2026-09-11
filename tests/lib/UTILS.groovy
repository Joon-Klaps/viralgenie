// Shared helpers for the pipeline-level nf-tests.
// Taken from nf-core/sarek/tests/lib/UTILS.groovy and modified for nf-core/viralmetagenome.
// A test file becomes a list of scenarios plus one line:
//
//     def scenarios = [ [ name: "-profile test", scope: "files", bam: true ] ]
//     scenarios.each { scenario -> test(scenario.name, UTILS.getTest(scenario)) }
//
// Scenario keys
// -------------
//   name          required. The nf-test test name.
//   scope         "files" | "multiqc" | "versions". Default "multiqc".
//                   files    - whole outdir: stable names + stable contents
//                   multiqc  - the multiqc/ subtree only
//                   versions - task count + versions.yml only (params matrices)
//   bam           true to add per-BAM samtools statistics (implies scope "files").
//   vcf           true to add per-VCF variants md5 (implies scope "files").
//   params        map merged into the `when { params { } }` block.
//   overview      true to snapshot the overview-tables CSV column names/values.
//                 A table the run did not produce (e.g. the contigs table under
//                 `--skip_assembly`) is snapshotted as "No <filename>" instead.
//   sort_overview true to sort overview column names (test_group did this).
//   stub          true to run with -stub.
//   failure       true if the run is expected to fail.
//   tag           extra nf-test tag.
//   ignoreFiles   extra glob(s) excluded from the stable-content assertion.
//   logs          true to snapshot the run's WARN + ERROR lines, or a list of
//                 levels to narrow it, e.g. logs: ["ERROR"]. Always on for
//                 `failure` scenarios.
//   snapshot_ignore extra substrings dropped from the captured log lines, on top
//                 of SNAPSHOT_IGNORE.

class UTILS {

    // Log lines that legitimately vary between runs, machines and Nextflow
    // versions. Snapshotting them produces failures that say nothing about the
    // pipeline.
    public static def SNAPSHOT_IGNORE = [
        "Access to undefined parameter",
        "Creating env using",
        "Downloading plugin",
        "Got an interrupted exception while taking agent result",
        "Pulling Singularity image",
        "publishDir path contains a variable with a null value",
        "Staging foreign file",
        "Unable to resume cached task",
        "Unable to stage foreign file",
    ]

    public static def VERSIONS_YML = "pipeline_info/nf_core_viralmetagenome_software_mqc_versions.yml"

    // Builds the ordered list handed to snapshot(). Keeping the order fixed here
    // is what makes snapshots comparable across test files.
    public static def getAssertions = { Map args ->
        def outdir = args.outdir
        def scenario = args.scenario ?: [:]
        def workflow = args.workflow

        def scope = scenario.scope ?: "multiqc"
        def ignoreFiles = scenario.ignoreFiles ? [scenario.ignoreFiles].flatten() : []
        def assertion = []

        // A failing run has no meaningful task count or versions file.
        if (!scenario.failure) {
            assertion.add(workflow.trace.succeeded().size())
            assertion.add(removeNextflowVersion("${outdir}/${VERSIONS_YML}"))
        }

        if (scope == "files") {
            assertion.add(
                getAllFilesFromDir(
                    outdir,
                    relative: true,
                    includeDir: true,
                    ignore: ['pipeline_info/*.{html,json,txt}'] + ignoreFiles
                )
            )

            // Stub runs produce empty placeholder files; hashing them is noise.
            if (!scenario.stub) {
                assertion.add(
                    getAllFilesFromDir(outdir, ignoreFile: 'tests/.nftignore', ignore: ignoreFiles)
                )

                if (scenario.bam) {
                    def bam_files = getAllFilesFromDir(outdir, include: ['**/*.bam'], ignore: ignoreFiles)
                    assertion.add(
                        bam_files.isEmpty()
                            ? 'No BAM files'
                            : bam_files.collect { file -> [file.getName(), bam(file.toString()).getStatistics()] }
                    )
                }

                if (scenario.vcf) {
                    def vcf_files = getAllFilesFromDir(outdir, include: ['**/*.vcf.gz'], ignore: ignoreFiles)
                    assertion.add(
                        vcf_files.isEmpty()
                            ? 'No VCF files'
                            : vcf_files.collect { file -> [file.getName(), path(file.toString()).vcf.getVariantsMD5()] }
                    )
                }
            }
        }
        else if (scope == "multiqc") {
            assertion.add(
                getAllFilesFromDir(outdir, include: ['multiqc/**'], relative: true, includeDir: false)
            )
            if (!scenario.stub) {
                assertion.add(
                    getAllFilesFromDir(
                        outdir,
                        include: ['multiqc/**'],
                        ignoreFile: 'tests/.nftignore',
                        ignore: ignoreFiles
                    )
                )
            }
        }
        // scope == "versions" adds nothing beyond task count + versions.yml.

        // The overview tables are a special case: we snapshot the column names and the values of the sample/index column, but not the other columns.
        if (scenario.overview) {
            ["samples_overview.tsv": "sample", "contigs_overview-with-iterations.tsv": "index"].each { name, idColumn ->
                def tsv = new File("${outdir}/overview-tables/${name}")
                if (!tsv.exists()) {
                    assertion.add("No ${name}" as String)
                    return
                }
                def table = path(tsv.toString()).csv(sep: "\t")
                assertion.add(scenario.sort_overview ? table.columnNames.sort() : table.columnNames)
                assertion.add(table.columns[idColumn].sort())
            }
        }

        // A `failure` scenario otherwise asserts only `workflow.failed`, which any
        // breakage satisfies -- a typo in a module would pass exactly as well as
        // the error the test was written for. Snapshotting the message is what
        // pins it to the *expected* failure, so logs are always on there.
        // Default on for `failure`; an explicit `logs: false` still opts back out.
        def wantLogs = scenario.containsKey('logs') ? scenario.logs : scenario.failure
        def logLevels = wantLogs instanceof List ? wantLogs : (wantLogs ? ["WARN", "ERROR"] : null)

        if (logLevels) {
            // filterNextflowOutput strips ANSI codes, timestamps, task and work-dir
            // hashes, run names, absolute paths and container-engine names, which is
            // what makes these lines stable enough to snapshot at all.
            def std = workflow.stderr + workflow.stdout
            assertion.add(
                filterNextflowOutput(
                    std,
                    include: logLevels,
                    ignore: SNAPSHOT_IGNORE + (scenario.snapshot_ignore ?: [])
                ) ?: "No ${logLevels.join('/')} messages"
            )
        }

        return assertion
    }

    // Returns the body of an nf-test `test(...)` block for one scenario.
    public static def getTest = { scenario ->
        return {
            tag "pipeline"
            if (scenario.tag) {
                tag scenario.tag
            }
            if (scenario.failure) {
                tag "failure"
            }
            if (scenario.stub) {
                options "-stub"
            }

            when {
                params {
                    outdir = "${outputDir}"
                    (scenario.params ?: [:]).each { key, value ->
                        delegate."$key" = value
                    }
                }
            }

            then {
                if (scenario.failure) {
                    assert workflow.failed
                }
                else {
                    assert workflow.success
                }

                assertAll(
                    {
                        assert snapshot(
                            *UTILS.getAssertions(
                                outdir: params.outdir,
                                scenario: scenario,
                                workflow: workflow
                            )
                        ).match()
                    }
                )
            }

            cleanup {
                // getenv returns a String, so the bare truth check that upstream uses fires on
                // NFT_CLEANUP=false as well. Treat the usual falsy spellings as "off".
                def nft_cleanup = System.getenv('NFT_CLEANUP')
                if (nft_cleanup && !(nft_cleanup.toLowerCase() in ['false', '0', 'no', 'off'])) {
                    println ""
                    println "CLEANUP"
                    println "Set NFT_CLEANUP to false to disable."
                    println "The following folders will be deleted:"
                    println "- ${workDir}"

                    new File("${workDir}").deleteDir()
                }
            }
        }
    }
}
