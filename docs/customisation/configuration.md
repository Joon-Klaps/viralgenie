# Custom configuration of modules

Within viralgenie, all modules (tools, e.g., `FASTP`, `FASTQC`) can be run with specific arguments. The pipeline has a default configuration but this can be overwritten by supplying a custom configuration file. This file can be provided to viralgenie using the `-c` Nextflow option.

To see which specific arguments or variables are used for a module or tool, have a look at the [`modules.config` file](https://github.com/Joon-Klaps/viralgenie/blob/dev/conf/modules.config). Here the arguments of a module are specified as follows:

```groovy hl_lines="3-6 9-13"
withName: IVAR_CONSENSUS {
    ext.args = [
        '-t 0.75',          // frequency to call consensus: 0.75 just the majority rule
        '-q 20',            // minimum quality score of base
        '-m 10',            // minimum depth to call consensus
        '-n N'              // Character to print in regions with less coverage
    ].join(' ').trim()
    ext.args2 = [
        '--count-orphans',  // Do not skip anomalous read pairs in variant calling.
        '--max-depth 0',    // Maximum number of reads to start to consider at each location, 0 means no limit
        '--min-BQ 20',      // Minimum base quality
        '--no-BAQ',         // Disable probabilistic realignment for the computation of base alignment quality
        '-aa',              // Output absolutely all positions, including unused reference sequences
    ].join(' ').trim()
    ...
}
```

In this example, the `IVAR_CONSENSUS` module is configured with the arguments `-q 20 -m 10` for the tool [`ivar consensus`](https://andersen-lab.github.io/ivar/html/manualpage.html#autotoc_md19) and `--ignore-overlaps --count-orphans --max-depth 0 --no-BAQ --min-BQ 0` for [`samtools mpileup`](https://www.htslib.org/doc/samtools-mpileup.html) as iVar uses the output of `samtools mpileup` directly.

!!! Tip
    The `ext.args` and `ext.args2` are used to specify the arguments for the tool. If unsure which tools use which arguments (`ivar:ext.args` and `samtools:ext.args2`), have a look at the nextflow module file directly! For example, at [`modules/nf-core/ivar/consensus.nf`](https://github.com/Joon-Klaps/viralgenie/blob/dev/modules/nf-core/ivar/consensus/main.nf), "$args" and "$args2" are used to specify the arguments for the tools:
    ```groovy hl_lines="5 10"

    """
    samtools \\
        mpileup \\
        --reference $fasta \\
        $args2 \\               // can be modified with ext.args2
        $bam \\
        $mpileup \\
        | ivar \\
            consensus \\
            $args \\            // can be modified with ext.args
            -p $prefix
    ...
    """
    ```

In case we do want to modify the arguments of a module, we can do so by providing a custom configuration file. The easiest way to do this would be to copy a segment from the modules.config and modify the arguments. This way, none of the other configurations will get lost or modified. For example, setting the minimum depth to call consensus to 5 and the minimum quality score of base to 30 for the `IVAR_CONSENSUS` module:
```groovy title='custom.config' hl_lines="5-6 12"
process {
    withName: IVAR_CONSENSUS {
        ext.args = [
            '-t 0.75',
            '-q 30',            // changed
            '-m 5',             // changed
            '-n N'
        ].join(' ').trim()
        ext.args2 = [
            '--count-orphans',
            '--max-depth 0',
            '--min-BQ 30',      // changed
            '--no-BAQ',
            '-aa',
        ].join(' ').trim()
        ...
    }
}
```

!!! Warning
    Make sure you include the `process{}` section.

Next, supply the file to viralgenie using the `-c` Nextflow option:
```bash
nextflow run Joon-Klaps/viralgenie \
    -profile docker \
    -c custom.config \
    --input samplesheet.csv ...
```

!!! Tip
    This guide not entirely clear? Also have a look at the [nf-core guide for customizing tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments).
