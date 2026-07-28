#!/usr/bin/env nextflow

/*
Merge *SJ.out.tab files across all samples and filters out:
  - Non-canonical junctions (col 5 == 0)
  - Junctions mapping to the mitochondrial chromosome
  - Junctions supported by fewer than the minimum number of allowed unique reads (min_unique_reads)

Outputs: a single merged junctions table
*/

process MERGE_SJ {

      tag "merging splice junction files..."
      conda 'conda-forge::gawk=5.3.0'

      publishDir "${params.output}/STAR_junctions", mode: 'copy'

      input:
      path sj_files
      val min_unique_reads

      output:
      path "merged_junction.SJ.out.tab", emit: merged_sj

      script:
      """
      # Concatenate all SJ files, then filter:
      #   col 5 > 0   : canonical junction motif (GT/AG, GC/AG, AT/AC)
      #   col 1 != MT : exclude mitochondrial chromosome
      #   col 7 >= N  : minimum uniquely mapping reads across all samples

      cat ${sj_files} \\
          | awk '(\$5 > 0 && \$1 != "MT" && \$7 >= ${min_unique_reads})' \\
          | sort -k1,1 -k2,2n -k3,3n \\
          | uniq \\
          > merged_junction.SJ.out.tab

      echo "Total junction after filtering: \$(wc -l < merged_junction.SJ.out.tab)"
      """

}
