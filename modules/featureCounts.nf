#!/usr/bin/env nextflow

// FeatureCounts
process READCOUNTS {

        tag "sample: $sample_name"

        publishDir "${params.output}/$sample_name/featureCounts", mode: 'copyNoFollow'

        input:
        tuple val(sample_name), path(bam_files)
        path annotation
        val strand

        output:
        path "*{exon,gene}.{stranded,unstranded}.counts.txt.summary", emit: FC_summary_files
        path "*{exon,gene}.{stranded,unstranded}.counts.txt*"

        script:

          if(params.prep == "single-end") {
              """
              echo 'Counting Reads on features ...stranded mode'
              featureCounts -T $task.cpus -t exon -a $annotation -s $strand -d 20 -D 4000 -o ${sample_name}.exon.stranded.counts.txt $bam_files

              echo 'Counting Reads on features ...unstranded mode'
              featureCounts -T $task.cpus -t exon -a $annotation -s 0 -d 20 -D 4000 -o ${sample_name}.exon.unstranded.counts.txt $bam_files

              echo 'Counting Reads on genes ...stranded mode'
              featureCounts -T $task.cpus -t gene -a $annotation -s $strand -d 20 -D 4000 -o ${sample_name}.gene.stranded.counts.txt $bam_files

              echo 'Counting Reads on genes ...unstranded mode'
              featureCounts -T $task.cpus -t gene -a $annotation -s 0 -d 20 -D 4000 -o ${sample_name}.gene.unstranded.counts.txt $bam_files
              """
          } else {
              """
              echo 'Counting Reads on features ...stranded mode'
              featureCounts -T $task.cpus -p --countReadPairs -B -t exon -a $annotation -s $strand -d 20 -D 4000 -o ${sample_name}.exon.stranded.counts.txt $bam_files

              echo 'Counting Reads on features ...unstranded mode'
              featureCounts -T $task.cpus -p --countReadPairs -B -t exon -a $annotation -s 0 -d 20 -D 4000 -o ${sample_name}.exon.unstranded.counts.txt $bam_files

              echo 'Counting Reads on genes ...stranded mode'
              featureCounts -T $task.cpus -p --countReadPairs -B -t gene -a $annotation -s $strand -d 20 -D 4000 -o ${sample_name}.gene.stranded.counts.txt $bam_files

              echo 'Counting Reads on genes ...unstranded mode'
              featureCounts -T $task.cpus -p --countReadPairs -B -t gene -a $annotation -s 0 -d 20 -D 4000 -o ${sample_name}.gene.unstranded.counts.txt $bam_files
              """

          }

}
