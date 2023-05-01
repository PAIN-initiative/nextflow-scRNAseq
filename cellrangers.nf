#!/usr/bin/ nextflow

nextflow.enable.dsl=1

/*
========================================================================================
   QC Report Generator Nextflow Workflow
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------
*/


// Replace this with the path to a directory containing raw fastq files
params.fastqs_dir = '/mnt/data0/projects/donglab/EWA_Ruifeng2023/data/fastqs/'
fastq_path = params.fastqs_dir


Channel
    .fromPath( params.samples_csv )
    .splitCsv( header: true, sep: ',' )
    .map { row ->  row.sample_id }
    .set { sample_id }

println """\
         RNA Seq - N F   P I P E L I N E
         ===================================
         Experiment       		     : ${params.experiment_id}
         Samplesheet        		   : ${params.in_key}
         CellRangersOuts Directory : ${params.cellrangers_outdir}
         QC Report input directory : ${params.qc_in_dir}
         QC Report Output directory: ${params.qc_output}
         """
         .stripIndent()

/* set ref directory e.g. (mouse,human,fly) */
ref = params.genomedir

process cellranger_count {

  publishDir (
  path: "${params.outdir}/cellrangersouts",
  mode: 'copy',
  overwrite: 'true',
    )

  input:
  each sample_id

  output:
  file("${sample_id}/*") into cellrangers_outs,cellranger,counts
  script: 
  """
  ${params.cellranger_software_path}/cellranger count --id=$sample_id \
                   --transcriptome=$ref \
                   --fastqs=$fastq_path/$sample_id \
                   --sample=$sample_id \
                   --localcores=30 \
                   --localmem=64
  """

}
