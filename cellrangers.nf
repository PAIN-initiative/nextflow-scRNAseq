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

Channel
    .fromPath( params.samples_csv )
    .splitCsv( header: true, sep: ',' )
    .map { row ->  row.sample_id }
    .set { sample_id_ch }

(sample_id) = sample_id_ch.into(1)



/*
Channel
    .fromPath( params.samples_csv )
    .splitCsv( header: true, sep: ',' )
    .map { row ->  row.sample_id }
    .set { sample_id }
*/
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
  file("${sample_id}/*") into cellrangers_outs
  script: 
  """
  ${params.cellranger_software_path}/cellranger count --id=$sample_id \
                   --transcriptome=${ref} \
                   --fastqs=${params.fastq_path}/$sample_id \
                   --sample=$sample_id \
                   --localcores=${params.cores} \
                   --localmem=${params.mem}
  """

}
