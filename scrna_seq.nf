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

(sample) = sample_id_ch.into(1)


println """\
         RNA Seq - N F   P I P E L I N E
         ===================================
         Experiment                : ${params.experiment_id}
         Samplesheet        	   : ${params.in_key}
         CellBender Directory : ${params.cellbender_dir}
         QC Report input directory : ${params.qc_in_dir}
         QC Report Output directory: ${params.qc_output}
         """
         .stripIndent()

process QC_Summary {

  publishDir (
        path: "${params.outdir}",
        mode: 'copy',
        overwrite: 'true',
  )	
        
    input:
    each sample 
    output: 
    
    	path("${sample}.h5") into qc_summary
    
    script:

  """
    Rscript ${baseDir}/scripts/rna_seq_pipeline_bwh/tenx_metadata_rna_adder.r \
    -i ${params.cellrangers_outs_dir}/${sample}/outs/filtered_feature_bc_matrix.h5   \
    -l ${params.cellrangers_outs_dir}/${sample}/outs/molecule_info.h5 \
    -s ${params.cellrangers_outs_dir}/${sample}/outs/metrics_summary.csv \
    -k ${params.in_key} \
    -j ${sample} \
  """

}

process QC_Report {

    publishDir params.qc_output 			// output dir
	
    input:
    
    file(qc_sum) from qc_summary.collect()
    output: 
    	
    script:
    """
    Rscript ${baseDir}/scripts/qcreporter/qc_batch_summary.r \
    	-e  ${params.experiment_id} \
    	-m  "snrna" \
    	-i  ${params.qc_in_dir} \
    	-z  ${params.cellbender_dir} \
    	-f  ${params.refdir} \
    	-k  ${params.in_key}   \
    	-d  ${params.qc_output} \
    	-o  ${params.qc_output}/${params.experiment_id}_rnaseq_sample_report.html \
	    -l  "Homo Sapiens" \
      -a  ${params.percent_ribo} \
      -j  ${params.resolution} \
      -b  ${params.filter_MALAT} \
      -c  ${params.percent_mito} \
      -u  ${params.filter_MITO} \
      -q  ${params.filter_RIBO} \
      -s  4
  """
}
