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

(sample_id,samples,sample) = sample_id_ch.into(3)


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

process scrublet {

  publishDir (
        path: "${params.outdir}/scrubletdir/",
        mode: 'copy',
        overwrite: 'true',
  )	 
    input:
    each samples 
    output: 
    
    path("${samples}_scrublet.{logic,score}") into scrublet_out
    
    script:

    """
	python3 ${baseDir}/scripts/scrublet_multi.py ${params.cellrangers_outs_dir} ${params.scrublet_SUFFIX} ${samples}
    """

}

process add_meta {

  publishDir (
        path: "${params.outdir}",
        mode: 'copy',
        overwrite: 'true',
  )	
        
    input:
    each sample 
    output: 
    
    	path("${sample}.h5") into meta_added
    
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
    
    file(in_h5) from meta_added.collect()
    file(scrub_dir) from scrublet_out.collect()
    output: 
    	
    script:
    """
    Rscript ${baseDir}/scripts/qcreporter/qc_batch_summary.r \
    	-e  ${params.experiment_id} \
    	-m  'scrna' \
    	-i  ${params.qc_in_dir} \
    	-z  ${params.cellrangers_outs_dir} \
    	-u  ${params.outdir}/scrubletdir/" \
    	-f  ${params.refdir} \
    	-k  ${params.in_key}   \
    	-d  ${params.qc_output} \
    	-o  ${params.qc_output}/${params.experiment_id}_rnaseq_sample_report.html \
	-c  ${params.percent_mito} \
	-j  ${params.resolution}
  """
}

