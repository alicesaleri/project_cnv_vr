#!/usr/bin/env nextflow

samples = '*.bam*'

process insurveyor_calling {
	cpus 2	
	maxForks 8
	module 'INSurVeyor'
	publishDir "insurveyor_vcf"
	input:
		tuple val(core), path(f)
	output:
		path("${core}.vcf.gz")
	script:
	out="${core}.vcf.gz"
	"""
	mkdir ${core}.vcf.gz
	chmod a+wrx ${core}.vcf.gz
	insurveyor.py --threads 4 ${f[0]} ./ ${params.genome}
	"""	
 }

workflow {
    input = Channel.fromFilePairs(params.input+"/"+samples)
   	  { fname -> fname.simpleName.replaceAll(".md.*", "")}

	 main:
		insurveyor_calling(input)
	        input.view()

}
                                              
