#!/usr/bin/env nextflow

samples = '*.bam*'

process svaba_calling {
	cpus 4	
	maxForks 8
	module 'svaba'
	memory '40GB'
	publishDir "svaba_vcf"
	input:
		tuple val(core), path(f)
	output:
		path("${core}.vcf.gz")
	script:
	out="${core}.vcf.gz"
	"""
	set -euxo pipefail
	svaba run -p 4 -G ${params.genome} -t ${f[0]} -z ./
	mv out.pass.vcf.gz ${out}
	"""	
 }

workflow {
    input = Channel.fromFilePairs(params.input+"/"+samples)
   	  { fname -> fname.simpleName.replaceAll(".md.*", "")}

	 main:
		svaba_calling(input)

}
                                              
