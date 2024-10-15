#!/usr/bin/env nextflow

samples = '*.bam*'

process insurveyor_calling {
	cpus 8
	maxForks 4
	module 'INSurVeyor'
	memory '40GB'
	publishDir "insurveyor_vcf"
	input:
		tuple val(core), path(f)
	output:
		path("${core}.vcf.gz")
	script:
	out="${core}.vcf.gz"
	"""
	set -euxo pipefail
	insurveyor.py --threads 8 ${f[0]} ./ ${params.genome}
	mv out.pass.vcf.gz ${out}
	"""	
 }

workflow {
    input = Channel.fromFilePairs(params.input+"/"+samples)
   	  { fname -> fname.simpleName.replaceAll(".md.*", "")}

	 main:
		insurveyor_calling(input)

}
                                              
