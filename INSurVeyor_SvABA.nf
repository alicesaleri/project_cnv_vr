#!/usr/bin/env nextflow

samples = '*.bam*'
samples = 'LP6005443-DNA_D01*bam*'

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

process svaba_calling {
	cpus 24
	maxForks 4
	module 'svaba'
	memory '40GB'
	publishDir "svaba_vcf"
	input:
		tuple val(core), path(f)
	output:
  		path("${core}.vcf.gz")
		path("${core}.vcf.gz.tbi")  
	script:
	out1="${core}.svaba.sv.vcf.gz"
	out2="${core}.svaba.indel.vcf.gz"
	out3="${core}.vcf.gz"
	"""
	set -euxo pipefail
	/usr/bin/time -o "resources.t" -f "%e %M" svaba run -p 24 -G ${params.genome} -t ${f[0]} -z 
	mv no_id.svaba.sv.vcf.gz ${out1}
	mv no_id.svaba.indel.vcf.gz ${out2}
	bcftools index -t ${out1}
	bcftools index -t ${out2}
	bcftools concat ${out1} ${out2} -a -D -Oz -o ${out3}
	tabix -p vcf ${out3} 
	"""	
 }


workflow {
    input = Channel.fromFilePairs(params.input+samples)
   	  { fname -> fname.simpleName.replaceAll(".md.*", "")}

	 main:
		//insurveyor_calling(input)
		svaba_calling(input)

}
                                              
