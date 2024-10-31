#!/usr/bin/env nextflow

samples = '*.bam*'
ensamples = '*.vcf*'
//samples = 'LP6005443-DNA_D01*bam*'

process insurveyor_calling {
	cpus 8
	maxForks 4
	module 'INSurVeyor'
	memory '40GB'
	publishDir "insurveyor_vcf"
	input:
		tuple val(core), path(f)
	output:
		path("${core}.insurveyor.vcf.gz")
		path("${core}.insurveyor.vcf.gz.tbi")
	script:
	out="${core}.insurveyor.vcf.gz"
	"""
	set -euxo pipefail
	insurveyor.py --threads 8 ${f[0]} ./ ${params.genome}
	mv out.pass.vcf.gz ${out}
	tabix -p vcf ${out}
	"""	
 }

process svaba_calling {
	errorStrategy 'finish'
	cpus 24
	maxForks 4
	module 'svaba'
	memory '40GB'
	publishDir "svaba_vcf"
	input:
		tuple val(core), path(f)
	output:
  		path("${core}.svaba.vcf.gz")
		path("${core}.svaba.vcf.gz.tbi")  
	script:
	out1="${core}.svaba.sv.vcf.gz"
	out2="${core}.svaba.indel.vcf.gz"
	out3="${core}.svaba.vcf.gz"
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

process fix_files{
	publishDir "svaba_fixed_vcf"
	input:
		tuple val (core), path(fs)
	output:
  		path(result)
		path("${result}.tbi")  
	script:
	result="${core}.svaba.fixed.vcf.gz"
	
	"""
	set -euxo pipefail
	zcat "${core}.svaba.vcf.gz" | sd "ID=GQ,Number=1,Type=String" "ID=GQ,Number=1,Type=Integer" | sed 's/ID=PL,Number=[^,]*/ID=PL,Number=G/' | bgzip > ${result}
	tabix -p vcf ${result}	
	"""
}

//	zcat "${core}.svaba.vcf.gz" | sd "ID=GQ,Number=1,Type=String" "ID=GQ,Number=1,Type=Float" | bgzip > ${result}


process truvari_ensamble {
	cpus 2
	maxForks 4
	module 'truvari'
	memory '20GB'
	publishDir "truvari_vcf"
	input:
                tuple val(core), path(files)
	output:
		path("${core}.truvari.vcf.gz")
		path("${core}.truvari.vcf.gz.tbi")
	script:
	out="${core}.vcf.gz"
	"""
	set -euxo pipefail
	bcftools merge -m none ${files[0]} ${files[2]} ${files[4]} ${files[6]} ${files[8]} -Oz -o ${core}.merge.vcf.gz --force-samples
	tabix -p vcf ${core}.merge.vcf.gz 
	truvari collapse -i ${core}.merge.vcf.gz -o ${core}.truvari.vcf.gz
	tabix -p vcf ${core}.truvari.vcf.gz 
	"""	
 }
//truvari bench --reference ${params.genome} --base ${files[0]} --comp ${files[2]} --comp ${files[4]} --comp ${files[6]} --comp ${files[8]} -o truBench 

process survivor_ensamble {
	cpus 2
	maxForks 4
	module 'SURVIVOR/1.0.7'
	memory '20GB'
	publishDir "survivor_vcf"
	input:
                tuple val(core), path(files)
	output:
		path("${core}.survivor.vcf.gz")
		path("${core}.survivor.vcf.gz.tbi")
	script:
	out="${core}.vcf.gz"
	"""
	set -euxo pipefail
	echo "${files[0]}" >> file.txt
	echo "${files[2]}" >> file.txt	
	echo "${files[4]}" >> file.txt
	echo "${files[6]}" >> file.txt
	echo "${files[8]}" >> file.txt

	SURVIVOR merge file.txt 1000 1 1 -1 -1 -1 ${core}.survivor.vcf.gz 
	tabix -p vcf ${core}.survivor.vcf.gz 

	"""	
 }





workflow {
	input = Channel.fromFilePairs(params.input+samples)
   		{ fname -> fname.simpleName.replaceAll(".md.*", "")}
	input2 = Channel.fromPath("/home/saleri/project_cnv_vr/scriptNF_running/svaba_vcf/*.vcf.gz").map {f -> 
		def core= f.getBaseName().toString().split("\\.")[0] 
		return tuple(core,f)}
	ensinput = Channel.fromFilePairs(params.directories.collect{it+ensamples}, size: 10) {fname -> fname.simpleName.replaceAll(".vcf", "") }
	

	main:
		//insurveyor_calling(input)
		//svaba_calling(input)
		//fix_files(input2)
		survivor_ensamble(ensinput)
		//truvari_ensamble(ensinput)		
		//println(params.directories.collect{it+ensamples})
		//ensinput.view()


}
                                       
                                              
