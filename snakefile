
rule fastq_to_ubam
	input: 
		"/scratch/mleukam/dlbcl/{sample}"
	output:
		"/scratch/mleukam/uBAMs/{sample}_unaligned.bam"
	output:
	shell:
		"java -Xmx4G -jar ${PICARD} FastqToSam F1={sample}_1.fq F2={sample}_2.fq O={output} SM={sample}"
