

rule trimming:
	input:
		"input/{sample}_1.fastq.gz",
		"input/{sample}_2.fastq.gz",
	output:
		"trim/{sample}_1_trim.fastq.gz",
		"trim/{sample}_1U_trim.fastq.gz",
		"trim/{sample}_2_trim.fastq.gz",
		"trim/{sample}_2U_trim.fastq.gz"
	threads: 1

	shell:
		"trimmomatic PE -threads {threads} {input} {output}  ILLUMINACLIP:/apps/trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10  SLIDINGWINDOW:4:25 MINLEN:20"

rule map_inital:
	input:
		FR="trim/{sample}_1_trim.fastq.gz",
		RR="trim/{sample}_2_trim.fastq.gz",
		mm_genome_idx="/nlustre/users/werner/UKZN/KRISP/Listeria/NC003210"
	output:
	        BAM="bam.initial/{sample}.bam"
	threads: 1
	shell:
		"minimap2 -t {threads} -ax sr {input.mm_genome_idx} -1 {input.FR} -2 {input.RR} | samtools view -f 1 -F 12 -bu | samtools sort -n -l 1 -o {output}; samtools index {output}"

rule kraken_id:
	input:
		bam_initial="bam.initial/{sample}.bam"
	output:
		kraken_classify="kraken.listeria/{sample}_filter.txt.gz",
		kraken_report="kraken.report/{sample}_report.txt"
	threads: 2
	shell:
		"samtools fastq {input} | kraken2 --db /apps/kraken2/standard --memory-mapping /dev/fd/0 --threads {threads} --report {output.kraken_report} | grep 'Listeria monocytogenes' | cut -f 2 | sed 's/\/[0-9]$//g' | gzip -c  > {output.kraken_classify}"

rule filter_bam:
	input:
		bam_initial="bam.initial/{sample}.bam",
		kraken_classify="kraken.listeria/{sample}_filter.txt.gz"
	output:
		bam_filtered="bam.filtered/{sample}.bam"
	threads: 1
	shell:
		"picard filterSamReads O={output.bam_filtered} I={input.bam_initial} FILTER=includeReadList READ_LIST_FILE={input.kraken_classify}"

rule spades:
	input:
		bam_filtered="bam.filtered/{sample}.bam"
	output:
		assembly_dir="assemblies/{sample}"
	threads: 4
	shell:
		"samtools fastq  {input} | gzip -1 -c - > test.fq.gz ; spades.py --12 test.fq.gz --careful -t {threads} -k 55,77,97,127 -o {output} --tmp-dir /tmp"

rule prokka:












