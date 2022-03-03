
rule all:
	input:
		"/gpfs/project/serramar/summed_altusccs_to_altus_filtered_unique_sequences_5kb.depth-coverage-cns.csv"

#align the sequencing reads to the unitigs
rule align:
	input: 
		nodes = "/gpfs/project/serramar/potato-hifi/altusHifi/altusGFA_270121/spud_altus_hifiasm.r_utg.nodes.fasta",
		sample1 = "/gpfs/project/projects/medbioinf/data/potato_data/altus_ccs/m54336U_201108_053012.ccs.fastq.gz", 
		sample2 = "/gpfs/project/projects/medbioinf/data/potato_data/altus_ccs/m54336U_201109_115807.ccs.fastq.gz", 
		sample3 = "/gpfs/project/projects/medbioinf/data/potato_data/altus_ccs/m54336U_201126_080917.ccs.fastq.gz" 			
	threads: 12	
	log:
		"/gpfs/project/serramar/altusccs_to_altusgfa.log"
	output:
		"/gpfs/project/serramar/altusccs_to_altusgfa.sam"
	shell:
		"/usr/bin/time -v /gpfs/project/serramar/software/minimap2-2.18_x64-linux/minimap2 -ax asm20 -t{threads} {input.nodes} {input.sample1} {input.sample2} {input.sample3} -o {output} |& tee {log}"
		
#sort sam file and convert to bam:		
rule sort_to_bam:		
	input:
		sam = "/gpfs/project/serramar/altusccs_to_altusgfa.sam"
	output:
		bam = "/gpfs/project/serramar/altusccs_to_altusgfa_sorted.bam"
	shell:
		"/usr/bin/time -v samtools view -@12 -bS {input.sam} | samtools sort -@12 - > {output.bam}; samtools index {output.bam};"
		
		
#filter bam file for MAPQ
rule filter:
	input:
		"/gpfs/project/serramar/altusccs_to_altusgfa_sorted.bam"
	threads: 4
	output:
		"/gpfs/project/serramar/altusccs_to_altusgfa_sorted_filtered.bam"
	shell:
		"/usr/bin/time -v samtools view -@{threads} -bq 60 {input} > {output}; samtools index -@{threads} {output};"

#compute sequencing depth:
rule compute_depth:
	input:
		bed = "/gpfs/project/serramar/spud_altus_hifiasm.r_utg.unique_sequences.bed",
		bams = "/gpfs/project/serramar/altusccs_to_altusgfa_sorted_filtered.bam",
	output:
		"/gpfs/project/serramar/altusccs_to_altus_filtered_unique_sequences.depth"	
	shell:
		"/usr/bin/time -v samtools depth -a -b {input.bed} -f {input.bams} -s -o {output}";

#sum up depths per node & window:
rule sum_depths:
	input:
		nodes = "/gpfs/project/serramar/spud_altus_hifiasm.r_utg.nodes",
		depth = "/gpfs/project/serramar/altusccs_to_altus_filtered_unique_sequences.depth"
	output:
		"/gpfs/project/serramar/summed_altusccs_to_altus_filtered_unique_sequences_5kb.depth"
	shell:
		"time python /gpfs/project/serramar/sum_node_cov_5kb.py {input.depth} {output} {input.nodes};"

#write copy number output file
rule plot_coverage:
	input:
		"/gpfs/project/serramar/summed_altusccs_to_altus_filtered_unique_sequences_5kb.depth"
	output:
		"/gpfs/project/serramar/summed_altusccs_to_altus_filtered_unique_sequences_5kb.depth-coverage-cns.csv" #how to add file endings here? Or give input as file path again and omit output if possible?
	shell:
		"python3 /home/rebecca/work/hifi-potato/wholegenome_kmercounts/write_dosages.py {input} {output}"		


		