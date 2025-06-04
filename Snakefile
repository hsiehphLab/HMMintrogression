import os, glob

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
shell.prefix("source %s/env.cfg; set -eo pipefail; " % SNAKEMAKE_DIR)

configfile:'config.yaml'

if not os.path.exists("log"):
	os.makedirs("log")

## parameters
callableSeq = config["callableSeq"]
REFfa = config["REFfa"]
AAfa = config["AAfa"]
INITIAL = config["INITIAL"]

dict_chrom_gmap = {}
with open(config["f_list_gmap"]) as fin:
	for line in fin:
		l = line.strip()
		p_chrom = re.compile("(chr[0-9]+)")
		result = p_chrom.search(l)
		try:
			dict_chrom_gmap[result.group(1)] = l 
		except AttributeError:
			continue

dict_chrom_outgrp = {}
with open(config["f_list_outgrpVCF"]) as fin:
	for line in fin:
		l = line.strip()
		p_chrom = re.compile("(chr[0-9]+)")
		result = p_chrom.search(l)
		try:
			dict_chrom_outgrp[result.group(1)] = l 
		except AttributeError:
			continue

dict_chrom_ingrp = {}
with open(config["f_list_ingrpVCF"]) as fin:
	for line in fin:
		l = line.strip()
		p_chrom = re.compile("(chr[0-9]+)")
		result = p_chrom.search(l)
		try:
			dict_chrom_ingrp[result.group(1)] = l 
		except AttributeError:
			continue

dict_pop_ind_ingrp = {}
for pop in config["HumanPopID"]:
	with open(config["HumanPopID"][pop]) as fin:
		if pop not in dict_pop_ind_ingrp:
			dict_pop_ind_ingrp[pop] = []
		for line in fin:
			ind = line.strip()
			dict_pop_ind_ingrp[pop].append(ind)
	

outgroup_json = config["outgroup_JSON"]

def _get_ingroup_json(wildcards):
	return config["JSON"][wildcards.POP]


def _get_pop_sampleID(wildcards):
	return config["HumanPopID"][wildcards.POP]

def _get_OUTGROUP_sampleID(wildcards):
	return config["OutgroupPopID"][wildcards.OUTGROUP]

def _get_ARC_sampleID(wildcards):
	return config["ARCPopID"][wildcards.ARC]

def largemem_attempt(wildcards, attempt):
	return attempt * 40

def largemem_attempt2(wildcards, attempt):
	return attempt * 40

def mem_attempt(wildcards, attempt):
	return attempt * 10

wildcard_constraints:
	POP="[a-zA-Z]+",
	OUTGROUP="[a-zA-Z]+",
	ARC="[a-zA-Z]+",
	CHROM="chr[0-9]+",
	IND = "[a-zA-Z0-9_]+"


def _expand_HMMoutByIND(wildcards):
	l_out = []
	pat = "HMMoutput/{{POP}}_vs_{{OUTGROUP}}_{{ARC}}/{{POP}}_{IND}_vs_{{OUTGROUP}}_{{ARC}}.{CHROM}.hmm.bed"
	for chrom in dict_chrom_ingrp:
		for ind in dict_pop_ind_ingrp[wildcards.POP]:
			l_out.append(pat.format(CHROM=chrom, IND=ind))
	return(l_out) 


rule all:
	input: expand("HMMoutput_all/{POP}_vs_{OUTGROUP}_{ARC}.hmm.bed.gz", POP=config["HumanPopID"].keys(), OUTGROUP=config["OutgroupPopID"], ARC=config["ARCPopID"].keys())

rule bcf_outgrp:
	input:
		listVCF_outgrp = config["f_list_outgrpVCF"]
	output:
		bcf = "hmmix/vcf.outgrp.bcf",
		idx = "hmmix/vcf.outgrp.bcf.csi"
	threads: 4
	resources:
		mem=mem_attempt,
		hrs=72,
	shell:
		"""
			bcftools concat --threads {threads} -f {input.listVCF_outgrp} -Ov | bcftools view -l 1 -O b - > {output.bcf} && bcftools index {output.bcf}
		"""

							
rule hmmix_outgroup1:
	input:
		bcf = "hmmix/vcf.outgrp.bcf",
		idx = "hmmix/vcf.outgrp.bcf.csi",
		AA = AAfa,
		REF = REFfa,
		callableSeq = callableSeq,
		json = outgroup_json 
	output:
		outgroup = "hmmix/outgroup/outgroup.txt",
	threads: 1
	resources:
		mem=largemem_attempt,
		hrs=72,
	shell:
		"""
			hmmix create_outgroup -ind={input.json} -vcf={input.bcf} -weights={input.callableSeq} -ancestral={input.AA} -refgenome={input.REF} -out={output.outgroup} 

		"""

rule hmmix_mutationRate:
	input:
		outgroup = "hmmix/outgroup/outgroup.txt",
		callableSeq = callableSeq,
	output:
		mutationRate = "hmmix/mutationRate/outgroup.mutationRate.bed"
	threads: 1
	resources:
		mem=mem_attempt,
		hrs=72,
	shell:
		"""
			hmmix mutation_rate -outgroup={input.outgroup} -weights={input.callableSeq} -window_size=1000000 -out={output.mutationRate}

		"""

rule bcf_ingrp:
	input:
		listVCF_ingrp = config["f_list_ingrpVCF"]
	output:
		bcf = "hmmix/ingroup/vcf.ingrp.{POP}.{IND}.bcf",
		idx = "hmmix/ingroup/vcf.ingrp.{POP}.{IND}.bcf.csi"
	threads: 4
	resources:
		mem=mem_attempt,
		hrs=72,
	shell:
		"""
			bcftools concat --threads {threads} -f {input.listVCF_ingrp} -Ov | bcftools view -l 1 -O b - > {output.bcf} && bcftools index {output.bcf}
		"""


rule hmmix_ingroup:
	input:
		bcf = "hmmix/ingroup/vcf.ingrp.{POP}.{IND}.bcf",
		AA = AAfa,
		callableSeq = callableSeq,
		json = _get_ingroup_json, 
		outgroup = "hmmix/outgroup/outgroup.txt",
	output:
		obs = "hmmix/ingroup/obs.{POP}.{IND}.txt"
	threads: 1
	resources:
		mem=largemem_attempt,
		hrs=72,
	shell:
		"""
			x=`basename {output.obs} .{wildcards.POP}.{wildcards.IND}.txt`; hmmix create_ingroup  -ind={wildcards.IND} -vcf={input.bcf} -weights={input.callableSeq} -outgroup={input.outgroup} -ancestral={input.AA} -out=hmmix/ingroup/${{x}}.{wildcards.POP} 

		"""


rule hmmix_train:
	input:
		obs = "hmmix/ingroup/obs.{POP}.{IND}.txt",
		callableSeq = callableSeq,
		mutationRate = "hmmix/mutationRate/outgroup.mutationRate.bed"
	output:
		train = "hmmix/train/trained.{POP}.{IND}.json"
	threads: 1
	resources:
		mem=largemem_attempt,
		hrs=72,
	shell:
		"""
			hmmix train -param={INITIAL} -obs={input.obs} -haploid -weights={input.callableSeq} -mutrates={input.mutationRate} -out={output.train}
		"""


rule VCF2HMMinput:
	input: 
			target = lambda wc: dict_chrom_ingrp[wc.CHROM],
			outgroup = lambda wc: dict_chrom_outgrp[wc.CHROM], 
			outgroupID = _get_OUTGROUP_sampleID,
			arcID = _get_ARC_sampleID,
	output: 
			"HMMinput/{POP}_vs_{OUTGROUP}_{ARC}/{POP}_{IND}_vs_{OUTGROUP}_{ARC}.{CHROM}.bgl.hardmask.wAA.vcf.hmmInput.txt"
	threads: 1
	resources:
		mem=4,
		hrs=72,
	run:
		if wildcards.POP != wildcards.OUTGROUP:
			shell(
				"""
					python scripts/VCF2HMM_input.py  --target-vcf {input.target}  --outgrp-vcf {input.outgroup}  --file-targetIDs <(echo {wildcards.IND})  --file-outgrpIDs {input.outgroupID}  --file-arcIDs {input.arcID} > {output}

				""")
		else:
			shell(""" touch {output} """)


rule HMM:
	input: 
			hmmIN="HMMinput/{POP}_vs_{OUTGROUP}_{ARC}/{POP}_{IND}_vs_{OUTGROUP}_{ARC}.{CHROM}.bgl.hardmask.wAA.vcf.hmmInput.txt",
			train="hmmix/train/trained.{POP}.{IND}.json"
	output: 
			"HMMoutput/{POP}_vs_{OUTGROUP}_{ARC}/{POP}_{IND}_vs_{OUTGROUP}_{ARC}.{CHROM}.hmm.bed"
	threads: 1
	benchmark:
		'benchmark/hmm.{POP}_{IND}_vs_{OUTGROUP}_{ARC}.{CHROM}.benchmarkds'
	resources:
		mem=60,
		hrs=72,
	run:
		if wildcards.ARC == "NDL":
			t, m, ht = 2000, 0.02, 0.5
		elif wildcards.ARC == "DNS":
			t, m, ht = 1400, 0.04, 0.5

		if wildcards.POP != wildcards.OUTGROUP:
			shell(
			"""
				python scripts/run_hmm.py -i {input.hmmIN} -o {output} -chr %s -arc %s -t %s -m %s -ht %s -js {input.train} """ % (wildcards.CHROM, wildcards.ARC, t, m, ht))
		else:
			shell(""" touch {output} """)



rule ConcatBED:
	input:	_expand_HMMoutByIND 
	output: 
			"HMMoutput_all/{POP}_vs_{OUTGROUP}_{ARC}.hmm.bed.gz"
	threads: 1
	resources:
		mem=4,
		hrs=72,
	run:
		shell(""" cat {input} |sort -k1,1 -k2,2n -k5,5 | bgzip -c > {output} && tabix -p bed {output}""")

