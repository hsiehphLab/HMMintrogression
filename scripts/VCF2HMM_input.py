import sys, os, re, vcf
import argparse


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--target-vcf", nargs="?", required=True, metavar='TARGET_POP_VCF', help='population of interest for introgression tests' )
	parser.add_argument("--outgrp-vcf", nargs="?", required=True, metavar='OUTGROUP_POP_VCF', help='outgroup population')
	parser.add_argument("--file-targetIDs", nargs="?", required=True, metavar='FILE_TARGET_SAMPLEID', help='a file of sample IDs for the target population')
	parser.add_argument("--file-outgrpIDs", nargs="?", required=True, metavar='FILE_OUTGROUP_SAMPLEID', help='a file of sample IDs for the outgroup population')
	parser.add_argument("--file-arcIDs", nargs="?", required=True, metavar='FILE_ARC_SAMPLEID', help='a file of ARC sample IDs')
	args = parser.parse_args()

	target_vcfReader = vcf.Reader(open(args.target_vcf, "rb"))
	outgrp_vcfReader = vcf.Reader(open(args.outgrp_vcf, "rb"))


	l_target, l_outgrp, l_ARC = [], [], []

	with open(args.file_targetIDs) as fin:
		for line in fin:
			l_target.append(line.strip())

	l_targetVCF_sampleIDs = target_vcfReader.samples
	l_targetVCF_sampleIDs_avail = [ s for s in l_targetVCF_sampleIDs if s in l_target]

	l_haploID_target = [ str(s)+"_hap"+str(i) for s in l_targetVCF_sampleIDs_avail for i in range(1,3)]
	out_header = ["pos","f_outgroup","f_ARC"]
	out_header.extend(l_haploID_target)
	print("\t".join([str(x) for x in out_header]))

	with open(args.file_outgrpIDs) as fin:
		for line in fin:
			l_outgrp.append(line.strip())

	l_outgrpVCF_sampleIDs = outgrp_vcfReader.samples
	l_outgrpVCF_sampleIDs_avail = [ s for s in l_outgrpVCF_sampleIDs if s in l_outgrp]

	with open(args.file_arcIDs) as fin:
		for line in fin:
			l_ARC.append(line.strip())

	l_targetVCF_ARCsampleIDs_avail = [ s for s in l_targetVCF_sampleIDs if s in l_ARC]


	for record in target_vcfReader:
		CHROM = record.CHROM
		POS = record.POS
		ID = record.ID
		l_REF_ALT = [record.REF]
		l_REF_ALT.extend([str(allele) for allele in record.ALT])

		try:
			AA = record.INFO["AA"]
		except KeyError:
			continue

		try:
			AA_genocode = l_REF_ALT.index(AA)
		except ValueError:
			continue

		l_geno_ARC = [int(g) for s in l_targetVCF_ARCsampleIDs_avail for g in re.split("[| /]", record.genotype(s)["GT"]) if g != "."]
		if len(l_geno_ARC) == 0:
			continue

		f_derived_ARC = float(len(l_geno_ARC) - l_geno_ARC.count(AA_genocode)) / len(l_geno_ARC)

		# pull out record from the outgroup pop
		obj_record_outgrp = outgrp_vcfReader.fetch(CHROM, POS-1, POS)
		l_record_outgrp = [x for x in obj_record_outgrp]
		for r in l_record_outgrp:
			if r.ID != ID:
				continue
			else:
#		if len(l_record_outgrp) != 1:
#			print (l_record_outgrp)
#			exit("%s:%s has more than one record! EXIT..." % (CHROM, POS))

				l_geno_outgrp = [int(g) for s in l_outgrpVCF_sampleIDs_avail for g in re.split("[| /]", l_record_outgrp[0].genotype(s)["GT"]) if g != "."]
				if len(l_geno_outgrp) == 0:
					continue

				f_derived_outgrp = float(len(l_geno_outgrp) - l_geno_outgrp.count(AA_genocode)) / len(l_geno_outgrp)
				
				l_geno_target = [g for s in l_targetVCF_sampleIDs_avail for g in re.split("[| /]", record.genotype(s)["GT"])]
				if len(l_geno_target) == 0:
					continue
				l_derivedGeno_target = [ "." if g == "." else 0 if int(g) == AA_genocode else AA_genocode if int(g) == 0 else g for g in l_geno_target]

				l_out = [POS, "%.5f" % f_derived_outgrp, "%.5f" % f_derived_ARC]
				l_out.extend(l_derivedGeno_target)
				print("\t".join([str(x) for x in l_out]))

