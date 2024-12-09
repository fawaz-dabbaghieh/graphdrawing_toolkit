import os
import sys
import gzip
from collections import defaultdict
import pdb


def process_info(info_dict, data):
	sv_type = information['SVTYPE']
	if sv_type == "SNV":
		sv_len = 1
	else:
		if "SVLEN" in information:
			sv_len = abs(int(information['SVLEN']))
		else:
			if 'END' in information:
				try:
					sv_len = abs(information['start1'] - information['END'])
				except:
					pbd.set_trace()

	if sv_type == "DEL" and sv_len < 50:
		data[chromosome]["DEL_s50"] += 1
		data[chromosome][sv_type] += 1
	elif sv_type == "DEL" and sv_len >= 50:
		data[chromosome]["DEL_b50"] += 1
		data[chromosome][sv_type] += 1
	elif sv_type == "INS" and sv_len < 50:
		data[chromosome]["INS_s50"] += 1
		data[chromosome][sv_type] += 1
	elif sv_type == "INS" and sv_len >= 50:
		data[chromosome]["INS_b50"] += 1
		data[chromosome][sv_type] += 1
	else:   
		data[chromosome][sv_type] += 1


def mix_tables(imp_data, data):
	for chrom in data:
		if chrom not in imp_data:
			imp_data[chrom] = data[chrom]
		else:
			for sv_type in data[chrom]:
				if sv_type not in imp_data[chrom]:
					imp_data[chrom][sv_type] = data[chrom][sv_type]
				else:
					data[chrom][sv_type] += imp_data[chrom][sv_type]
	

if len(sys.argv) < 3:
	print("You need to give a prefix and the input VCF file")
	sys.exit()

prefix = sys.argv[1]
in_vcf = sys.argv[2]
data = defaultdict(lambda: defaultdict(int))
imp_data = defaultdict(lambda: defaultdict(int))

if in_vcf.endswith(".gz"):
	open_file = gzip.open(in_vcf, "rt")
else:
	open_file = open(in_vcf, "r")


for l in open_file:
	if l.startswith("##"):
		continue
	elif l.startswith("#"):
		l = l.strip().split()
		chrom = l.index("#CHROM")
		info = l.index("INFO")
		# example of info column from PAV output
		# ID=chr1-131-DEL-1;SVTYPE=DEL;SVLEN=-1;TIG_REGION=ptg000217l_length_870075_coverage_23:865050-865050;QUERY_STRAND=-;HOM_REF=0,1;HOM_TIG=0,1
		# as other callers might have different info in this column, I just need to check for SVTYPE and SVLEN
		# i need the second and 3rd parts here, the type and len, with len as absolute value to get rid of the negative
	else:
		l = l.strip().split()
		chromosome = l[chrom]
		information = dict()
		information['start1'] = int(l[1])
		for c in l[info].split(";"):
			c = c.split("=")
			if len(c) == 2:
				if c[0] == "SVLEN":
					information[c[0]] = abs(int(c[1]))
				if c[0] == "END":
					information[c[0]] = abs(int(c[1]))
				else:
					information[c[0]] = c[1]

		if l[info].startswith("IMPRECISE"):
			process_info(information, imp_data)
		else:
			process_info(information, data)

# pdb.set_trace()
# output precise stats
to_keep = [0,""]
delim = "\t"
for chrom in data:
	if len(data[chrom]) > to_keep[0]:
		to_keep = [len(data[chrom]), f"CHROM{delim}{delim.join(list(data[chrom].keys()))}\n"]

with open(f"{prefix}_precise_stats.tsv", "w") as outfile:
#	outfile.write(f"CHROM{delim}{delim.join(list(data['chr1'].keys()))}\n")
	outfile.write(to_keep[1])
	for chrm in data.keys():
		outfile.write(f"{chrm}{delim}{delim.join([str(x) for x in data[chrm].values()])}\n")


# output imprecise stats
with open(f"{prefix}_imprecise_stats.tsv", "w") as outfile:
#	outfile.write(f"CHROM{delim}{delim.join(list(imp_data['chr1'].keys()))}\n")
	outfile.write(to_keep[1])
	for chrm in imp_data.keys():
		outfile.write(f"{chrm}{delim}{delim.join([str(x) for x in imp_data[chrm].values()])}\n")

# output total as well
totals = defaultdict(int)
with open(f"{prefix}_totals.txt", "w") as outfile:
	for chrm in data.keys():
		for key in data[chrm]:
			totals[key] += data[chrm][key]
	for key in totals:
		outfile.write(f"{key}{delim}{totals[key]}\n")
