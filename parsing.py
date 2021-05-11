from Bio import SeqIO
from Bio.SeqIO import parse 
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq
from Bio import Entrez

import re
import io
import os.path
import sys


def progress(count, total, status=''):

	sys.stdout.write('\x1b[1A')
	sys.stdout.write('\x1b[2K') #retire ligne précédente

	bar_len = 60
	filled_len = int(round(bar_len * count / float(total)))

	percents = round(100.0 * count / float(total), 1)
	bar = '=' * (filled_len) + '>' + '.' * (bar_len - filled_len)

	sys.stdout.write('[%s] %s%s\t|\t%s\r' % (bar, percents, '%', status))
	sys.stdout.flush()

def write_seq(file, entity_id):
	old_stdout = sys.stdout
	new_stdout = io.StringIO()
	sys.stdout = new_stdout
	function_group = []
	with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=entity_id) as handle:
		gb_record = SeqIO.read(handle, "gb")
	for f in gb_record.features:
		if f.type == "CDS" or f.type == "centromere" or f.type == "intron" or f.type == "mobile_element" or f.type == "ncRNA" or f.type == "rRNA" or f.type == "telomere" or f.type == "tRNA" or f.type == "3'UTR" or f.type == "5'UTR":
			function_group.append(str(f.type))
			loc = str(f.location)
			loc_ = loc
			loc = loc.replace(":", " ")
			loc = loc.replace("[", "")
			loc = loc.replace("]", "")
			loc = loc.replace("(", "")
			loc = loc.replace(")", "")
			loc = loc.replace("-", "")
			loc = loc.replace("+", "")
			loc = loc.replace("join", "")
			loc = loc.replace("complement", "")
			loc = loc.replace("{", "")
			loc = loc.replace("}", "")
			loc = loc.replace(",", "")
			loc = loc.replace(">", "")
			loc = loc.replace("<", "")
			loc = loc.replace("order", "")
			coord = loc.split(' ')
			int_coord = []
			for num in coord:
				int_coord.append(int(num))
			if not "join" in loc_ and not "complement" in loc_:
				inf, sup = min(int_coord), max(int_coord)
				print("%s\t%d..%d" % (f.type, inf, sup))
				print(gb_record.seq[inf:sup] + '\n')
			if "complement" in loc_ and not "join" in loc_:
				inf, sup = min(int_coord), max(int_coord)
				print("%s\tcomplement(%d..%d)" % (f.type, inf, sup))
				print(gb_record.seq[inf:sup] + '\n')
			if "join" in loc_:
				join(int_coord, gb_record.seq, f, file)
	output = new_stdout.getvalue()
	print(output, file=file)
	sys.stdout = old_stdout
	return list(dict.fromkeys(function_group))


def create_file(index, split_string, entity_id, ids):
	# Vérifier si le dossier existe ou non : inutile si les lignes sont distinctes
	path = 'Results'+'/'+split_string[1]+'/'+split_string[2]+'/'+split_string[3]+'/'+split_string[0]
	function_group = ""
	if not os.path.isdir(path):
		os.makedirs(path)

	if not os.path.isfile(path + '/' + entity_id + '.txt'):
		file = open(path + '/' + entity_id + '.txt', "a")
	
		ids = ids + ',' + entity_id
		for type_ in write_seq(file, entity_id):
			function_group = function_group + '\t' + type_
		file.close()
		index.write(path + '/' + entity_id + '.txt' + function_group + '\n')

	return (path, ids)

def find_ids(file, entity_name):
	# iterate over lines, and print out line numbers which contain
	# the word of interest.
	for line in file:
		split_string = line.split("\t")
		if entity_name in split_string: # or word in line.split() to search for full words
			return split_string[1]

def init():

	if os.path.isfile('index.txt'):
		os.remove('index.txt')
	# Ouvrir le fichier en lecture seule
	overview = open('GENOME_REPORTS/overview.txt', "r")
	archaea_ids = open('GENOME_REPORTS/IDS/Archaea.ids', "r")
	bacteria_ids = open('GENOME_REPORTS/IDS/Bacteria.ids', "r")
	eukaryota_ids = open('GENOME_REPORTS/IDS/Eukaryota.ids', "r")
	viruses_ids = open('GENOME_REPORTS/IDS/Viruses.ids', "r")
	# utiliser readlines pour lire toutes les lignes du fichier
	# La variable "lignes" est une liste contenant toutes les lignes du fichier
	overview_lines = overview.readlines()
	archaea_lines = archaea_ids.readlines()
	bacteria_lines = bacteria_ids.readlines()
	eukaryota_lines = eukaryota_ids.readlines()
	viruses_lines = viruses_ids.readlines()
	index = open('index.txt', "a")
	ids_files = [archaea_lines, bacteria_lines, eukaryota_lines, viruses_lines]
	# fermez le fichier après avoir lu les lignes
	overview.close()
	archaea_ids.close()
	bacteria_ids.close()
	eukaryota_ids.close()
	viruses_ids.close()

	nb_lines = len(overview_lines)
	ids = ''
	count=0
	# Itérer sur les lignes sauf la première


	for num, line in enumerate(overview_lines[1:], 1):
		#print(line.strip())
		progress(num, nb_lines, status='Creating directories')

		split_string = line.split("\t")
		if split_string[1] == 'Archaea':
			entity_id = find_ids(ids_files[0], split_string[0])
			if not entity_id == None and 'NC' in entity_id:
				(path, ids) = create_file(index, split_string, entity_id, ids)
				count=count+1
		if split_string[1] == 'Bacteria':
			entity_id = find_ids(ids_files[1], split_string[0])
			if not entity_id == None and 'NC' in entity_id:
				(path, ids) = create_file(index, split_string, entity_id, ids)
				count=count+1
		if split_string[1] == 'Eukaryota':
			entity_id = find_ids(ids_files[2], split_string[0])
			if not entity_id == None and 'NC' in entity_id:
				(path, ids) = create_file(index, split_string, entity_id, ids)
				count=count+1
		if split_string[1] == 'Viruses':
			entity_id = find_ids(ids_files[3], split_string[0])
			if not entity_id == None and 'NC' in entity_id:
				(path, ids) = create_file(index, split_string, entity_id, ids)
				count=count+1

	index.close()

def join(coord, sequence, f, file=''):
	old_stdout = sys.stdout
	new_stdout = io.StringIO()
	sys.stdout = new_stdout
	l = len(coord)
	loc = str(f.location)
	loc = loc.replace(":", "..")
	loc = loc.replace("[", "")
	loc = loc.replace("]", "")
	loc = loc.replace("(", "")
	loc = loc.replace(")", "")
	loc = loc.replace("-", "")
	loc = loc.replace("+", "")
	loc = loc.replace("{", "(")
	loc = loc.replace("}", ")")
	loc = loc.replace(">", "")
	loc = loc.replace("<", "")
	loc = loc.replace("order", "")
	string = ""
	for i in range(round(l/2)):
		inf = coord[2*i]
		sup = coord[2*i+1]
		print("%s\t%s\t%d..%d" % (f.type, loc, inf, sup))
		print(sequence[inf:sup] + '\n')
	output = new_stdout.getvalue()
	print(output, file=file)
	sys.stdout = old_stdout


## --------------------------------------------------------------------------- ##

Entrez.email = "thomas18199@hotmail.fr"
init()



