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

	bar_len = 60
	filled_len = int(round(bar_len * count / float(total)))

	percents = round(100.0 * count / float(total), 1)
	bar = '=' * (filled_len) + '>' + '.' * (bar_len - filled_len)

	sys.stdout.write('[%s] %s%s\t|\t%s\r' % (bar, percents, '%', status))
	sys.stdout.flush()

def write_seq(file, gb_record):
	old_stdout = sys.stdout
	new_stdout = io.StringIO()
	sys.stdout = new_stdout
	function_group = []
	features = gb_record.features
	for f in features:
		if f.type == "CDS" or f.type == "centromere" or f.type == "intron" or f.type == "mobile_element" or f.type == "ncRNA" or f.type == "rRNA" or f.type == "telomere" or f.type == "tRNA" or f.type == "3'UTR" or f.type == "5'UTR":
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

def date_convert(date):
	date = date.replace("JAN", "01")
	date = date.replace("FEB", "02")
	date = date.replace("MAR", "03")
	date = date.replace("APR", "04")
	date = date.replace("MAY", "05")
	date = date.replace("JUN", "06")
	date = date.replace("JUL", "07")
	date = date.replace("AUG", "08")
	date = date.replace("SEP", "09")
	date = date.replace("OCT", "10")
	date = date.replace("NOV", "11")
	date = date.replace("DEC", "12")

	date = date.split("-")
	date = date[2] + date[1] + date[0]
	return date

def update(old_date, new_date):
	olddate = date_convert(old_date)
	newdate = date_convert(new_date)
	return olddate < newdate

def create_file(split_string, entity_id, ids, gb_record, path):
	# Vérifier si le dossier existe ou non : inutile si les lignes sont distinctes
	function_group = ""
	new_date = gb_record.annotations.get("date")
	if not os.path.isdir(path):
		os.makedirs(path)

	if not os.path.isfile(path + '/' + entity_id + '.txt'):
		date_file = open(path + '/date.dat', "a")
		date_file.write(new_date)
		date_file.close()
		file = open(path + '/' + entity_id + '.txt', "a")
	
		ids = ids + ',' + entity_id
		write_seq(file, gb_record)
		file.close()

	date_file = open(path + '/date.dat', "r")
	old_date = date_file.readlines()[0]
	date_file.close()

	if update(old_date, new_date):
		print("changing date")
		date_file = open(path + '/date.dat', "w")
		date_file.write(new_date)
		date_file.close()
		file = open(path + '/' + entity_id + '.txt', "w")
	
		ids = ids + ',' + entity_id
		write_seq(file, gb_record)
		file.close()


def filter(index, filtre, entity_id, path):

	with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=entity_id) as handle:
		gb_record = SeqIO.read(handle, "gb")

	if os.path.isfile(path + '/date.dat'):
		new_date = gb_record.annotations.get("date")
		date_file = open(path + '/date.dat', "r")
		old_date = date_file.readlines()[0]
		date_file.close()
		if not update(old_date, new_date):
			return (False, gb_record)

	if not os.path.isfile(path + '/date.dat'):

		new_date = gb_record.annotations.get("date")
		date_file = open(path + '/date.dat', "a")
		date_file.write(new_date)
		date_file.close()
		return (False, gb_record)

	function_group = []
	functiongroup = ''
	features = gb_record.features
	for f in features:
		if f.type == "CDS" or f.type == "centromere" or f.type == "intron" or f.type == "mobile_element" or f.type == "ncRNA" or f.type == "rRNA" or f.type == "telomere" or f.type == "tRNA" or f.type == "3'UTR" or f.type == "5'UTR":
			function_group.append(str(f.type))
	function_group = list(dict.fromkeys(function_group))
	for type_ in function_group:
			functiongroup = functiongroup + '\t' + type_

	index.write(path + '/' + entity_id + '.txt' + functiongroup + '\n')

	if filtre[0] == '':
		return (True, gb_record)
	#for f in gb_record.features:
	record = str(gb_record.features)
	for group in filtre:
		if not group in record:
			return (False, gb_record)
	return (True, gb_record)

def find_ids(file, entity_name):
	# iterate over lines, and print out line numbers which contain
	# the word of interest.
	for line in file:
		split_string = line.split("\t")
		if entity_name in split_string: # or word in line.split() to search for full words
			return split_string[1]

def init(filtre=['']):

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
		progress(num, nb_lines, status='Creating directories')

		split_string = line.split("\t")
		path = 'Results'+'/'+split_string[1]+'/'+split_string[2]+'/'+split_string[3]+'/'+split_string[0]
		if split_string[1] == 'Archaea':
			entity_id = find_ids(ids_files[0], split_string[0])
			if not entity_id == None and 'NC' in entity_id:
				(state, gb_record) = filter(index, filtre, entity_id, path)
				if state:
					create_file(split_string, entity_id, ids, gb_record, path)
					count=count+1
		if split_string[1] == 'Bacteria':
			entity_id = find_ids(ids_files[1], split_string[0])
			if not entity_id == None and 'NC' in entity_id:
				(state, gb_record) = filter(index, filtre, entity_id, path)
				if state:
					create_file(split_string, entity_id, ids, gb_record, path)
					count=count+1
		if split_string[1] == 'Eukaryota':
			entity_id = find_ids(ids_files[2], split_string[0])
			if not entity_id == None and 'NC' in entity_id:
				(state, gb_record) = filter(index, filtre, entity_id, path)
				if state:
					create_file(split_string, entity_id, ids, gb_record, path)
					count=count+1
		if split_string[1] == 'Viruses':
			entity_id = find_ids(ids_files[3], split_string[0])
			if not entity_id == None and 'NC' in entity_id:
				(state, gb_record) = filter(index, filtre, entity_id, path)
				if state:
					create_file(split_string, entity_id, ids, gb_record, path)
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



