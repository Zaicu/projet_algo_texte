from Bio import SeqIO
from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import Entrez

import re
import urllib.request
import io
import os.path
import sys
import platform
import datetime # C'est à enlever des qu'on a résolu le problème de date

if platform.system() == "Windows": SEP = "\\"
else: SEP = "/"

def download_file(url, dir, logs):
	if not os.path.isdir(dir):
		os.makedirs(dir)
	file_name = os.path.join(dir, url.split('/')[-1])
	u = urllib.request.urlopen(url)
	f = open(file_name, 'wb')
	meta = u.info()
	file_size = int(meta.get_all("Content-Length")[0])
	logs.write("Downloading: %s Bytes: %s" % (file_name, file_size))
	#print("Downloading: %s Bytes: %s" % (file_name, file_size))

	file_size_dl = 0
	block_sz = 8192
	while True:
		buffer = u.read(block_sz)
		if not buffer:
			break
		file_size_dl += len(buffer)
		f.write(buffer)
		status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
		status = status + chr(8)*(len(status)+1)
		sys.stdout.write('%s\n' % (status))
		sys.stdout.flush()

	f.close()

def parse_location(location):
	loc = location.replace(":", " ")
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
	return loc


def write_seq(file_path, gb_record):
	#old_stdout = sys.stdout
	#new_stdout = io.StringIO()
	#sys.stdout = new_stdout
	file = open(file_path, "w")
	function_group = []
	features = gb_record.features

	for f in features:
		if f.type == "CDS" or f.type == "centromere" or f.type == "intron" or f.type == "mobile_element" or f.type == "ncRNA" or f.type == "rRNA" or f.type == "telomere" or f.type == "tRNA" or f.type == "3'UTR" or f.type == "5'UTR":
			loc = str(f.location)
			coord = parse_location(loc)
			coord = coord.split(' ')
			int_coord = []
			for num in coord:
				int_coord.append(int(num))
			if not "join" in loc and not "complement" in loc:
				inf, sup = min(int_coord), max(int_coord)
				file.write("%s\t%d..%d\n" % (f.type, inf, sup))
				file.write(str(gb_record.seq[inf:sup]) + '\n\n')
			if "complement" in loc and not "join" in loc:
				inf, sup = min(int_coord), max(int_coord)
				file.write("%s\tcomplement(%d..%d)\n" % (f.type, inf, sup))
				file.write(str(gb_record.seq[inf:sup]) + '\n\n')
			if "join" in loc:
				join(int_coord, gb_record.seq, f, file)
	#output = new_stdout.getvalue()
	#print(output)
	#sys.stdout = old_stdout
	file.close()

def date_convert(date): #mettre la date dans un format exploitable afin de les comparer numériquement
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

def date_compare(old_date, new_date): #compare les dates et renvoie un booléen correspondant au devoir de mise à jour
	if old_date == "":
		return True
	olddate = date_convert(old_date)
	newdate = date_convert(new_date)
	return olddate < newdate

def get_record(entity_id):
	handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=entity_id)
	record = SeqIO.read(handle, "gb")
	handle.close()
	return record

def filters(index, filtre, entity_id, path):

	path_date      = os.path.join(path, 'date.dat')
	path_entity_id = os.path.join(path, entity_id + '.txt')


	with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=entity_id) as handle:
		gb_record = SeqIO.read(handle, "gb")

	if os.path.isfile(path_date):
		new_date = gb_record.annotations.get("date")
		date_file = open(path_date, "r")
		old_date = date_file.readlines()[0]
		date_file.close()
		if not update(old_date, new_date):
			return (False, gb_record)

	if not os.path.isfile(path_date):

		new_date = gb_record.annotations.get("date")
		date_file = open(path_date, "a")
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

	index.write(path_entity_id + functiongroup + '\n')

	if filtre[0] == '':
		return (True, gb_record)
	#for f in gb_record.features:
	record = str(gb_record.features)
	for group in filtre:
		if not group in record:
			return (False, gb_record)
	return (True, gb_record)

def join(coord, sequence, f, file=''):
	#old_stdout = sys.stdout
	#new_stdout = io.StringIO()
	#sys.stdout = new_stdout
	l = len(coord)
	loc = "(" + parse_location(str(f.location)).replace(" ", ",") + ")"
	string = ""
	for i in range(round(l/2)):
		inf = coord[2*i]
		sup = coord[2*i+1]
		file.write("%s\t%s\t%d..%d\n" % (f.type, loc, inf, sup))
		file.write(str(sequence[inf:sup]) + '\n\n')
	#output = new_stdout.getvalue()
	#print(output)
	#sys.stdout = old_stdout

def update(path):
	path           = path.replace(':','_').replace(' ','_')

	if not os.path.isdir(path): # Je crois qu'on est jamais dans ce cas là
		os.makedirs(path)

	# Vérifier si le dossier existe ou non : inutile si les lignes sont distinctes
	for files in os.listdir(path):
		if not files == 'date.dat':
			entity_id = files.replace(".txt", "")

	path_date      = os.path.join(path, 'date.dat')
	path_entity_id = os.path.join(path, entity_id + '.txt')

	record = get_record(entity_id)
	print(path_entity_id)

	new_date = record.annotations.get("date")

	if not os.path.isfile(path_entity_id):
		date_file = open(path_date, "a")
		date_file.write(new_date)
		date_file.close()
		file = open(path_entity_id, "a")

		ids = ids + ',' + entity_id
		write_seq(file, record)
		file.close()

	if not os.path.isfile(path_date):
		date_file = open(path_date, "a")
		date_file.write(new_date)
		date_file.close()
		file = open(path_entity_id, "w")
		write_seq(file, record)
		file.close()

	date_file = open(path_date, "r")
	old_date = date_file.readlines()[0]
	date_file.close()

	if date_compare(old_date, new_date):
		date_file = open(path_date, "w")
		date_file.write(new_date)
		date_file.close()
		file = open(path_entity_id, "w")
											# Comprend pas là
		ids = ids + ',' + entity_id
		write_seq(file, record)
		file.close()

def parse(filtre):
	if platform.system() == "Windows": index = open('index_windows.txt', "r")
	else: index = open('index.txt', "r")
	index_lines = index.readlines()
	index.close()

	print(index_lines)
	for line in index_lines:
		BOOL = 1
		for group in filtre:
			if not group in line:
				BOOL = 0
		if BOOL == 1:
			print(line.split("\t")[0])
			update(line.split("\t")[0])












def progress(count, total, status=''):

	bar_len = 45
	filled_len = int(round(bar_len * count / float(total)))

	percents = round(100.0 * count / float(total), 1)
	bar = '=' * (filled_len) + '>' + '.' * (bar_len - filled_len)

	sys.stdout.write('[%s] %s%s | %s\r' % (bar, percents, '%', status))
	sys.stdout.flush()

def associate(ids, paths, dates):
	ids    = ids.split(',')
	paths  = paths.split(',')
	dates  = dates.split(',')
	lenght = len(ids)

	today      = open("today", "r")
	date_ids   = today.readline().rstrip('\n')
	today.close()
	today_date = datetime.date.today().strftime("%d-%b-%Y").upper()

	today_ids = []

	print(date_ids)
	print(today_date)
	if date_ids == today_date:
		print('a')
		today = open("today", "r")
		today_ids = [line.rstrip('\n') for line in today.readlines()[1:]]
		today.close()
	else:
		print('b')
		today = open("today", "w")
		today.write(today_date + '\n')
		today.close()

	today = open("today", "a")

	i = 0
	#handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=ids)
	fgroup = ""
	for entity_id, path, date in zip(ids, paths, dates):

		i = i + 1
		print(f'{i}/{lenght}')

		if entity_id in today_ids:
			continue

		seq_record = get_record(entity_id)

		if not date_compare(date, seq_record.annotations.get("date")): #Normalement c'est bien le cas dans lequel y a pas besoin de mettre à jour
			today.write(entity_id + '\n')
			continue

		date_file = open(os.path.join(path, "date.dat"), "w") # mettre à jour le date.dat
		date_file.write(seq_record.annotations.get("date"))
		date_file.close()


		function_group = []
		for f in seq_record.features:
			if f.type == "CDS" or f.type == "centromere" or f.type == "intron" or f.type == "mobile_element" or f.type == "ncRNA" or f.type == "rRNA" or f.type == "telomere" or f.type == "tRNA" or f.type == "3'UTR" or f.type == "5'UTR":
				function_group.append(str(f.type))
		function_group = list(dict.fromkeys(function_group)) # retire les doublons dans la liste function_group
		print(function_group)

		print(os.path.join(path, entity_id))
		write_seq(os.path.join(path, entity_id), seq_record)
		today.write(entity_id + '\n')

		#fgroup = fgroup + functiongroup + ","

		#print(paths[i] + '\t' +entity_id + functiongroup, file=index)

	today.close()
	#tab_group = fgroup.split(',')
	#paths = path.split(',')

	#for i in range(len(tab_group)-1):
		#print(paths[i] + '\t' + tab_group[i], file=index)

def find_ids(file, entity_name):
	# iterate over lines, and print out line numbers which contain
	# the word of interest.
	for line in file:
		split_string = line.split("\t")
		if entity_name in split_string: # or word in line.split() to search for full words
			return split_string[1]

def create_tree(overview_lines, ids_files):

	ids   = ""
	paths = ""
	dates = ""

	all_entity_id = []
	r = re.compile(".*txt")
	if os.path.isfile('index.txt'):
		index = open('index.txt', 'r')
		for line in index.readlines():
			#append le NC dans all_entity_id
			entity_id = line.split('\t')[1].split('\n')[0]
			#print(os.listdir(path))
			all_entity_id.append(entity_id)
		index.close()

	nb_lines = len(overview_lines)
	for num, line in enumerate(overview_lines[1:], 1):
		progress(num, nb_lines, status='Creating directories')
		entity_id = "None" # init à chaque boucle

		split_string = line.split("\t")
		new_name = split_string[0].replace("[", "")
		new_name = new_name.replace("]", "")
		new_name = new_name.replace(":", "_")
		new_name = new_name.replace("'", "")
		new_name = new_name.replace("(", "")
		new_name = new_name.replace(")", "")
		new_name = new_name.replace(".", "")

		"""if root_dir == "GENOME_REPORTS":
			path = ('Results'+SEP+split_string[1]+SEP+split_string[2]+SEP+split_string[3]+SEP+split_string[0]).replace(' ','_').replace(':','_')
		else :
			path = root_dir.replace(SEP+"GENOME_REPORTS",'').replace("GENOME_REPORTS",'')+( SEP+'Results'+SEP+split_string[1]+SEP+split_string[2]+SEP+split_string[3]+SEP+split_string[0]).replace(' ','_').replace(':','_')"""

		# Si la ligne n'est pas conforme on ne la traite pas
		if len(split_string) < 4:
			continue

		path = (os.path.join('Results', split_string[1], split_string[2], split_string[3], new_name)).replace(':','_').replace(' ','_')

		if split_string[1] == 'Archaea':
			kingdom = 0
		if split_string[1] == 'Bacteria':
			kingdom = 1
		if split_string[1] == 'Eukaryota':
			kingdom = 2
		if split_string[1] == 'Viruses':
			kingdom = 3

		if kingdom == 0 or kingdom == 1 or kingdom == 2 or kingdom == 3:
			entity_id = str(find_ids(ids_files[kingdom], split_string[0]))

		if not entity_id == "None":
			if not os.path.isdir(path):
				os.makedirs(path)

			path_date      = os.path.join(path, 'date.dat')
			path_entity_id = os.path.join(path, entity_id + '.txt')

			if os.path.isfile(path_date):
				date_file = open(path_date, 'r')
				date      = date_file.readline()
				date_file.close()
			else:
				date = ""

			file  = open(path_entity_id, "a") #crée le fichier NC # Je sais pas si il faudrait pas faire ca plus tard dans l'algo (je pense que si)
			file.close()

			ids   = ids   + entity_id + ","
			paths = paths + path      + ","
			dates = dates + date      + ","


	return (ids[:-1], paths[:-1], dates[:-1])

def download(logs):
	dirPath = "GENOME_REPORTS"
	dirIds  = os.path.join(dirPath, "IDS")
	download_file("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt"     ,dirPath, logs)
	download_file("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/IDS/Archaea.ids"  , dirIds, logs)
	download_file("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/IDS/Bacteria.ids" , dirIds, logs)
	download_file("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/IDS/Eukaryota.ids", dirIds, logs)
	download_file("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/IDS/Viruses.ids"  , dirIds, logs)


def init(logs, filtre=['']):

	download(logs)
	# Ouvrir le fichier en lecture seule
	overview        = open('GENOME_REPORTS/overview.txt'     , "r")
	archaea_ids     = open('GENOME_REPORTS/IDS/Archaea.ids'  , "r")
	bacteria_ids    = open('GENOME_REPORTS/IDS/Bacteria.ids' , "r")
	eukaryota_ids   = open('GENOME_REPORTS/IDS/Eukaryota.ids', "r")
	viruses_ids     = open('GENOME_REPORTS/IDS/Viruses.ids'  , "r")

	# utiliser readlines pour lire toutes les lignes du fichier
	# La variable "lignes" est une liste contenant toutes les lignes du fichier
	overview_lines  = overview.readlines()
	archaea_lines   = archaea_ids.readlines()
	bacteria_lines  = bacteria_ids.readlines()
	eukaryota_lines = eukaryota_ids.readlines()
	viruses_lines   = viruses_ids.readlines()
	ids_files       = [archaea_lines, bacteria_lines, eukaryota_lines, viruses_lines]

	# fermez le fichier après avoir lu les lignes
	overview.close()
	archaea_ids.close()
	bacteria_ids.close()
	eukaryota_ids.close()
	viruses_ids.close()

	nb_lines = len(overview_lines)
	ids      = ''
	count    = 0
	kingdom  = -1

	(ids, paths, dates) = create_tree(overview_lines, ids_files)
	# Itérer sur les lignes sauf la première

	if not ids == "":
		print("ok")
		associate(ids, paths, dates)


## --------------------------------------------------------------------------- ##

Entrez.email = "thmslpn@gmail.com"
mail = ['test@gmail.com','test1@gmail.com','test2@gmail.com','test3@gmail.com','test4@gmail.com','test5@gmail.com','test6@gmail.com','test7@gmail.com','test8@gmail.com','test9@gmail.com',]
#init()
#init()
#parse(['tRNA'])
