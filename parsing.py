from Bio import SeqIO
from Bio.SeqIO import parse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import Entrez
import numpy as np
from multiprocessing import Pool
import threading

import re
import urllib.request
import io
import os.path
import sys
import platform
import datetime # C'est à enlever des qu'on a résolu le problème de date

if platform.system() == "Windows": SEP = "\\"
else: SEP = "/"

lock = threading.Lock()

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
	handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=entity_id) # https://www.ncbi.nlm.nih.gov/books/NBK25499/
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

class myThread (threading.Thread):
	def __init__(self, threadID, name, data):
		threading.Thread.__init__(self)
		self.threadID = threadID
		self.name = name
		self.data = data
	def run(self):
		print("Starting " + self.name)
		parsing(self.data)
		print("Exiting " + self.name)

class Coordinate:
	def __init__(self, coord, len_record):
		coord_parse = re.sub('\[([0-9]+):([0-9]+)\]\((\+|\-)\)', "\g<1>\t\g<2>\t\g<3>", coord)

		# Vérifier que le min et max sont bien dans la longueur de la séquence, et que min < max
		if coord_parse != coord:
			coord_parse = coord_parse.split("\t")
			min = int(coord_parse[0])
			max = int(coord_parse[1])
			if coord_parse[2] == "+":
				complement = False
			else:
				complement = True
			if min >= max or max > len_record:
				min        = None
				max        = None
				complement = None

		else:
			min        = None
			max        = None
			complement = None
		# reverse_complement(), regarder le tutorial sur le site internet https://www.tutorialspoint.com/biopython/biopython_advanced_sequence_operations.htm

class Location:
	def __init__(self, loc, len_record):
		loc_parse = re.sub('join{(.*)}', "\g<1>", loc)

		if loc_parse == loc:
			coordinate = Coordinate(loc, len_record)
			join       = False

		else:
			coord_tab  = loc_parse.replace(" ","").split(",")
			coordinate = []
			for coord in coord_tab:
				coordinate.append(Coordinate(coord, len_record))
			join       = True

		#si join == True, alors vérifier que les séquences du join ne s'entrecroisent pas


def write_seq(file_path, gb_record, group):
	#old_stdout = sys.stdout
	#new_stdout = io.StringIO()
	#sys.stdout = new_stdout
	file           = open(file_path, "w")
	function_group = []
	features       = gb_record.features
	len_record     = len(gb_record.seq)

	for f in features:
		if f.type == group:
			loc = str(f.location)
			print(loc)
			location = Location(loc, len_record)
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

def parsing(data):
	reduced_ids   = data[0]
	reduced_paths = data[1]
	reduced_dates = data[2]
	today_ids     = data[3]
	list_group    = data[4]

	i = 0
	print("reduced_ids : ", len(reduced_ids))
	print("reduced_paths : ", len(reduced_paths))
	print("reduced_dates : ", len(reduced_dates))
	lenght = len(reduced_ids)
	#handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=ids)
	fgroup = ""
	for entity_id, path, date in zip(reduced_ids, reduced_paths, reduced_dates):

		i = i + 1
		print(f'{i}/{lenght}')

		if entity_id in today_ids: #faut changer le truc des dates et faire par group aussi
			continue

		organism_name = os.path.basename(path)

		seq_record = get_record(entity_id)

		if not date_compare(date, seq_record.annotations.get("date")):
			lock.acquire(blocking=True)
			today = open("today", "a")
			today.write(entity_id + '\n')
			today.close()
			lock.release()
			continue

		date_file = open(os.path.join(path, "date.dat"), "w") # mettre à jour le date.dat
		date_file.write(seq_record.annotations.get("date"))
		date_file.close()


		function_group = []
		for f in seq_record.features:
			if f.type in list_group:
				function_group.append(str(f.type))
		function_group = list(dict.fromkeys(function_group)) # retire les doublons dans la liste function_group
		print(function_group)

		for group in function_group:
			print(os.path.join(path, group + " " + organism_name + " " + entity_id + ".txt"))
			write_seq(os.path.join(path, group + " " + organism_name + " " + entity_id + ".txt"), seq_record, group) #le seq_record il faut le réduire pour qu'il corresponde au group uniquement (pour le moment je le garde en entier et je traite dans write_seq)

		lock.acquire(blocking=True)
		today = open("today", "a")
		today.write(entity_id + '\n')
		today.close()
		lock.release()

		#fgroup = fgroup + functiongroup + ","

		#print(paths[i] + '\t' +entity_id + functiongroup, file=index)
	return

def parallelize(reduced_ids, reduced_paths, reduced_dates, today_ids, list_group, func):
	cores = 2 #Number of CPU cores on your system
	partitions = cores #Define as many partitions as you want
	# export CLOUDSDK_PYTHON=/usr/bin/python3.7

	reduced_ids_split   = np.array_split(reduced_ids  , partitions)
	reduced_paths_split = np.array_split(reduced_paths, partitions)
	reduced_dates_split = np.array_split(reduced_dates, partitions)
	data_split = [(reduced_ids, reduced_paths, reduced_dates, today_ids, list_group) for (reduced_ids, reduced_paths, reduced_dates) in zip(reduced_ids_split, reduced_paths_split, reduced_dates_split)]

	# Create new threads
	thread1 = myThread(1, "Thread-1", data_split[0])
	thread2 = myThread(2, "Thread-2", data_split[1])
	#thread3 = myThread(3, "Thread-3", data_split[2])

	# Start new Threads
	thread1.start()
	thread2.start()
	#thread3.start()

	# Wait for the threads
	thread1.join()
	thread2.join()
	#thread3.join()

	return

def associate(ids, paths, dates, directory_parsing, list_group):
	ids    = ids.split(',')
	paths  = paths.split(',')
	dates  = dates.split(',')

	reduced_ids   = []
	reduced_paths = []
	reduced_dates = []

	today_ids     = []

	for id, path, date in zip(ids, paths, dates):
		if directory_parsing in path:
			reduced_ids.append(id)
			reduced_paths.append(path)
			reduced_dates.append(date)

	lenght = len(reduced_ids)

	today_date = datetime.date.today().strftime("%d-%b-%Y").upper()

	if os.path.isfile('today'):
		today      = open("today", "r")
		date_ids   = today.readline().rstrip('\n')
		today.close()
	else:
		date_ids = ""	# Si today n'est pas créer, on le crée dans le prochain else

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

	#today = open("today", "a")

	#Parsing
	parallelize(np.array(reduced_ids), np.array(reduced_paths), np.array(reduced_dates), today_ids, list_group, parsing)

	#today.close()
	#tab_group = fgroup.split(',')
	#paths = path.split(',')

	#for i in range(len(tab_group)-1):
		#print(paths[i] + '\t' + tab_group[i], file=index)

# Inutile maintenant, version de backup au cas où
def find_ids(file, entity_name):
	# iterate over lines, and print out line numbers which contain
	# the word of interest.
	for line in file:
		split_string = line.split("\t")
		if entity_name in split_string: # or word in line.split() to search for full words
			return split_string[1]
	#print(split_string[1])

def name_to_id(ids_files):
	dict = {}
	for file in ids_files:
		for line in file:
			split_string = line.split("\t")
			dict[split_string[5]] = split_string[1]
	return dict

def create_tree(overview_lines, ids_files, logs, prgss):

	ids   = ""
	paths = ""
	dates = ""

	name_to_id_dict = name_to_id(ids_files)

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
	prgss.setVisible(True)

	for num, line in enumerate(overview_lines[1:], 1):
		prgss.signal_update.emit(100*num/nb_lines)
		#prgss.setValue(100*num/nb_lines)
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
			#entity_id = str(find_ids(ids_files[kingdom], split_string[0])) # J'imagine que c'est ça qui prend du temps
			if split_string[0] in name_to_id_dict:
				entity_id = name_to_id_dict[split_string[0]]

		if not entity_id == "None":
			if not os.path.isdir(path):
				os.makedirs(path)

			path_date      = os.path.join(path, 'date.dat')
			path_entity_id = os.path.join(path, entity_id)

			if os.path.isfile(path_date):
				date_file = open(path_date, 'r')
				date      = date_file.readline()
				date_file.close()
			else:
				date = ""

			#file  = open(path_entity_id, "a") #crée le fichier NC # Je sais pas si il faudrait pas faire ca plus tard dans l'algo (je pense que si)
			#file.close()

			ids   = ids   + entity_id + ","
			paths = paths + path      + ","
			dates = dates + date      + ","

	prgss.signal_update.emit(100)

	return (ids[:-1], paths[:-1], dates[:-1])

def download(logs):
	dirPath = "GENOME_REPORTS"
	dirIds  = os.path.join(dirPath, "IDS")
	download_file("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt"     ,dirPath, logs)
	download_file("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/IDS/Archaea.ids"  , dirIds, logs)
	download_file("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/IDS/Bacteria.ids" , dirIds, logs)
	download_file("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/IDS/Eukaryota.ids", dirIds, logs)
	download_file("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/IDS/Viruses.ids"  , dirIds, logs)


def init(logs, prgss, filtre=['']):

	download(logs)
	# Ouvrir le fichier en lecture seule
	overview        = open(os.path.join("GENOME_REPORTS", "overview.txt"        ), "r")
	archaea_ids     = open(os.path.join("GENOME_REPORTS", "IDS", "Archaea.ids"  ), "r")
	bacteria_ids    = open(os.path.join("GENOME_REPORTS", "IDS", "Bacteria.ids" ), "r")
	eukaryota_ids   = open(os.path.join("GENOME_REPORTS", "IDS", "Eukaryota.ids"), "r")
	viruses_ids     = open(os.path.join("GENOME_REPORTS", "IDS", "Viruses.ids"  ), "r")

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

	(ids, paths, dates) = create_tree(overview_lines, ids_files, logs, prgss)
	logs.write("Tree done")

	# A retirer après

	associate(ids, paths, dates, os.path.join("Results"), ["CDS", "centromere", "intron", "mobile_element", "ncRNA", "rRNA", "telomere", "tRNA", "3'UTR", "5'UTR"])
	# os.path.join("Results", "Bacteria", "Terrabacteria_group")
	# os.path.join("Results","Viruses","Other","Geminiviridae") crash


	return (ids, paths, dates)

## --------------------------------------------------------------------------- ##

Entrez.email = "thmslpn@gmail.com"

if __name__ == "__main__":
	#init()
	directory_parsing = os.path.join("Results")
	(ids, paths, dates) = init(sys.stdout)
	#path => reduire mes listes en fonction de ce qui est selectionné
	if not ids == "":
		print("ok")
		list_group = ["CDS", "centromere", "intron", "mobile_element", "ncRNA", "rRNA", "telomere", "tRNA", "3'UTR", "5'UTR"]
		associate(ids, paths, dates, directory_parsing, list_group)
	#parse(['tRNA'])
