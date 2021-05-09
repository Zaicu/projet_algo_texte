#ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/
#https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/
#un des deux sites pour le parsing

#site API
#https://biopython.org/docs/1.75/api/Bio.GenBank.html

#ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank cool aussi
#prendre intron et CDS (voir dans genomic.gff.gz)
#Results\Kingdom\Group\SubGroup\Organism\Organism intron.txt
#test https://www.ncbi.nlm.nih.gov/nuccore/NC_009925.1/

import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO


import urllib.request
import meta
import gzip
import shutil
import threading

import requests
from bs4 import BeautifulSoup

def download_file(url):
	try:
		os.mkdir("data")
	except:
		print("déjà crée")
	file_name = "data/" + url.split('/')[-1]
	u = urllib.request.urlopen(url)
	f = open(file_name, 'wb')
	meta = u.info()
	file_size = int(meta.get_all("Content-Length")[0])
	print("Downloading: %s Bytes: %s" % (file_name, file_size))

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
		print(status)

	f.close()
	with gzip.open(file_name, 'rb') as f_in:
		with open(file_name[:-3], 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
	os.remove(file_name)

	#temporaire
	os.remove(file_name[:-3])

def print_gen(path):
	# get all sequence records for the specified genbank file
	recs = [rec for rec in SeqIO.parse(path, "genbank")]

	# print the number of sequence records that were extracted
	print(len(recs))

	# print annotations for each sequence record
	for rec in recs:
		print(rec.annotations)

	# print the CDS sequence feature summary information for each feature in each
	# sequence record
	for rec in recs:
		feats = [feat for feat in rec.features if feat.type == "CDS"]
		for feat in feats:
			print(feat)

def listFD(url):
	page = requests.get(url).text
	#print(page)
	soup = BeautifulSoup(page, 'html.parser')
	return ([(url + '/' + node.get('href'))[:-1] for node in soup.find_all('a') if (node.get('href').endswith("/") and not "//genomes" in (url + '/' + node.get('href')))], [(url + '/' + node.get('href')) for node in soup.find_all('a') if (node.get('href').endswith("gbff.gz"))])

def go_trought_directory(url):
	print(url[37:])
	try:
		os.mkdir(url[37:])
	except:
		print("déjà crée")
	threads = []
	list = listFD(url)
	for dir in list[0]:
		go_trought_directory(dir)
		#t = threading.Thread(target=go_trought_directory, args=(dir,))
		#threads.append(t)
		#t.start()
	#for t in threads:
		#t.join()
	for files in list[1]:
		open(url.split('/')[-1],"w+")
		download_file(files)

for file in listFD("https://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi")[0]:
	print(file)

url = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/ANME-1_cluster_archaeon_AG-394-G06/all_assembly_versions/GCA_009903405.1_ASM990340v1/GCA_009903405.1_ASM990340v1_genomic.gbff.gz"
download_file(url)

#print_gen("data/GCA_009903405.1_ASM990340v1_genomic.gbff")

url = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank"
go_trought_directory(url)
