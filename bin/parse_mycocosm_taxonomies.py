#!/usr/bin/env python3
import os
import re
import sys
import pandas as pd
from urllib.request import urlopen
from mechanize import Browser
from bs4 import BeautifulSoup
import warnings

ids = open(sys.argv[1]).readlines()

file_out = open(sys.argv[2],'w')
file_out.write('taxon_oid\tJGI_name\tncbi_taxonomy\n')
excluded_file = open('parse_mycocosm_excluded_ids.txt','w')

for id in ids:
	id = id.strip()
	id_myco = id.replace('-cluster','')


	url = "https://mycocosm.jgi.doe.gov/"+id_myco+"/"+id_myco+".home.html"
	try:
		page = urlopen(url)
	except:
		warnings.warn('ID '+id+' NO EXIST ON MYCOCOSM WEBSITE')
		excluded_file.write(id+"\tno available on mycocosm database.\n")
		continue
	for li in page:
		li = str(li)
		if "organismName" in li:

			species_name_with_strain = li.split('</div>')[0].split('&bull;')[1].strip()

			species_name = " ".join(li.split('</div>')[0].split('&bull;')[1].strip().split(' ')[0:3])

			if not "(" in species_name:
				species_name = " ".join(species_name.split()[0:2])
			elif "(" in species_name.split()[1]:
				species_name = " ".join([species_name.split()[0],species_name.split()[2]])
			elif "(" in species_name.split()[2]:
				species_name = " ".join([species_name.split()[0],species_name.split()[1]])
	browser = Browser()
	uri = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi"
	browser.open(uri)
	browser.select_form(name="form")
	browser.form['name'] = species_name
	response = browser.submit()
	soup = BeautifulSoup(response)

	taxos = soup.find_all('a', {'alt': ['kingdom','phylum','class','order','family','genus']})
	taxos = [tax.text for tax in taxos]
	species = soup.find_all('a', {'title' : 'species'})
	if len(species) > 0:
		for specie in species:
			taxos.append(specie.text)
	else:
		taxos.append(species_name)

	if len(taxos) == 0:
		warnings.warn('TAXO NO FOUND FOR '+id+' ORGANISM.')
		excluded_file.write(id+"\ttaxo not found on taxonomy browser.\n")
		continue

	else:
		file_out.write(id.strip()+"\t"+species_name_with_strain+"\t"+";".join(taxos)+'\n')

