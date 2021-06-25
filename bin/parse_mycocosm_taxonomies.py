#!/usr/bin/env python3
import os
import re
import sys
import pandas as pd
from urllib.request import urlopen
from mechanize import Browser
from bs4 import BeautifulSoup

ids = open(sys.argv[1]).readlines()

file_out = open(sys.argv[2],'w')
file_out.write('taxon_oid\tJGI_name\tncbi_taxonomy\n')

for id in ids:
	id_myco = id.strip().replace('-cluster','')


	url = "https://mycocosm.jgi.doe.gov/"+id_myco+"/"+id_myco+".home.html"
	try:
		page = urlopen(url)
	except:
		print(id_myco, id)

	for li in page:
		li = str(li)
		if "organismName" in li:

			species_name_with_strain = li.split('</div>')[0].split('&bull;')[1].strip()

			species_name = " ".join(li.split('</div>')[0].split('&bull;')[1].strip().split(' ')[0:3])

			if not "(" in species_name:
				species_name = " ".join(species_name.split()[0:2])
			else:
				species_name = " ".join([species_name.split()[0],species_name.split()[2]])
	browser = Browser()
	uri = "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi"
	browser.open(uri)
	browser.select_form(name="form")
	browser.form['name'] = species_name
	response = browser.submit()
	soup = BeautifulSoup(response)

	taxos = soup.find_all('a', {'alt': ['kingdom','phylum','class','order','family','genus','species']})
	taxos = [tax.text for tax in taxos]
	species = soup.find_all('a', {'title' : 'species'})
	for specie in species:
		taxos.append(specie.text)
	file_out.write(id.strip()+"\t"+species_name_with_strain+"\t"+";".join(taxos)+'\n')


# response = urlopen(urljoin(uri, "/gender/genie.php"))
# forms = ParseResponse(response, backwards_compat=False)
# form = forms[0]

# #print form

# form['text'] = 'cheese'
# form['genre'] = ['fiction']

# print(urlopen(form.click()).read())

