#!/usr/bin/env python3
#
# Copyright (C) 2018 INRA
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Vincent Darbot INRAE'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'frogs-support@inrae.fr'
__status__ = 'dev'

import os
import sys
import argparse
import pandas as pd
import scipy.stats as stats
from collections import defaultdict
from statsmodels.sandbox.stats.multicomp import multipletests

def read_input_table(table, ranks, pvalue_deseq):

	rank_to_observ_to_count = {}



	for rank in ranks:
		rank_to_observ_to_count[rank] = {}

		for categorie in ['all','over','under']:
			rank_to_observ_to_count[rank][categorie] = defaultdict(int)
	df = pd.read_csv(table, header=0, delimiter='\t')
	NB_TOTAL = len(df.index)

	df_over = df[df["padj"] < pvalue_deseq]
	df_over = df_over[df_over['log2FoldChange'] > 0]
	NB_OVER = len(df_over.index)
	df_under = df[df["padj"] < pvalue_deseq]
	df_under = df_under[df_under['log2FoldChange'] < 0]
	NB_UNDER = len(df_under.index)

	for rank in ranks:
		for observ in df_over[rank]:
			rank_to_observ_to_count[rank]['over'][observ] += 1

	for rank in ranks:
		for observ in df_under[rank]:
			rank_to_observ_to_count[rank]['under'][observ] += 1

	for rank in ranks:
		for observ in df[rank]:
			rank_to_observ_to_count[rank]['all'][observ] += 1

	return NB_TOTAL, NB_OVER, NB_UNDER, rank_to_observ_to_count


def process_stats(nb_total, nb_over, nb_under, rank_to_observ_counts, pvalue, padj_filter):

	to_adj = []
	sign_results = {}
	cur = 0

	for rank, categorie in rank_to_observ_counts.items():
		categorie['all'] = dict(categorie['all'] )
		categorie['under'] = dict(categorie['under'] )
		categorie['over'] = dict(categorie['over'] )

		
		for sens in ["under","over"]:
			if sens == "under":
				nb_sens = nb_under
			else:
				nb_sens = nb_over
			for observ, count in categorie[sens].items():

				ratio, p = stats.fisher_exact([[nb_total, categorie['all'][observ]], [nb_sens, categorie[sens][observ]]])
				if p < pvalue:

					sign_results[cur] = [rank, sens, observ , nb_total, categorie['all'][observ], nb_sens, categorie[sens][observ], ratio, p]
					to_adj.append(p)
					cur+=1


	p_adjusted = multipletests(to_adj, method='bonferroni', alpha=pvalue)

	for j in range(len(list(p_adjusted[1]))):

		if padj_filter:
			padj = list(p_adjusted[1])[j]
			if padj < pvalue:
				sign_results[j].append(padj)

			else:
				del sign_results[j]
				continue
		else:
			sign_results[j].append(padj)

		if padj < 0.001 :
			sign = "***"
		elif 0.001 < padj < 0.01:
			sign = "**"
		elif 0.01 < padj < 0.05:
			sign = "*"
		elif 0.05 < padj < 0.1:
			sign = "."
		else:
			sign = " "
		
		sign_results[j].append(sign)
	
	return sign_results

def write_output(sign_results, output):

	FH_out = open(output, 'wt')

	FH_out.write('Taxonomic_rank\tDESeq2_Log2fc_sens\tTaxonomic_name\tNb_total_clusters\tNb_observ_taxonomic_name_from_total\tNb_Deseq2_differentially_abundant(OVER_or_UNDER)\tNb_observ_taxonomic_names_from_diff_abundant\tRatio_between_total_and_diff_abundant\tPvalue\tPadj\tSignificance\n')

	for id, signs in sign_results.items():
		FH_out.write("\t".join(map(str,signs))+"\n")

###################################################################################################################
###											  MAIN														   ###
###################################################################################################################
if __name__ == "__main__":
	# Manage parameters
	parser = argparse.ArgumentParser(description="Test which taxonomies are significantly overrepresented in differentially expressed clusters (DESeq2 output)")
	parser.add_argument( '-d', '--debug', default=False, action='store_true', help="Keep temporary files to debug program.")
	parser.add_argument( '-v', '--version', action='version', version=__version__)
	# Inputs
	group_input = parser.add_argument_group('Inputs')
	group_input.add_argument('-t', '--table', required=True, help='Table output of DESeq2 (!! padj cutoff must be set to 1.')
	group_input.add_argument('-r', '--ranks', default='Kingdom,Phylum,Class,Order,Family,Genus', help='Taxonomy rank of interest between : "Kingdom",Phylum", "Class", "Order", "Family", "Genus", "Species". To choose differents ranks, separate them by commas (ie kingdom)')
	group_input.add_argument('--pvalue_deseq', default=0.05, type=float, help='pvalue cutoff from DESeq2 to considere a cluster to be differentially significant. [Default: %(default)s]')
	group_input.add_argument('-p','--pvalue', default=0.05, type=float, help='pvalue to considere a taxonomie observation to be differentially significant. [Default: %(default)s]')
	group_input.add_argument('--padj_filter', default=False, action='store_true', help='if flag added, only observations significant after bonferroni correction will be written [Default: False]')
	# Outputs
	group_output = parser.add_argument_group('Outputs')
	group_output.add_argument('-o', '--output', default='taxo_enrichment.tsv', help='BIOM file with added affiliation annotations from blast/needleall and/or RDPtools. [Default: %(default)s]')
	args = parser.parse_args()


	if "," in args.ranks:
		ranks_input = args.ranks.split(',')
	else:
		ranks_input = [args.ranks]

	ranks = []
	for rank in ranks_input:
		if rank in ["Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species"] and rank not in ranks:
			ranks.append(rank)

		else:
			raise Exception('--ranks paramater values : "Kingdom",Phylum", "Class", "Order", "Family", "Genus", "Species"')


	nb_total, nb_over, nb_under, observ_table = read_input_table(args.table,ranks, args.pvalue_deseq)

	sign_results = process_stats(nb_total, nb_over, nb_under, observ_table, args.pvalue, args.padj_filter)

	write_output(sign_results, args.output)