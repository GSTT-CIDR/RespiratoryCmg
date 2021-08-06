__credits__ = "Theory of finding correct taxonomy match taken from Natalie Groves' supernatant script."


from collections import defaultdict
import argparse
import time
import pandas as pd


    def argument_parser():
        parser = argparse.ArgumentParser()
        parser.add_argument('--centrifuge',
                            required = False,
                            action='store',
                            help = "path to centrifuge classification results file")
        parser.add_argument('--tree',
                            action='store',
                            default='/Users/alderc/1-projects/12-GSST/1-Project/2-metagenomics/index_hpvc/hpvc.tree',
                            help = "path to taxonomy tree file (Used to build centrifuge index)")  # path to taxonomy tree
        parser.add_argument('--names',
                            action='store',
                            default='/Users/alderc/1-projects/12-GSST/1-Project/2-metagenomics/index_hpvc/hpvc.names',
                            help="path to taxonomy name (Used to build centrifuge index)")
        parser.add_argument("--metag",
                            action="store_true",
                            help="Use to activate function relevant to metagenomics pipeline")
        return parser


class Taxonomy(object):
    """Finds the full taxonomy rank details for a given identifier.

    Args:
        taxid (int): taxonomy identifier from reference data

    Yields:
        taxonomy (dict): e.g. {'species': 'E.coli', 'genus': 'Escherichia', 'family': 'Enterobacteriaeceae' ...}
        rank (str): closest desired rank for identifier (e.g. subspecies -> species if subspecies not in "tax_ranks")

    """

    def __init__(self, taxid):
        self.taxid = int(taxid)
        self.taxonomy = dict()
        self.rank = None

    def get_taxonomy(self, ref_data, tax_ranks):
        # Check taxonomy identifier is valid
        if self.taxid not in ref_data:
            raise KeyError('Taxonomy identifier %s not found in reference data.' % self.taxid)

        # Check reference taxonomy contains a root
        root = None
        for key, value in ref_data.items():
            if value['name'] == 'root':
                root = key
                break
        if not root:
            raise ValueError("Reference taxonomy must contain a root.")

        # Find result at each rank until reach the top of the taxonomy
        rank = ref_data[self.taxid]['rank']
        taxid = self.taxid
        while taxid > root:
            if rank in tax_ranks:
                self.taxonomy[rank] = ref_data[taxid]['name']
                if not self.rank:
                    self.rank = rank
            taxid = ref_data[taxid]['parent']
            if taxid in ref_data:
                rank = ref_data[taxid]['rank']
            else:
                raise KeyError("Incomplete reference data: could not find %s." % taxid)

        # Add in blank results for missing ranks
        for rank in tax_ranks:
            if rank not in self.taxonomy:
                self.taxonomy[rank] = None


class CentrifugeParser(object):
    """Parses the Centrifuge classification output to find number of reads per taxonomy match.

    Args:
        qscore_threshold (int): threshold for including/excluding reads. Default is 300.
        taxonomy_ranks (list): ranks to restrict results to.

    """
    def __init__(self, qscore_threshold=300, taxonomy_ranks=None):
        if taxonomy_ranks is None:
            taxonomy_ranks = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
        self.qthresh = qscore_threshold
        self.ranks = taxonomy_ranks

    @staticmethod
    def read_file(filepath):
        """Converts tab-delimited file into list of lines.

        Args:
            filepath (str): path to file to read

        Returns:
            lines (list): list of lines split by tab

        """
        with open(filepath, 'r') as file_obj:
            lines = [x.strip().split('\t') for x in file_obj.readlines()]
        return lines

    @staticmethod
    def parse_ref_data(tree_data, name_data):
        """Parses data from taxonomy tree and name table files into the ref_data dictionary.

        Args:
            tree_data (list): e.g. [['123', '|', '12', '|', 'species']]
            name_data (list): e.g. [['123', '|', 'E.coli']]

        Returns:
            ref_data (dict): e.g. {123: {'parent': 12, 'rank': 'species', 'name': 'E.coli'}}
            name_to_rank (dict): e.g. {'E.coli': 'species'}

        """

        ref_data = defaultdict(dict)
        name_info = defaultdict(dict)
        for entry in tree_data:
            taxid = int(entry[0])
            parent_taxid = int(entry[2])
            rank = entry[4]
            ref_data[taxid]['parent'] = parent_taxid
            ref_data[taxid]['rank'] = rank
        for entry in name_data:
            taxid = int(entry[0])
            name = entry[2]
            ref_data[taxid]['name'] = name
            name_info[name] = {'rank': ref_data[taxid]['rank'], 'taxid': taxid}

        return ref_data, name_info

    def find_passed_hits_per_read(self, results):
        """Returns hits which have a valid taxonomy ID and are above the score threshold.

        Args:
            results (list): parsed centrifuge classification file

        Returns:
            passed_hits (dict): taxonomy ids per read e.g. {'read1': {123, 234, 345}}
        """

        passed_hits = defaultdict(set)
        for result in results:
            score = int(result[3])
            taxid = int(result[2])
            if score > self.qthresh and taxid != 0:
                read_id = result[0]
                passed_hits[read_id].add(taxid)

        return passed_hits

    def match_low_specificity(self, taxonomy_list, rank_idx):
        """Works up the taxonomy tree until all hits are in agreement.

        Example:
            ['Ecoli', 'Eferg'] -> 'Escherichia'

        Args:
            taxonomy_list (list): list of Taxonomy class objects for each hit in the results
            rank_idx (int): self.ranks index for rank of first hit

        Returns:
            name (str): name of the agreed result e.g. Escherichia
        """

        agreement = False
        name = None
        while not agreement:
            rank_idx += 1
            if rank_idx > len(self.ranks) - 1:
                return None
            rank = self.ranks[rank_idx]
            if all(t.taxonomy[rank] == taxonomy_list[0].taxonomy[rank] for t in taxonomy_list):
                name = taxonomy_list[0].taxonomy[rank]
                agreement = True

        return name

    def match_high_specificity(self, taxonomy_list, rank_idx):
        """Works down the taxonomy tree to find the most specific hit agreed by all.

        Example:
            ['Ecoli', 'Escherichia'] -> 'Ecoli'

        Args:
            taxonomy_list (list): list of Taxonomy class objects for each hit in the results
            rank_idx (int): self.ranks index for rank of first hit

        Returns:
            name (str): name of most specific hit found e.g. E.coli
        """

        current_rank = self.ranks[rank_idx]
        name = taxonomy_list[0].taxonomy[current_rank]
        most_specific_rank = min([self.ranks.index(t.rank) for t in taxonomy_list])
        rank_idx -= 1
        while rank_idx > most_specific_rank - 1:
            new_list = [n for n in taxonomy_list if self.ranks.index(n.rank) <= rank_idx]
            rank = self.ranks[rank_idx]
            if all(t.taxonomy[rank] == new_list[0].taxonomy[rank] for t in new_list):
                name = new_list[0].taxonomy[rank]
                rank_idx -= 1
            else:
                break

        return name

    def find_taxonomic_match(self, hits, ref_data):
        """Compares all hits for one read and returns single agreed hit.

        Args:
            hits (list): taxonomy ids for hits e.g. [123, 234, 345]
            ref_data (dict): reference dictionary produced from parse_ref_data

        Returns:
            name (str): name of hit e.g. 'E.coli'
        """

        # Get taxonomy info for each hit
        taxonomy_list = []
        for taxid in hits:
            tax_obj = Taxonomy(taxid)
            tax_obj.get_taxonomy(ref_data, self.ranks)
            taxonomy_list.append(tax_obj)

        # If no hits, return None
        if len(taxonomy_list) == 0:
            name = None

        # If one hit, get name of hit
        elif len(taxonomy_list) == 1:
            rank_idx = 0
            name = None
            while not name:
                rank = self.ranks[rank_idx]
                name_at_rank = taxonomy_list[0].taxonomy[rank]
                if name_at_rank:
                    name = name_at_rank
                else:
                    rank_idx += 1

        # If all hits have the same rank, find the least specific result agreed by all
        elif all(t.rank == taxonomy_list[0].rank for t in taxonomy_list):
            rank_idx = self.ranks.index(taxonomy_list[0].rank)
            name = self.match_low_specificity(taxonomy_list, rank_idx)

        # If hits have varying ranks
        else:
            rank_idx = max([self.ranks.index(t.rank) for t in taxonomy_list])
            rank = self.ranks[rank_idx]
            # if all agree at maximum rank, work down to find most specific result
            if all(t.taxonomy[rank] == taxonomy_list[0].taxonomy[rank] for t in taxonomy_list):
                name = self.match_high_specificity(taxonomy_list, rank_idx)
            # else work up to find least specific result
            else:
                name = self.match_low_specificity(taxonomy_list, rank_idx)

        return name

    def collate_cfg_results(self, cfg_result, tree_file, name_file):
        """Converts raw Centrifuge result into a dictionary of reads per result.

        Args:
            cfg_result (str): path to file containing raw Centrifuge classification output.
            tree_file (str): path to file containing taxonomy tree used to generate Centrifuge database.
            name_file (str): path to file containing name table used to generate Centrifuge database.

        Returns:
            classified_count (int): total number of classified reads
            results (dict): e.g. key is (result, rank) and value is number of reads e.g. {('Ecoli', 'species'): 5}
        """

        tree_data = self.read_file(tree_file)
        name_data = self.read_file(name_file)
        ref_data, name_info = self.parse_ref_data(tree_data, name_data)

        all_hits = self.read_file(cfg_result)[1:]
        passed_hits = self.find_passed_hits_per_read(all_hits)
        classified_count = 0
        results = defaultdict(int)
        for read, hits in passed_hits.items():
            result = self.find_taxonomic_match(hits, ref_data)
            if result:
                rank = name_info[result]['rank']
                taxid = name_info[result]['taxid']
                classified_count += 1
                results[(result, taxid, rank)] += 1

        return classified_count, results

    def meta_cfg_output(self, cfg_result, tree_file, name_file):
        """Converts raw Centrifuge result into a dictionary of reads per result.

        Args:
            cfg_result (str): path to file containing raw Centrifuge classification output.
            tree_file (str): path to file containing taxonomy tree used to generate Centrifuge database.
            name_file (str): path to file containing name table used to generate Centrifuge database.

        Returns:
            classified_count (int): total number of classified reads
            results (dict): e.g. key is (result, rank) and value is number of reads e.g. {('Ecoli', 'species'): 5}
        """

        tree_data = self.read_file(tree_file)
        name_data = self.read_file(name_file)
        ref_data, name_info = self.parse_ref_data(tree_data, name_data)

        all_hits = self.read_file(cfg_result)[1:]
        passed_hits = self.find_passed_hits_per_read(all_hits)
        classified_count = 0
        results = defaultdict(int)
        reads = defaultdict(list)
        for read, hits in passed_hits.items():
            result = self.find_taxonomic_match(hits, ref_data)
            if result:
                rank = name_info[result]['rank']
                taxid = name_info[result]['taxid']
                classified_count += 1
                results[(result, taxid, rank)] += 1
                reads[read] = [result, taxid, rank]

        return classified_count, results, reads


def print_results(read_count, result_dict):
    """Prints results to stdout."""
    print('Classified reads:\t%s' % read_count)
    for k, v in result_dict.items():
        fraction = float(v) / float(read_count)
        pct = '%.1f' % (fraction * 100)
        print('{result}\t{taxid}\t{rank}\t{reads}\t{pct}'.format(result=k[0], taxid=k[1], rank=k[2], reads=v, pct=pct))


def write_results(read_count, result_dict, out_res):
    """Writes results to file"""
    df_list = list()
    for k, v in result_dict.items():
        fraction = float(v) / float(read_count)
        pct = '%.1f' % (fraction * 100)
        df_list.append([k[0], k[1], k[2], v, pct])
    df = pd.DataFrame(df_list, columns=["Tax_name","Tax_ID", "Rank", "Reads", "Percentage"])
    df.to_csv(out_res, sep="\t", index=False)
    return None
    # for k, v in result_dict.items():
    #     fraction = float(v) / float(read_count)
    #     pct = '%.1f' % (fraction * 100)
    #     print('{result}\t{taxid}\t{rank}\t{reads}\t{pct}'.format(result=k[0], taxid=k[1], rank=k[2], reads=v, pct=pct))

# def metag_df(read_dict):
#     """
#     Creates df of centrifuge output assigning reads with multiple hits with LCA i.e One hit per read.
#     Score is excluded but all reads were above the minimum threshold (default: 300)
#     Parameters
#     ----------
#     read_dict : dict
#         dict with format read_id:(Tax_name, Tax_id, Rank)
#
#     Returns
#     -------
#     df
#         centrifuge style df of read id and taxonomic identification
#     """
#     df = pd.DataFrame.from_dict(read_dict, orient="index", columns = ["Name", "Tax_ID", "Rank"])






if __name__ == '__main__':
    args = argument_parser().parse_args()
    cfg_parser = CentrifugeParser()
    params = snakemake.config["parameters"]["centrifuge"]
    classified_reads, reads_per_result, reads_df = cfg_parser.meta_cfg_output(snakemake.input[0],
                                                                              params["tree"],
                                                                              params["names"])
    write_results(classified_reads, reads_per_result, snakemake.output.table)
