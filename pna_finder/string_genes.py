# Module that allows PNA Finder to access STRING database

import sys
import urllib


def STRING_int(genes, species):
    """
    Uses the STRING database to test whether a gene network has more interactions than expected. Parts of this function
    are provided by the STRING database at https://string-db.org/cgi/help.pl
    :param genes: list of genes in interaction network
    :param species: NCBI Taxonomy ID for species
    :return:
    """

    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "ppi_enrichment"

    # Construct the request

    request_url = string_api_url + "/" + output_format + "/" + method + "?"
    request_url += "identifiers=" + "%0d".join(genes)
    request_url += "&" + "species=" + species

    # Call STRING

    try:
        response = urllib.request.urlopen(request_url)
    except urllib.error.HTTPError as err:
        error_message = err.read()
        print(error_message)
        sys.exit()

    # Read and parse the results

    result = response.readline()

    while result:
        print(result)
        result = response.readline()


def STRING_net(gene, species, max_nodes=500, required_score=400, exclude_tm=False):
    """
    Provides a gene interaction network for a given gene in a given species. Parts of this function are provided by the
    STRING database at https://string-db.org/cgi/help.pl
    :param gene: central gene of the interaction network
    :param species: NCBI Taxonomy ID for species
    :param max_nodes: The maximum number of nodes (proteins) to be included in the interaction network
    :param required_score: The minimum score (0-1000) for an interaction to be included in the network
    :param exclude_tm: Option to exclude text mining-based interaction results
    :return:
        net_genes: the genes within the interaction network
        nodes: the number of genes within the network
        edges: the number of interactions between the genes in the network
    """

    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"

    # Construct the request

    request_url = string_api_url + "/" + output_format + "/" + method + "?"
    request_url += "identifiers=%s" % gene
    request_url += "&" + "species=" + str(species)
    request_url += "&" + "add_nodes=" + str(max_nodes) + "&" + "required_score=" + str(required_score)

    try:
        response = urllib.request.urlopen(request_url)
    except urllib.error.HTTPError as err:
        error_message = err.read()
        print(error_message)
        return None, None, None

    # Read and parse the results

    line = response.readline().decode('utf-8')
    interaction_list = []
    direct_interaction_list = []
    net_genes = []
    edges = 0

    while line:
        line_list = line.strip().split("\t")
        g1, g2 = line_list[2], line_list[3]
        if exclude_tm:
            prior = 0.041
            score_list = [float(line_list[ii]) for ii in range(6, len(line_list) - 1) if line_list[ii] != '0']
            score_list_nop = [(float(score) - prior)/(1 - prior) for score in score_list]
            nop_prod = 1
            for score_nop in score_list_nop:
                if score_nop != 0:
                    nop_prod = nop_prod * (1 - score_nop)
            s_tot_nop = 1 - nop_prod
            s_tot = s_tot_nop + prior * (1 - s_tot_nop)

            if s_tot < float(required_score)/1000:
                line = response.readline().decode('utf-8')
                continue

            interaction_list += [line_list]
            if gene in [g1, g2]:
                direct_interaction_list += [g1, g2]

        else:
            interaction_list += [line_list]
            if gene in [g1, g2]:
                direct_interaction_list += [g1, g2]

        line = response.readline().decode('utf-8')

    pairs = []
    for interaction in interaction_list:
        g1, g2 = interaction[2], interaction[3]

        if [g1, g2] in pairs or [g2, g1] in pairs:
            continue
        else:
            pairs += [[g1, g2]]

        if g1 in direct_interaction_list and g2 in direct_interaction_list:
            if g1 not in net_genes:
                net_genes.append(str(g1))
            if g2 not in net_genes:
                net_genes.append(str(g2))
            edges += 1
        else:  # For clustering coefficient calculation
            interaction_list.remove(interaction)

    nodes = len(net_genes)

    return net_genes, nodes, edges
