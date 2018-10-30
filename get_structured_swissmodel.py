import sys, requests, json

# similar to get_structured_pdb.py, this script defines a function that gets the
# best structure for a givet uniprot ID. Here we will use swissmodel to get a
# homology model, instead of a solved experimental structure.


def get_structured_swissmodel(uniprotac):
    '''This function takes a uniprot accession number and fetches the 'best'
    homology model from available from the swissmodel repository.
    The model will be downloaded to the ./homology_models/ folder'''

    # for this we will use the swissmodel REST API.
    # for more information on that see: https://swissmodel.expasy.org/docs/repository_help#smr_api

    requestURL = 'https://swissmodel.expasy.org/repository/uniprot/' + uniprotac + '.json'

    r = requests.get(requestURL, headers={"Accept": "application/json"})
    # we check that this is an entry at swismodel. Otherwise this will raise a 404 error.
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    # read the text as a json object
    # this will return a dictionary with 2 keys: 'query' and 'result'
    structure_info = json.loads(r.text)
    # the value of result is again a dictionary with:
    # dict_keys(['sequence_length', 'uniprot_entries', 'structures', 'sequence', 'md5']
    # the value for structures is a list of structures. Each being a dictionary with
    # dict_keys(['similarity', 'gmqe', 'oligo-state', 'crc64', 'coverage', 'alignment', 'md5', 'from',/
    # ' qmean_norm', 'coordinates', 'to', 'identity', 'template', 'provider', 'qmean', 'method'])

    # Now picking the best structure from this list is not trivial.
    # do we go by coverage, qmean or sequence identity?
    # a good place to start is to pick the one with the highest sequence identity

    # start with the first one
    best_structure = structure_info['result']['structures'][0]
    # check the other guys
    for structure in structure_info['result']['structures']:
        # if the current structure has higher seq id than the previous best
        if float(structure['identity']) > float(best_structure['identity']):
            # update the best
            best_structure = structure

    # This structure is associated with an URL, found in the dictionary
    model_url = best_structure['coordinates']
    # donwload the structure. Note that models can easily dissapear from the
    # swissmodel servers. Or be updated
    r = requests.get(model_url)
    # this way the model will be saved as uniprotAC_template_provider.pdb
    with open('homology_models/' + uniprotac + '_' + best_structure['template'] + '_' + best_structure['provider'] + '.pdb', 'w') as f:
        f.write(r.text)


if __name__ == '__main__':
    print(get_structured_swissmodel(sys.argv[1]))
