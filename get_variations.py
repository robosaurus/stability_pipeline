import json
import requests, sys

# This file is for testing the uniprot RESTful API
# I will try to fetch the variation data


def get_variations(uniprotid):
    requestURL = "https://www.ebi.ac.uk/proteins/api/variation?offset=0&size=100&accession=" + uniprotid

    r = requests.get(requestURL, headers={"Accept": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    responseBody = r.text

    # the response is a json file, with the variation information
    # read it as a python dict

    variation_dictionary = json.loads(responseBody)
    # There will be some parsing at some point.
    # maybe make some dictionaries for exac and clinvar etc..
    print(json.dumps(variation_dictionary, indent=4))


if __name__ == '__main__':
    get_variations(sys.argv[1])
