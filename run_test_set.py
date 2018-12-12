from predict_stability import predict_stability_for_ac

output_folder = '/groups/sbinlab/software/stability_pipeline/output/'

for uniprot_accession in ['P43246', 'P40692', 'P35557', 'P51587', 'P60484', 'Q8NEA6']:
        predict_stability_for_ac(uniprot_accession, output_folder)
