from predict_stability import predict_stability_for_ac

output_folder = '/groups/sbinlab/software/stability_pipeline/output/'

for uniprot_accession in ['O00631']:
        predict_stability_for_ac(uniprot_accession, output_folder)
