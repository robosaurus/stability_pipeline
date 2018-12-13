# stability\_pipeline

This is a pipeline for using Rosetta to predict protein stability changes from single amino acid substitutions (AA subs).
Starting from a uniprot accession number (uniprot ac), the pipeline finds a list of suitable experimental structures for maximum coverage.
(or a homology model if there are no experimental structures).
it will pre-relax the structures and submit rosetta jobs for all possible AA subs for each position in the structure.
it will also get all the single AA sub variants associated with the uniprot ac from Clinvar and from Exac.
all the ddgs and clinvar/exac variants will be written to a summary file, for each structure.

the pipeline was written for Professor Kresten Lindorff-Larsen, by me (lasse.nygaard.biochemist@gmail.com) in november/december 2018.
feel free to write me with any questions about the scripts.

## Usage

You can invoke the pipeline for a single uniprot/_ac with:

python3 /path/to/stability_pipeline/predict_stability.py uniprot_ac full//path/to/output/
for example:
python3 /groups/sbinlab/stability_pipeline/predict\_stability.py O00631 /groups/sbinlab/robotron/stabilities
(make sure the output folder exists)

what happens:
-the gene name for the uniprot\_ac is determined, and used to look for single AA variants in exac and clinvar. The varaints are saved as json files at output/uniprot\_ac/clinvar\_sAA\_variants.json and exac\_sAA\_variants.json.   

-using the SIFTS mappings (Structure integration with function, taxonomy and sequence, not to be confused with the variant prediction tool)
it finds experimental structures for the uniprot accession, the structural coverage is written to a file (output/uniprot_ac/experimental\_coverage\_map.csv)
and from these a list of structures is selected, to give maximum coverage of the uniprot sequence with as few structures as possible.
(if no experimental structures exist, a single homology is fetched from the swissmodel repository instead)

-Each of these structures are cleaned and the relevant chain isolated, with the clean\_pdb rosetta script. And mutiles are written detailing saturation mutagenesis.
an alignment between the struture sequence and the uniprot sequence is performed using MUSCLE. (this allows us to determine the uniprot numbering for the structure sequence)
3 sbatch scripts are written to and submitted from output/uniprot\_ac/rosetta\_runs/uac\_structure\_name/ each dependent on the previous one finishing
rosetta\_relax.sbatch (pre-relaxation of the structure)
rosetta\_cartesian\_saturation_mutagenesis.sbatch (the actual ddg predictions, they are submitted with nice 100 (since there will be a lot of jobs))
parse\_ddgs.sbatch (parses the results, collects the clinvar and exac variants and produces the prediction file)

the prediction file is written to output/uniprot\_ac/predictions/

If you are worried about submitting too much to slurm, or if you are just testing some stuff you can use predict\_stability\_no\_submission.py
It behaves the same way as predict\_stability, except it only writes the sbatch files, and does not submit them.

If you have a list of uniprot\_acs you can put the list at the top of the run\_test\_set.py

Set the paths to Rosetta, and utility scripts and programs in the file rosetta\_paths.py.
It is possible to change the command line flags for rosetta, using the flagfiles present the rosetta\_parameters/ folder.

### Prerequisites

The pipeline is written for python3. and relies on some modules. But nothing too exotic. everything is available through conda.
also feel free to use my python3 interpreter at /groups/sbinlab/robotron/software/miniconda3/bin/python3 if you want.

The pipeline needs some utility scripts that come bundled with rosetta (clean\_pdb.py). set the path in rosetta\_paths.py.

It also uses the MUSCLE to align the structures. The path to muscle should also be set in rosett\_paths.py
(Edgar, R.C. (2004) MUSCLE: multiple sequence alignment with high accuracy and high throughput, Nucleic Acids Res. 32(5):1792-179)

I set all the path, and installed muscle. So it should work now. unless something moves.
