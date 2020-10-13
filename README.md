# Synthetic CanCOGen variantsets for testing

Generates synthetic CanCOGen variantsets for ingest into the CanDIG v1
variant service, matching the patientIds of datasets generated in
phenopackets format by [CanDIG/CanCOGen_synthetic_data](https://github.com/CanDIG/CanCOGen_synthetic_data).

Lists of variants of interest are provided in a tsv file (here: [var_by_chrom.tsv](./var_by_chrom.tsv)),
and those variants are assigned to a random number of patients, with patient and sample ids matching
those provided in the synthetic data above.

## Use

Check out the repo with --recursive, so that the synthetic phenopackets data will also be available:

```
git clone --recursive https://github.com/CanDIG/CanCOGen_synthetic_variant_data.git
```

Then create a virtual environment with the requirements installed:

```
virtualenv env
source env/bin/activate
pip install -r requirements.txt
```

Now usage information can be gotten with the help command

```
./generate_variantsets.py --help
usage: generate_variantsets.py [-h] [--variantsfile VARIANTSFILE] [--patientsfile PATIENTSFILE] [--min_num_occurances MIN_NUM_OCCURANCES] [--max_num_occurances MAX_NUM_OCCURANCES] output_dir

Variantset generator.

positional arguments:
  output_dir            Directory for VCF files and ingest script

optional arguments:
  -h, --help            show this help message and exit
  --variantsfile VARIANTSFILE, -v VARIANTSFILE
                        path for list of variants in TSV format
  --patientsfile PATIENTSFILE, -p PATIENTSFILE
                        path for list of patients in phenopackets JSON
  --min_num_occurances MIN_NUM_OCCURANCES, -m MIN_NUM_OCCURANCES
                        Minimum number of occurances of each variant across sets
  --max_num_occurances MAX_NUM_OCCURANCES, -M MAX_NUM_OCCURANCES
                        Maximum number of occurances of each variant across sets
```

Files in the release were generated with:

```
./generate_variantsets.py vcfs -v var_by_chrom.tsv -p cancogen_phenopackets.json vcfs
```

Setting parameters apporpriately in [`variantset_ingest.sh`](./variantset_ingest.sh) and
running will register the dataset into a CanDIGv1 variant service.
