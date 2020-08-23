#!/usr/bin/env python
"""
Takes a list of variants in TSV format, and outputs VCFs
"""
import argparse
import collections
import csv
import os
import json
import random
import vcfpy

def get_variants_as_records(varfile):
    """
    Arguments:
    varfile -- list of desired variants as CSV

    Returns:
    recordlist -- list of those variants in vcfpy.Record form
    """
    recordlist = []
    reader = csv.DictReader(varfile, delimiter='\t')
    for row in reader:
        ref, alt = row['REF'], row['ALT']

        if len(alt) == len(ref):
            alttype = vcfpy.SNV
        else:
            alttype = vcfpy.INDEL

        altrecord = vcfpy.Substitution(alttype, alt)

        record = vcfpy.Record(row['CHROM'], int(row['POS']), [], row['REF'], [altrecord], None,
                              ['PASS'], {}, None, None)

        recordlist.append(record)
    return recordlist


def get_patientlist(patientsfile):
    """
    Arguments:
    patientsfile -- list of synthetic patients as phenopackets JSON

    Returns:
    id_list -- list of patient ids
    """
    metadata_records = json.load(patientsfile)
    id_list = [item['subject']['id'] for item in metadata_records]
    return id_list


def generate_variant_assignments(npatients, nvariants, min_occ, max_occ):
    """
    Randomly assigns variants to samples

    Arguments:
    npatients - number of patients
    nvariants - number of variants to be assigned
    min_occ, max_occ - minimum and maximum number of times each variant
                       is to be assigned to a patient

    Returns:
    assignments - a dictionary of lists of variants; assinments[3] is
                  the list of variant indexes assigned to the fourth patient
    """
    assignments = collections.defaultdict(list)
    patientlist = list(range(npatients+1))

    for variant_num in range(nvariants):
        n_occ = random.randint(min_occ, max_occ)
        patient_nums = random.choices(patientlist, k=n_occ)
        for patient_num in patient_nums:
            assignments[patient_num].append(variant_num)

    return assignments

def write_vcf(vcffilename, sample_name, records):
    """
    Generate a VCF with the given records and randomly generated genotypes

    Arguments:
    vcffilename - path to generated file
    records - list of vcfpy.Record describing the variants
    """
    lengths = [249250621, 243199373, 198022430, 191154276, 180915260,
               171115067, 159138663, 146364022, 141213431, 135534747,
               135006516, 133851895, 115169878, 107349540, 102531392,
               90354753,  81195210,  78077248,  59128983,  63025520,
               48129895,  51304566]

    samples = vcfpy.SamplesInfos([sample_name])
    header = vcfpy.Header(samples=samples)
    header.add_line(vcfpy.HeaderLine("fileformat", "VCFv4.3"))
    header.add_line(vcfpy.HeaderLine("fileDate", "20200901"))
    for chrom, length in enumerate(lengths):
        header.add_contig_line({"ID": str(chrom), "assembly": "GRCh37", "length": length})
    header.add_format_line({"ID":"GT", "Number":1, "Type":"String", "Description": "Genotype"})

    with open(vcffilename, 'wb') as vcffile:
        writer = vcfpy.Writer.from_stream(vcffile, header, samples, use_bgzf=True)
        for record in records:
            genotype = random.choice(['0/0', '0/1', '1/1'])
            newrecord = vcfpy.Record(record.CHROM,
                                     record.POS,
                                     record.ID,
                                     record.REF,
                                     record.ALT,
                                     record.QUAL,
                                     record.FILTER,
                                     record.INFO,
                                     ["GT"],
                                     calls=[vcfpy.record.Call(sample_name, {"GT": genotype})])
            writer.write_record(newrecord)
        writer.close()


def main(outdir, varfile, patientsfile, min_occur, max_occur):
    """
    Reads in the variantsfile and patients file, generates the
    random assignment of variants to patients, and outputs the
    vcfs and a script to ingest the variantsets.
    """
    var_records = get_variants_as_records(varfile)
    patient_ids = get_patientlist(patientsfile)

    nvars, npatients = len(var_records), len(patient_ids)
    assignments = generate_variant_assignments(npatients, nvars, min_occur, max_occur)
    for patient_num, patient_id in enumerate(patient_ids):
        sample_name = f"{patient_id}_sample_1"
        outfilename = os.path.join(outdir, f"{patient_id}.vcf.gz")
        write_vcf(outfilename, sample_name, [var_records[i] for i in assignments[patient_num]])

    # now generate the json file for candig-ingest
    candig_patientlist = [ { "Patient": {"patientId": patient_id },
                             "Sample": {"patientId": patient_id, "sampleId": f"{patient_id}_sample_1"} }
                           for patient_id in patient_ids ]
    with open(os.path.join(outdir, "ingest.json"), 'w') as jsonfile:
        json.dump({"metadata": candig_patientlist}, jsonfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Variantset generator.')
    parser.add_argument('output_dir', type=str, default='variantsets',
                         help='Directory for VCF files and ingest script')
    parser.add_argument('--variantsfile', '-v',
                        type=argparse.FileType('r'), default="var_by_chrom.tsv",
                        help='path for list of variants in TSV format')
    parser.add_argument('--patientsfile', '-p',
                        type=argparse.FileType('r'), default="cancogen_phenopackets.json",
                        help='path for list of patients in phenopackets JSON')
    parser.add_argument('--min_num_occurances', '-m',
                        type=int, default=1,
                        help="Minimum number of occurances of each variant across sets")
    parser.add_argument('--max_num_occurances', '-M',
                        type=int, default=10,
                        help="Maximum number of occurances of each variant across sets")
    args = parser.parse_args()

    try:
        os.mkdir(args.output_dir)
    except FileExistsError:
        pass

    main(args.output_dir, args.variantsfile, args.patientsfile,
         args.min_num_occurances, args.max_num_occurances)
