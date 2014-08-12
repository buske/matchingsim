To run the main patient generation code, use patients/randompatients/generate_patient_pairs.py.

The parameters for running this script are as follows:\n
--data_path PATH, -d PATH: Directory from which to grab required data (hgmd, orphanet, hpo)\n
--vcf_path PATH: If you are also generating infected vcfs, use this flag to specify the directory where the original vcfs are found. Note there must be at least 2 files to generate pairs.\n
--out_path PATH, -o PATH: Output directory for infected vcf files and corresponding hpo\n
--generate {PATIENTS, PAIRS}: Specify if you are generating individual patients, or pairs. The default is pairs.\n
-N num: Number of samples to generate (either number of PATIENTS or number of PAIRS)\n
-I {AD, AR}: Which inheritance patterns are allowed for diseases being sampled. At least one is required, but both can be given\n
-D default_freq: Default frequency for phenotypes if frequency info is not found (default is 1.0)\n
--drop_intronic: Drop intronic variants from HGMD\n
--imprecision: Add imprecision to sampled phenotypes (i.e., randomly push the phenotypes up the hpo)\n
--noise: Add phenotypic noise (random phenotypes)\n
-V: When picking which disease to infect a patient with, sample disease weighted by the number of variants, rather than uniformly over diseases which is the default.\n
--logging{DEBUG,INFO,WARNING,ERROR,CRITICAL}: logging level\n
