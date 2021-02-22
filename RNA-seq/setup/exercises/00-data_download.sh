#!/bin/bash
set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

mkdir -p ../data/fastq
cd ../data/fastq
for srr_iter in {1652983..1653014}
do
	last_digit=${srr_iter:${#srr_iter} - 1}
	ascp -QT -l 1000m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR165/00${last_digit}/SRR${srr_iter} .
done
