#!/bin/bash

fasta=$1

makeblastdb -in $fasta -parse_seqids -blastdb_version 5 -dbtype prot