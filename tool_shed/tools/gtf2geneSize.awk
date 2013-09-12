#!/bin/sh
## create a gene size in bp file (Name\tTotalSize) from a gtf file


#Usage gtf2geneSize.awk gtf($1) output($2)

awk ' BEGIN {
	FS="\t"
	na="pre"
	st=1000000000000000
	en=0
}
{
	if ($3 == "exon" || $3 == "start_codon") {
		gsub (" ", "",$9)
		x=split($9,col,";")
		for (i=1 ; i<= x ; i++) {
			split(col[i],info,"\"")
			if ("gene_id" == info[1]) {
				geN=info[2]
			}
		}
		if (geN != na) {
			if (na != "pre") {
				print na "\t" en-st
			}
			na=geN
			st=$4
			en=$5
		} else {
			if ($4 < st) {
				st=$4
			}
			if ($5 > en) {
				en=$5
			}
		}
	}
}
END {
	if (na != "pre") {
		print na "\t" en-st
	}
}
' $1 > $2