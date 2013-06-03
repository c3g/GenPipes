#!/bin/sh
## create a initial raw count matrix from a gtf file


#Usage gtf2tmpMatrix.awk gtf($1) output($2)


awk ' BEGIN{
	FS="\t"
}
{ 
	x=split($9,col,";")
	na="-1"
	for (i = 1 ; i <= x ; i++) {
		split(col[i],info,"\"")
		if (info[1] == " gene_id " || info[1] == "gene_id ") {
			ens=info[2]
		} else if (info[1] == " gene_name ") {
                        na=info[2]
                }
	}
	if (na == "-1") {
		na=ens
	}
	print ens "\t" na
} ' $1 | sort -u > $2 
