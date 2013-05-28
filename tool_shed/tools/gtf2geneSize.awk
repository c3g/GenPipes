#!/bin/sh
## create a gene size in bp file (Name\tTotalSize) from a gtf file


#Usage gtf2geneSize.awk gtf($1) output($2)

awk ' BEGIN {
	na="pre"
	si=0
}
{
	if ($3 == "exon" || $3 == "start_codon") {
		split($10,info,"\"")
		if (na != info[2]) {
			if (na != "pre") {
				print na "\t" si
			}
			na=info[2]
			si=$5-$4
		} else {
			si=si+($5-$4)
		}
	}
}
END {
	if (na != "pre") {
		print na "\t" si
	}
}
' $1 > $2
