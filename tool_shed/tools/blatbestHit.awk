awk ' BEGIN { 
	cpU=1
}
NR == 6 {
	old=$10
	nameQ=$10
	nameS=$14
	lenM=$1
	lenQ=$11
	lenQs=$12
	lenQe=$13
	lenS=$15
	lenSs=$16
	lenSe=$17
	other="-"
}
NR >7 {
	if ($10 != old) {
		print nameQ "\t" nameS "\t" lenM "\t" lenQ "\t" lenQs "\t" lenQe "\t" lenS "\t" lenSs "\t" lenSe "\t" other
		old=$10
		nameQ=$10
		nameS=$14
		lenM=$1
		lenQ=$11
		lenQs=$12
		lenQe=$13
		lenS=$15
		lenSs=$16
		lenSe=$17
		other="-"
	} else {
		if ($1 > lenM) {
			ns="new"
		} else if ($1 == lenM) {
			ns="add"
		} else {
			ns="no"
		}
		if (ns == "new") {
			old=$10
			nameQ=$10
			nameS=$14
			lenM=$1
			lenQ=$11
			lenQs=$12
			lenQe=$13
			lenS=$15
			lenSs=$16
			lenSe=$17
			other="-"
		} else if (ns == "add") {
			if (other == "-") {
				other=nameS 
			} else {
				other=other ":" nameS
			}
		}
	}
}
END {
	print nameQ "\t" nameS "\t" lenM "\t" lenQ "\t" lenQs "\t" lenQe "\t" lenS "\t" lenSs "\t" lenSe "\t" other
} ' $1