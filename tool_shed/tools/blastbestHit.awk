awk ' BEGIN { 
	cpU=1
}
NR == 1 {
	old=$1
	nameQ=$1
	nameS=$2
	idP=$3
	lenH=$4
	eval=$11
}
NR >1 {
	if ($1 != old) {
		print nameQ "\t" nameS
		old=$1
		nameQ=$1
		nameS=$2
		idP=$3
		lenH=$4
		eval=$11
	} else {
		if ($11 < eval) {
			ns="new"
		} else if ($11 == eval) {
			if ($3 > idP) {
				ns="new"
			} else if ($3 == idP) {
				if ($4 > lenH) {
					ns="new"
				} else if ($4 == lenH) {
					ns="add"
				} else {
					ns="no"
				}
			} else {
				ns="no"
			}
		} else {
			ns="no"
		}
		if (ns == "new") {
			nameS=$2
			idP=$3
			lenH=$4
			eval=$11
		} else if (ns == "add") {
			nameS=nameS ":" $2
		}
	}
}
END {
	print nameQ "\t" nameS
} ' $1