awk -v out=$2 -v logF=$3 ' {
	if ($0 != "") {
		if ($1 == "#") {
			outB=out
			for (i = 2; i <= NF ; i++) {
				gsub("/","",$i)
				if ($i != ":") {
					outB=outB "." $i
				}
			}
			outB=outB ".csv"
			print outB > logF
		} else {
			print $0 > outB
		}
	}
} ' $1