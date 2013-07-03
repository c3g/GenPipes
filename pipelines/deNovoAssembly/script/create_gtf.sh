grep "^>" $1 | perl -ne '$_ =~ /^>((comp\d+_c\d+)_seq\d+).*len=(\d+).*/;print $1,"\tprotein_coding\texon\t1\t".$3."\t.\t+\t.\tgene_id \"".$2."\"; transcript_id \"".$1."\";\n";' > $2
