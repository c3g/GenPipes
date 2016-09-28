##
ROOT="$MUGQIC_INSTALL_HOME/genomes/greengenes_db/"; mkdir -p $ROOT ; cd $ROOT

## The last version is 13/8
wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
tar -zxvf gg_13_8_otus.tar.gz 
mv gg_13_8_otus 138
rm -f gg_13_8_otus.tar.gz


