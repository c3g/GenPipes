SHELL = /bin/sh

PROGRAMS = snpEff1 \
           tophat1

snpEff1: VERSION = 3_3
snpEff1: SOFTWARE_ARCHIVE = snpEff_v$(VERSION).zip
snpEff1: SOFTWARE_DIR = snpEff_$(VERSION)
snpEff1: MODULE_DESC = "Adds snpEff1 to your environment"

define snpEff1_download_extract
wget http://sourceforge.net/projects/snpeff/files/snpEff_v$(VERSION)_core.zip
unzip snpEff_v$(VERSION)_core.zip
endef

tophat1: VERSION = 2.0.4
tophat1: SOFTWARE_ARCHIVE = tophat-$(VERSION).Linux_x86_64.tar.gz
tophat1: SOFTWARE_DIR = tophat-$(VERSION).Linux_x86_64
tophat1: MODULE_DESC = "Adds tophat1 to your environment"

###############################################################################
# Everything below this line is generic and should NOT be modified
###############################################################################

SOFTWARE = $@
#SOFTWARE_ROOT = $(MUGQIC_INSTALL_HOME)/software/$(SOFTWARE)
SOFTWARE_ROOT = $(HOME)/software/$(SOFTWARE)
ARCHIVE = $(SOFTWARE_ROOT)/archive
TMP_DOWNLOAD = $(SOFTWARE_ROOT)/tmp
#MODULE_ROOT = $(MUGQIC_INSTALL_HOME)/modulefiles/mugqic/$(SOFTWARE)
MODULE_ROOT = $(HOME)/modulefiles/mugqic
MODULE_DIR = $(MODULE_ROOT)/$(SOFTWARE)

define init-install
	echo $@
	mkdir -p $(ARCHIVE)
	mkdir -p $(TMP_DOWNLOAD)
	cd $(TMP_DOWNLOAD)
endef

modulefile="\#%Module1.0\n\
proc ModulesHelp { } {\n\
	     puts stderr \"\tMUGQIC - $(MODULE_DESC) \"\n\
}\n\
module-whatis \"MUGQIC - $(MODULE_DESC) \"\n\
\n\
set             root               \$$::env(MUGQIC_INSTALL_HOME)/software/$(SOFTWARE)/$(SOFTWARE)-$(VERSION).Linux_x86_64\n\
prepend-path    PATH               \$$root"

versionfile="\#%Module1.0\n\
set ModulesVersion \"$(VERSION)\""

define install-module
	mkdir -p $(MODULE_DIR)
	echo -e $(subst \n ,\n,$(modulefile)) > $(MODULE_DIR)/$(VERSION)
	echo -e $(subst \n ,\n,$(versionfile)) > $(MODULE_DIR)/$(SOFTWARE).version
endef

#.PHONY: all

all: $(PROGRAMS)

$(PROGRAMS): $@_$(SOFTWARE_DIR)
#$(PROGRAMS):
#echo "software_dir" $(SOFTWARE_DIR)

$(addprefix $(MODULE_ROOT)/,$(PROGRAMS)):
	echo "at: " $@
	echo "PROGRAMS: " $(PROGRAMS)
	$(init-install)
	$(install-module)
