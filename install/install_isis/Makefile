###############################################################################                                                                                
#                                                                             #
# Makefile for system configuration and installation of default software.     #
#                                                                             #
# Author: Fritz-Walter Schwarm <schwarm@smokingscience.com>                   #
#                                                                             #
# Directories:                                                                #
# 	doc:  Documentation and frequently used commands                      #
# 	etc:  System configuration files                                      #
# 	home: User configuration files                                        #
#                                                                             #
# Targets:                                                                    #
# 	etc:      Create symbolic links in /etc for each file or symlink in   #
# 	          $(root)/etc. If the file exists backup it to .bak file.     #
# 	home:     Create symbolic links in ${HOME} for each file or symlink   #
# 	          in $(root)/home. If the file exists backup it to .bak file. #
# 	software: Compile and install additional software.                    #
#                                                                             #
# Variables:                                                                  #
# 	$(root):     Root directory of system repository (this directory)     #
# 	$(etc):      System configuration directory (usually /etc)            # 
# 	$(home):     User home directory                                      # 
# 	$(softdir): Software installation directory                           #
#                                                                             # 
#                                                                             #
###############################################################################


root := $(realpath $(dir $(firstword $(MAKEFILE_LIST))))
ubuntu_version := $(shell cat /etc/debian_version | sed 's/.*=\([\.0-9]\)/\1/g')

etc := /etc
home := ${HOME}
var := /var
dirs := etc home var

files += $(shell find $(dirs) -type f)
files += $(shell find $(dirs) -type l)

# edit installation directory
softdir := /sw
# edit what you want to install.
# at least required for ISIS: headas, slang, isis, isisscripts in that order
software := headas slang slxfig jed slirp isis isisscripts

# install some important libraries
packages_file := packages_for_isis.txt
packages := `cat $(packages_file) | grep ^[a-zA-Z0-9] | tr '\n' ' '`

# look up the current versions of software for every software required above except headas (need to download yourself)
# and ISIS (will automatically download the correct software)
# version info important for script, which automatically fetches these software versions from online archives
# 
# URLs to check:
#
# HEADAS:                        download to $(softdir)/headas yourself from
#                                http://heasarc.gsfc.nasa.gov/FTP/software/lheasoft/release/
# S-Lang / Jed / SLxfig / slirp: Verion numbers found at http://www.jedsoft.org/snapshots/
# heasoft-6.26.1src_no_xspec_modeldata.tar.gz     
headas_version := 6.26.1# edit
headas_packages := libreadline6-dev libncurses5-dev xorg-dev gcc g++ gfortran perl-modules python-dev
# packages += $(headas_packages)
headas_dir := $(softdir)/headas
headas_current := $(headas_dir)/heasoft-$(headas_version)
headas := $(headas_current)/$(shell ls -t1 $(headas_current) | grep `arch` | head -n 1)

# edit
slang_version := pre2.3.3-9
slang_dir := $(softdir)/slang

# edit
slxfig_version := pre0.2.0-121
slxfig_dir := $(softdir)/slxfig

# edit
jed_version := pre0.99.20-137
jed_dir := $(softdir)/jed

# edit
slirp_version := pre2.0.0-34
slirp_dir := $(softdir)/slirp

isis_dir := $(softdir)/isis

isisscripts_dir := $(softdir)/isisscripts

.PHONY:	info all $(files) $(dirs) apt packages software $(software)

info:
	@echo
	@echo -n "ETC: $(etc)\n\n"
	@echo -n "HOME: $(home)\n\n"
	@echo -n "SOFTDIR: $(softdir)\n\n"
@echo -n "FILES: $(files)\n\n"
	@echo -n "DIRECTORIES: $(dirs)\n\n"
	@echo -n "PACKAGES: $(packages)\n\n"
	@echo -n "SOFTWARE: $(software)\n\n"
	@echo -n "HEADAS: $(headas)\n\n"
	@echo -n "USAGE: make [<FILE>|<DIR>|packages|software]\n\n"

all:
	packages headas slang slxfig isis jed

$(files):
	@$(eval folder = $(firstword $(subst /, ,$@)))
	@$(eval file = $(patsubst $(folder)/%,$($(folder))/%,$@))
	@echo Installing $@ to $(file)
	@mkdir -p $(dir $(file))
	@[ -f $(file) ] || [ -d $(file) ] && [ ! -L $(file) ] && mv $(file) $(file).bak || rm -f $(file)
	@ln -s $(root)/$@ $(file)

$(dirs):
	@$(eval dir_files = $(filter $@/%,$(files)))
	@make $(dir_files)

apt:	$(filter etc/apt/%,$(files))
	sudo mkdir -p /etc/apt/sources.list.d
	sudo rm /etc/apt/sources.list.d/owncloud.list
	sudo sh -c "echo 'deb http://download.opensuse.org/repositories/isv:/ownCloud:/community/xUbuntu_$(ubuntu_version)/ /' > /etc/apt/sources.list.d/owncloud.list"
	-wget -nc -O apt/owncloud_Release.key http://download.opensuse.org/repositories/isv:ownCloud:community/xUbuntu_$(ubuntu_version)/Release.key
#	@mv apt/Release.key apt/owncloud_Release.key
	sudo apt-key add - < apt/owncloud_Release.key

packages:	
	apt-get update
	apt-get upgrade
	apt-get install $(packages)

software:
	$(software)


headas:	
	apt-get update
	apt-get upgrade
	apt-get install $(headas_packages)
	@mkdir -p $(headas_dir) 
	wget -nc -P $(headas_dir) http://heasarc.gsfc.nasa.gov/FTP/software/lheasoft/release/heasoft-$(headas_version)src_no_xspec_modeldata.tar.gz
	@-tar xkvzf $(headas_dir)/heasoft-$(headas_version)src_no_xspec_modeldata.tar.gz -C $(headas_dir)
	@( cd $(headas_dir)/heasoft-$(headas_version)/BUILD_DIR ; ./make clean | tee clean.log 2>&1 )
	@( cd $(headas_dir)/heasoft-$(headas_version)/BUILD_DIR ; ./configure | tee config.out 2>&1 )
	@( cd $(headas_dir)/heasoft-$(headas_version)/BUILD_DIR ; make | tee build.log 2>&1 )
	@( cd $(headas_dir)/heasoft-$(headas_version)/BUILD_DIR ; make install | tee install.log 2>&1 )
	@ln -sfn $(headas) $(headas_dir)/current

slang:
	@mkdir -p $(slang_dir)
	wget -nc -P $(slang_dir) http://www.jedsoft.org/snapshots/slang-$(slang_version).tar.gz
	@tar xvzf $(slang_dir)/slang-$(slang_version).tar.gz -C $(slang_dir)
	@( cd $(slang_dir)/slang-$(slang_version) ; ./configure --prefix=$(softdir) | tee config.out 2>&1 )
	@( cd $(slang_dir)/slang-$(slang_version) ; make | tee build.log 2>&1 )
	@( cd $(slang_dir)/slang-$(slang_version) ; make install | tee install.log 2>&1 )
	@ln -sfn $(softdir)/bin/slsh /usr/bin/


slxfig:	slang
	@mkdir -p $(slxfig_dir)
	wget -nc -P $(slxfig_dir) http://www.jedsoft.org/snapshots/slxfig-$(slxfig_version).tar.gz
	@tar xvzf $(slxfig_dir)/slxfig-$(slxfig_version).tar.gz -C $(slxfig_dir)
	@( cd $(slxfig_dir)/slxfig-$(slxfig_version) ; ./configure --prefix=$(softdir) --with-slang=$(softdir) | tee config.out 2>&1 )
	@( cd $(slxfig_dir)/slxfig-$(slxfig_version) ; make | tee build.log 2>&1 )
	@( cd $(slxfig_dir)/slxfig-$(slxfig_version) ; make install | tee install.log 2>&1 )

jed:	slang
	@mkdir -p $(jed_dir)
	wget -nc -P $(jed_dir) http://www.jedsoft.org/snapshots/jed-$(jed_version).tar.gz
	@tar xvzf $(jed_dir)/jed-$(jed_version).tar.gz -C $(jed_dir)
	@( cd $(jed_dir)/jed-$(jed_version) ; ./configure --prefix=$(softdir) --with-slang=$(softdir) | tee config.out 2>&1 )
	@( cd $(jed_dir)/jed-$(jed_version) ; make clean | tee clean.log 2>&1 )
	@( cd $(jed_dir)/jed-$(jed_version) ; make | tee build.log 2>&1 )
	@( cd $(jed_dir)/jed-$(jed_version) ; make install | tee install.log 2>&1 )
	@ln -sfn $(softdir)/bin/jed /usr/bin/
	@ln -sfn $(softdir)/bin/xjed /usr/bin/
	@ln -sfn $(softdir)/bin/jed-script /usr/bin/


slirp:	slang
	@mkdir -p $(slirp_dir)
	wget -nc -P $(slirp_dir) http://www.jedsoft.org/snapshots/slirp-$(slirp_version).tar.gz
	@tar xvzf $(slirp_dir)/slirp-$(slirp_version).tar.gz -C $(slirp_dir)
	@( cd $(slirp_dir)/slirp-$(slirp_version) ; ./configure --prefix=$(softdir) --with-slang=$(softdir) | tee config.out 2>&1 )
	@( cd $(slirp_dir)/slirp-$(slirp_version) ; make | tee build.log 2>&1 )
	@( cd $(slirp_dir)/slirp-$(slirp_version) ; make install | tee install.log 2>&1 )
	@ln -sfn $(softdir)/bin/slirp /usr/bin/

isis:	
	@mkdir -p $(isis_dir)
	git clone https://github.com/houckj/isis $(isis_dir)	
	@( cd $(isis_dir) ; ./configure --prefix=$(softdir) --with-slang=$(softdir) --with-headas=$(headas) | tee config.out 2>&1 )
	@( cd $(isis_dir) ; make | tee build.log 2>&1 )
	@( cd $(isis_dir) ; make install | tee install.log 2>&1 )
	@( cd $(isis_dir) ; make clean | tee clean.log 2>&1 )
	@ln -sfn $(softdir)/bin/isis /usr/bin/
	@ln -sfn $(softdir)/bin/isis-script /usr/bin/
	@touch ${HOME}/.isisrc
	@echo 'Minimum_Stat_Err=1.e-30;'>> ${HOME}/.isisrc
	@echo 'prepend_to_isis_load_path("/ata-school/tutorials/mwl_tutorial");'>> ${HOME}/.isisrc
	@echo 'prepend_to_isis_module_path("/ata-school/tutorials/mwl_tutorial");'>> ${HOME}/.isisrc

isisscripts:
	@mkdir -p $(isisscripts_dir)
	git clone http://www.sternwarte.uni-erlangen.de/git.public/isisscripts $(isisscripts_dir)
	@mkdir -p $(softdir)/share/isis
	@ln -sfn $(isisscripts_dir)/isisscripts.* $(softdir)/share/isis/

                                                                                                                                                                 
