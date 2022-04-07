#!/bin/bash

# Shell script allowing to install simpiTB's dependencies
#
# authors : David Couvin

# current directory:
CURDIR=`pwd` 

# create log file:
LOGFILE=$CURDIR/installation.log
if [ -e $LOGFILE ]; then rm $LOGFILE ; fi

sudo chmod 755 .
sudo chmod 755 *

# begin install
sudo apt-get update -qq
sudo apt-get -y install locales sudo make unzip zlib1g-dev cpanminus gcc bzip2 libncurses5-dev libncursesw5-dev libssl-dev r-base libxml-libxml-perl libgd-gd2-perl bioperl >> $LOGFILE


# python2 and python3
sudo apt-get -y install python2 >> $LOGFILE
sudo apt-get -y install python3 >> $LOGFILE
sudo apt-get -y install python3-pip >> $LOGFILE

# git, blast, emboss, pandas, roary, fast-lineage-caller, fasttree, curl, wget,
sudo apt-get -y install git >> $LOGFILE
sudo apt-get -y install ncbi-blast+ >> $LOGFILE
sudo apt-get -y install emboss >> $LOGFILE
sudo apt-get -y install roary >> $LOGFILE
sudo apt-get -y install prodigal >> $LOGFILE
#sudo apt-get -y install samtools >> $LOGFILE
sudo apt-get -y install fasttree curl wget >> $LOGFILE
pip3 install pandas >> $LOGFILE
pip3 install statistics >> $LOGFILE
pip3 install fast-lineage-caller >> $LOGFILE
pip3 install grapetree >> $LOGFILE

pip3 install tabulate biopython cgecore gitpython python-dateutil >> $LOGFILE
cpan -f -i Bio::SeqIO

# spades
sudo apt-get -y install spades >> $LOGFILE

# snippy
git clone https://github.com/tseemann/snippy.git 
export PATH=${CURDIR}/snippy/bin/snippy:$PATH

wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -vxjf samtools-1.9.tar.bz2
cd samtools-1.9
make
cd .. 
export PATH="$PATH:/usr/bin/samtools-1.9"
sudo cp /usr/bin/samtools-1.9 /usr/bin/samtools

# clone simpiTB
git clone https://github.com/dcouvin/simpiTB.git >> $LOGFILE
cd simpiTB

# create bin folder:
#echo "create $CURDIR/bin folder..." >> $LOGFILE
#if [ ! -d $CURDIR/bin ];then mkdir $CURDIR/bin; fi

export PATH=$CURDIR:$PATH

#cd bin

# prokka
sudo apt-get install libdatetime-perl libxml-simple-perl libdigest-md5-perl git default-jre bioperl
sudo cpan -f -i Bio::Perl
git clone https://github.com/tseemann/prokka.git

# export $PATH prokka ($CURDIR)
export PATH=${CURDIR}/prokka/bin/:$PATH
export PATH=${CURDIR}/prokka/bin/prokka:$PATH

${CURDIR}/prokka/bin/prokka --setupdb

git clone https://github.com/dcouvin/SpolLineages.git >> $LOGFILE
git clone https://github.com/xiaeryu/SpoTyping.git >> $LOGFILE
git clone https://github.com/phglab/MIRUReader.git >> $LOGFILE
git clone https://git@bitbucket.org/genomicepidemiology/resfinder.git >> $LOGFILE
cd resfinder
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git db_resfinder >> $LOGFILE
git clone https://git@bitbucket.org/genomicepidemiology/pointfinder_db.git db_pointfinder >> $LOGFILE
cd ..
#cd ..

# conda installation (if it does not exist)
if ! command -v conda &> /dev/null
then
 cd $HOME
 curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
 bash ./Miniconda3-latest-Linux-x86_64.sh -b -p ${HOME}/miniconda3
 eval "$(/opt/miniconda3/bin/conda shell.bash hook)"

 export PATH=${HOME}/miniconda3/bin:$PATH # Change to match installation location, if not default.
 export PATH=${HOME}/conda/bin:$PATH

 conda upgrade -c defaults --override-channels conda
 conda config --add channels conda-forge
 conda config --add channels defaults
 conda config --add channels r
 conda config --add channels bioconda

 conda install -c bioconda tb-profiler=2.8.1
 conda install -c bioconda mykrobe
 conda install -c bioconda grapetree
 conda install -c conda-forge tabulate
 conda install -c conda-forge biopython
 conda install -c bioconda cgecore
 conda install -c conda-forge gitpython
 conda install -c conda-forge python-dateutil
fi
