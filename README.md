# SVfitMEM
Version of SVfit based on full gg->H, H->tautau matrix element, documented in the paper arXiv:1603.05910 .

# Installation instructions
The SVfitPerformanceStudies package has been developed and tested with CMSSW 7_6_x.
It depends on the following other packages:
	TauAnalysis/SVfitTF

In order to install the code, execute:
	cmsrel CMSSW_7_6_3
	cd CMSSW_7_6_3/src
	cmsenv
	git clone https://github.com/veelken/SVfitMEM TauAnalysis/SVfitMEM
	git clone https://github.com/veelken/SVfitTF TauAnalysis/SVfitTF

Do *not* try to compile the code yet.

The SVfitPerformanceStudies package requires the installation of an external library (VAMP) for numeric integration, which is not included in CMSSW by default.
In order to install this library, execute:
	In order to install this library, execute:
	cd $CMSSW_BASE
	mkdir VAMP
	cd VAMP
	wget http://www.hepforge.org/archive/whizard/vamp-2.2.7.tar.gz
	tar -xzvf vamp-2.2.7.tar.gz
	cd vamp-2.2.7
	mkdir $CMSSW_BASE/VAMP/vamp-2.2.7/prefix
	mkdir $CMSSW_BASE/VAMP/vamp-2.2.7/prefix/share/
	mkdir $CMSSW_BASE/VAMP/vamp-2.2.7/prefix/share/doc
	mkdir $CMSSW_BASE/VAMP/vamp-2.2.7/prefix/share/doc/vamp
	./configure --prefix=$CMSSW_BASE/VAMP/vamp-2.2.7/prefix
	make
	make install
	cp $CMSSW_BASE/VAMP/vamp-2.2.7/prefix/lib/* $CMSSW_BASE/lib/$SCRAM_ARCH

You now need to declare the VAMP library as external library in CMSSW.
In order to do this, execute:
	cp /afs/cern.ch/user/v/veelken/public/vamp.xml $CMSSW_BASE/config/toolbox/slc6_amd64_gcc491/tools/selected/
and edit the vamp.xml file. 
Replace the lines containing "/afs/cern.ch/user/v/veelken/" by the directory where you installed the VAMP library (execute 'echo $CMSSW_BASE/VAMP' to see the full path).
Then execute:	
	scram setup vamp

Now compile the CMSSW part of the code:
	cd $CMSSW_BASE/src
	scram b -j 4

In case of compilation problems, check that the VAMP library has been installed correctly.
In order to do this, execute:
	scram tool info vamp
The output that you should get is:
	Name : vamp
	Version : 2.2.6
	++++++++++++++++++++
	SCRAM_PROJECT=no
	LIB=vamp
	LIBDIR=...../vamp/vamp-2.2.6/src/.libs
	INCLUDE=...../vamp/vamp-2.2.6/src
	USE=f77compiler
In case you do not get this output, it means that CMSSW does not know that you have installed the VAMP library.
Have a look on the twiki https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideScram to see hot to fix it.

In case the VAMP library has been correctly installed and the code still does not compile, please contact christian.veelken AT cern.ch .

