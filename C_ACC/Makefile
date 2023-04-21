SHELL=/bin/sh
CLASS=A
VERSION=
SFILE=config/suite.def

all:
	cd BT; $(MAKE) CLASS=$(CLASS) VERSION=$(VERSION)
clean:
	- rm -f core 
	- rm -f *~ */core */*~ */*.o */npbparams.h */*.obj */*.exe
	- rm -f sys/setparams sys/setparams.h
	- rm -f {DC/,}ADC.{logf,view,dat,viewsz,groupby,chunks}.* 