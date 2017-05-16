objects=brent.o find.o modal.o secular.o getmod.o psv.o global_data.o modrt.o inv2x2.o zhr.o solve.o flat.o getpfile.o
#objects=brent.o test.o modal.o secular.o getmod.o psv.o global_data.o inv4x4.o modrt.o inv2x2.o genrt.o
#CFLAG= -lm -L/lib -L/usr/lib
compile = gfortran -c -O3 -B/usr/lib/x86_64-linux-gnu
link = gfortran -unsharedf95
dsp_jx: $(objects)
	$(link) -o $@ $+
#secular:$(objects)
#	${compile} -O $(CFLAG) -o secular $(objects)

#test.o:test.f90 secular.o getmod.o global_data.mod
#	${compile}  test.f90
find.o:find.f90 modal.o zhr.o flat.o global_data.mod getpfile.o
	${compile} find.f90
getmod.o:getmod.f90 global_data.mod
	${compile}  getmod.f90
flat.o:flat.f90 global_data.mod
	${compile} flat.f90
psv.o:psv.f90 global_data.mod
	${compile}  psv.f90
solve.o:global_data.mod
	${compile}  solve.f90
inv2x2.o:inv2x2.f90
	${compile} inv2x2.f90
modrt.o:modrt.f90 global_data.mod
	${compile}  modrt.f90
secular.o:modrt.f90 psv.o global_data.mod
	${compile}  secular.f90  psv.f90 modrt.f90
zhr.o:modrt.f90 psv.o global_data.mod 
	${compile}  zhr.f90
brent.o:brent.f90 global_data.mod
	${compile} brent.f90
getpfile.o:getpfile.f90 global_data.mod
	${compile} getpfile.f90
modal.o:modal.f90 brent.f90 global_data.mod
	${compile} modal.f90
global_data.mod:global_data.f90
	${compile}  global_data.f90
clean:
	-rm $(objects) global_data.mod 
