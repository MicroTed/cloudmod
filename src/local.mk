NETCDF = /Users/Shared/opt/local/netcdf460ser
OUTPUTINC = -I$(NETCDF)/include
OUTPUTLIB = -L$(NETCDF)/lib
OUTPUTOPT = -DNETCDF -DNCFPLUS
LINKOPTS  =  -Wl,-w  -lnetcdf -lnetcdff -lhdf5_hl -lhdf5 -lm -lz 

FC   = mpif90
#OPTS = -O2 -xHost -ip -assume byterecl -fp-model precise -ftz -no-fma
OPTS = -O2 -ip -assume byterecl -fp-model precise -ftz
#OPTS = -O2 -g -ftrapuv -traceback -fpe0 -ip -assume byterecl -fp-model precise -ftz
#OPTS = -O0 -C -g -ftrapuv -traceback -fpe0 -ip -assume byterecl -fp-model precise -ftz
#OPTS = -O1 -g -ftrapuv -traceback -fpe0 -ip -assume byterecl -fp-model precise -ftz
CPP  = cpp -C -P -traditional -Wno-invalid-pp-token -ffreestanding
DM   = -DMPI
