include ../../example.mk
include ../../common.mk

OBJ = main.o

# Add the TinyXML2 library to the LIBS variable
LIBS += -ltinyxml2

sph_dlb: OPT := $(filter-out -DTEST_RUN,$(OPT))
sph_dlb: $(OBJ)
	$(CC) $(OPT) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

sph_dlb_test: $(OBJ)
	$(CC) $(OPT) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

all: sph_dlb sph_dlb_test

run: sph_dlb
	mpirun -np 4 ./sph_dlb


.PHONY: clean all run run_test

clean:
	rm -f *.o *~ core sph_dlb sph_dlb_test *.vtp *.pvtp *.txt
