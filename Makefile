include ../../example.mk
include ../../common.mk

OBJ = main.o

sph_dlb: OPT := $(filter-out -DTEST_RUN,$(OPT))
sph_dlb: $(OBJ)
	$(CC) $(OPT) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

sph_dlb_test: $(OBJ)
	$(CC) $(OPT) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

all: sph_dlb sph_dlb_test

run: sph_dlb
	mpirun --oversubscribe -np 2 ./sph_dlb

run_test: sph_dlb_test
	mpirun --oversubscribe -np 2 ./sph_dlb_test

.PHONY: clean all run run_test

clean:
	rm -f *.o *~ core sph_dlb sph_dlb_test *.vtp *.pvtp
