
OFILES=mpfit.o
LIBFILE=libmpfit.a

all: $(LIBFILE) testmpfit

clean:
	rm -f $(OFILES) testmpfit $(LIBFILE)

mpfit.o: mpfit.c mpfit.h
	$(CC) -O3 -fPIC -c -o $@ $< $(CFLAGS)

$(LIBFILE): $(OFILES)
	$(AR) r $@ $(OFILES)

testmpfit: testmpfit.c libmpfit.a
	$(CC) -O3 -o $@ $(CFLAGS) testmpfit.c -L. -lmpfit -lm
