# Makefile for TSQR project on Seagull (Intel oneAPI icx + MKL)

CC      = icx
CFLAGS  = -O3 -std=c11 -Wall -Wextra
LDFLAGS = -qmkl -lm

# Sources
TSQR_SRC   = case1tsqr.c
TEST_SRC   = testtsqr.c
SCALE_SRC  = stqrscaling.c

# Targets
TEST_EXE   = tsqr_test1
SCALE_EXE  = scaling_tsqr

.PHONY: all clean run_test run_scale

all: $(TEST_EXE) $(SCALE_EXE)

$(TEST_EXE): $(TSQR_SRC) $(TEST_SRC) tsqr.h
	$(CC) $(CFLAGS) $(TSQR_SRC) $(TEST_SRC) $(LDFLAGS) -o $@

$(SCALE_EXE): $(TSQR_SRC) $(SCALE_SRC) tsqr.h
	$(CC) $(CFLAGS) $(TSQR_SRC) $(SCALE_SRC) $(LDFLAGS) -o $@

run_test: $(TEST_EXE)
	./$(TEST_EXE)

run_scale: $(SCALE_EXE)
	./$(SCALE_EXE)

clean:
	rm -f $(TEST_EXE) $(SCALE_EXE) *.o scaling.csv
