# Research probe only; does not change the solver or its rank predictor.
# Build from the repository root:
#   make -f experiments/Makefile build/dixon_rank_probe
include Makefile
build/dixon_rank_probe: experiments/dixon_rank_probe.c src/dixon/dixon_flint.c $(DIXON_SHARED_LIB)
	$(CC) $(ALL_CFLAGS) -UNDEBUG -o $@ $< -L. -ldrsolve $(FLINT_LIBS) $(SYSTEM_LIBS) $(LDFLAGS) $(RPATH_FLAGS)
