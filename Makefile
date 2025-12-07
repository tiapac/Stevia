EXEC    = Stevia
FC      = gfortran
FFLAGS  = -O3 -ffree-line-length-0
LDFLAGS =
LIBS    =
BUILD_DIR = build

SRCS = src/constants.f90 \
	src/parameters.f90 \
	src/interpolation.f90 \
	src/probe.f90 \
	src/main.f90

OBJS = $(patsubst src/%.f90,$(BUILD_DIR)/%.o,$(SRCS))

$(EXEC): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

$(BUILD_DIR)/%.o: src/%.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(EXEC)