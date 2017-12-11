include env.in

MODULES   := tools
SRC_DIR   := $(addprefix src/, $(MODULES))
SRC       := $(foreach sdir, $(SRC_DIR), $(wildcard $(sdir)/*.cpp))
OBJ       := $(patsubst src/%.cpp,build/%.o, $(SRC))

INCLUDES  := -I./src

BUILD_DIR := $(addprefix build/, $(MODULES) )
BUILD_APP_DIR := $(addprefix build/, apps)

vpath %.cpp $(SRC_DIR)
vpath %.cpp apps

define make-genergal
$1/%.o: %.cpp
	$(MPICC) $(INCLUDES) -c $$< -o $$@
endef

.PHONY: all checkdirs clean

all: checkdirs build/dsvd.mpi

build/dsvd.mpi: build/apps/dsvd.o $(OBJ)
	$(MPICC) $^ -o $@ $(SCALAPACK)

checkdirs: $(BUILD_DIR) $(BUILD_APP_DIR)

$(BUILD_DIR):
	@mkdir -p $@

$(BUILD_APP_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(BUILD_DIR) $(BUILD_APP_DIR)

$(foreach bdir, $(BUILD_DIR), $(eval $(call make-genergal, $(bdir))))
$(foreach bdir, $(BUILD_APP_DIR), $(eval $(call make-genergal, $(bdir))))
