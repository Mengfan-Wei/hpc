# OPTIONS
CC     := icc
CFLAGS := -O3 -xSKYLAKE-AVX512 -qopt-zmm-usage=high
# FILES
EXEC := bmatmult
OBJ  := bmatmult.o
# RULES
$(EXEC) : $(OBJ)
	$(CC) -o $@ $^

%.o: %.c
	$(CC) -c $<

# USEFUL PHONY TARGET
.PHONY : clean
clean:
	rm -rf $(EXEC) $(OBJ) 
