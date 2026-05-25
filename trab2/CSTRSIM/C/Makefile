CFLAGS  = -lm -Wall
#DEBUG = -g
D = util
SOURCES_EULER = cstrsim.c $(D)/params.h $(D)/nrutil.c $(D)/ludcmp.c $(D)/gaussrand.c \
	$(D)/newtrap.c $(D)/controle.c $(D)/lista.c $(D)/arquivo.c $(D)/faults.c
SOURCES_RK = $(D)/runkut.c $(D)/params.h $(D)/nrutil.c $(D)/odeint.c $(D)/ludcmp.c $(D)/gaussrand.c \
	$(D)/newtrap.c $(D)/controle.c $(D)/lista.c $(D)/arquivo.c $(D)/faults.c
SOURCES_ORGN = $(D)/orgn.c $(D)/params.h $(D)/nrutil.c $(D)/ludcmp.c $(D)/gaussrand.c \
	$(D)/newtrap.c $(D)/controle.c $(D)/lista.c $(D)/arquivo.c $(D)/faults.c
TARGET = simulator

#Default
all: euler

#Euler (default) - Solve concentrations and temperature with Euler method
euler:
	@for i in $(SOURCES_EULER);do \
	    echo "Compiling object: $$i..."; \
	    gcc -c $(DEBUG) $$i; \
	done
	@echo "Linking objects: $(TARGET)..."
	@gcc $(DEBUG) *.o $(TARGET).c -o $(TARGET) $(CFLAGS)
	make clean
	
clean:
	@rm -fv *.o *~ util/*.gch

allclean: clean
	@rm -fv $(TARGET)

help:
	@echo "usage: make"
	@echo "       [original | runge_kutta]\n"

#run: build
#	./$(TARGET)

#Original  - Solve concentrations and temperature with original equations
original:
	@for i in $(SOURCES_ORGN);do \
	    echo "Compiling object: $$i..."; \
	    gcc -c $(DEBUG) $$i; \
	done
	@echo "Linking objects: $(TARGET)..."
	@gcc *.o $(DEBUG) $(D)/sim_orgn.c -o $(TARGET) $(CFLAGS)
	make clean

#Runge-Kutta - Solve concentrations and temperature with Rung-Kutta method
runge_kutta:
	@for i in $(SOURCES_RK);do \
	    echo "Compiling object: $$i..."; \
	    gcc -c $(DEBUG) $$i; \
	done
	@echo "Linking objects: $(TARGET)..."
	@gcc *.o $(DEBUG) $(D)/sim_rk.c -o $(TARGET) $(CFLAGS)
	make clean
