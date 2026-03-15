
all: implicito explicito

implicito:
	@echo "=> Compilando Solver Implicito..."
	$(MAKE) -C src/solver_implicito

explicito:
	@echo "=> Compilando Solver Explicito..."
	$(MAKE) -C src/solver_explicito

clean:
	@echo "=> Limpando binários..."
	$(MAKE) -C src/solver_implicito clean
	$(MAKE) -C src/solver_explicito clean
