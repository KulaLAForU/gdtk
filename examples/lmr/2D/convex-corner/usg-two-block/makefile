LMR_DIR := lmrsim

help:
	@echo "Main targets available are:"
	@echo -e "   prep"
	@echo -e "   run"
	@echo -e "   post"

prep: prep-gas grid init

prep-gas: gm-ideal-air.inp
	prep-gas gm-ideal-air.inp ideal-air.gas

grid:
	lmr prep-grid

init:
	lmr prep-sim
	lmr prep-mapped-cells

run:
	lmr run

post: vtk

vtk:
	lmr snapshot2vtk --all

test:
	PYTHONDONTWRITEBYTECODE=1 pytest -v

test-gzip:
	LMR_FILE_FORMAT=gziptext PYTHONDONTWRITEBYTECODE=1 pytest -v
clean:
	- rm -rf $(LMR_DIR)
	- rm ideal-air.gas
