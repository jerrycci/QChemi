.PHONY: run fmt lint

PY ?= python3
PYTHONPATH := $(PYTHONPATH):$(shell pwd)/src
export PYTHONPATH

run:
	$(PY) examples/bace1_minimal/run_example.py

fmt:
	@echo "Add your formatter (black/ruff) here."

lint:
	@echo "Add your linter (ruff/flake8) here."
