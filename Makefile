# ─────────────────────────────────────────────────────────────────────────────
# Top-level Makefile for the FEM course
# Prerequisites: conda activate fem-env
# ─────────────────────────────────────────────────────────────────────────────

CONDA_ENV   ?=
ifdef CONDA_ENV
  RUN := conda run -n $(CONDA_ENV) --no-capture-output
else
  RUN :=
endif

NB_DIR   := notebooks
BOOK_DIR := book

NOTEBOOKS := $(wildcard $(NB_DIR)/*.ipynb)

.PHONY: all \
        notebooks notebooks-execute notebooks-html notebooks-pdf \
        book book-html book-pdf book-execute \
        slides slides-execute \
        docs \
        clean clean-notebooks clean-book clean-slides clean-docs purge

# ── Default ───────────────────────────────────────────────────────────────────
## all : build notebooks (HTML + PDF), book (HTML + PDF), slides, and docs
all: notebooks book slides docs

# ── Notebooks (Quarto book) ───────────────────────────────────────────────────
# type:book → _output/index.html + per-chapter HTML + one combined PDF.
# Quarto handles execution so symbolic_fem_workbench resolves via the active env.

## notebooks : render notebook book to HTML and PDF (executes if needed)
notebooks: notebooks-html notebooks-pdf

## notebooks-html : render HTML book website → notebooks/_output/index.html
notebooks-html:
	$(RUN) quarto render $(NB_DIR) --to html

## notebooks-pdf : render single combined PDF → notebooks/_output/FEM-Companion-Notebooks.pdf
notebooks-pdf:
	$(RUN) quarto render $(NB_DIR) --to pdf

## notebooks-execute : force re-execution of all notebooks (ignore freeze), both formats
notebooks-execute:
	$(RUN) quarto render $(NB_DIR) --execute

# ── Book (Quarto) ──────────────────────────────────────────────────────────────
## book : render Quarto book to both HTML and PDF
book: book-html book-pdf

## book-html : render HTML website only
book-html:
	$(RUN) quarto render $(BOOK_DIR) --to html

## book-pdf : render PDF only
book-pdf:
	$(RUN) quarto render $(BOOK_DIR) --to pdf

## book-execute : re-run all Quarto notebooks (ignore freeze cache), both formats
book-execute:
	$(RUN) quarto render $(BOOK_DIR) --execute

# ── Slides (revealjs) ─────────────────────────────────────────────────────────
## slides : render each chapter as a revealjs presentation → book/_slides/
slides:
	$(MAKE) -C $(BOOK_DIR) slides RUN="$(RUN)"

## slides-execute : re-run notebooks then render slides
slides-execute:
	$(MAKE) -C $(BOOK_DIR) slides-execute RUN="$(RUN)"

# ── Docs (symbolic_fem_workbench) ────────────────────────────────────────────
## docs : render the package documentation website
docs:
	$(RUN) quarto render docs

## clean-docs : remove docs output
clean-docs:
	rm -rf docs/_output

# ── Clean ─────────────────────────────────────────────────────────────────────
## clean : remove all generated output
clean: clean-notebooks clean-book clean-slides clean-docs

## clean-notebooks : remove notebook output
clean-notebooks:
	rm -rf $(NB_DIR)/_output $(NB_DIR)/_pdf

## clean-book : remove Quarto book output
clean-book:
	rm -rf $(BOOK_DIR)/_output

## clean-slides : remove rendered slides
clean-slides:
	rm -rf $(BOOK_DIR)/_slides

## purge : clean everything including freeze caches
purge: clean
	rm -rf $(NB_DIR)/_freeze $(BOOK_DIR)/_freeze

# ── Help ──────────────────────────────────────────────────────────────────────
help:
	@grep -E '^## ' Makefile | sed 's/## /  /'
