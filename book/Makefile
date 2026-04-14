# ─────────────────────────────────────────────────────────────────────────────
# Build targets for the FEM lecture book
# Prerequisites: conda activate fem-env  (or prepend CONDA_RUN below)
# ─────────────────────────────────────────────────────────────────────────────

QUARTO        := quarto
BOOK_DIR      := .
SLIDES_DIR    := _slides

CONDA_ENV ?=
ifdef CONDA_ENV
  RUN := conda run -n $(CONDA_ENV) --no-capture-output
else
  RUN :=
endif

.PHONY: all html pdf slides slides-execute execute html-execute pdf-execute clean purge

## all          : HTML + PDF book, plus slides
all: html pdf slides

## html         : render the HTML book only
html:
	$(RUN) $(QUARTO) render $(BOOK_DIR) --to html

## pdf          : render the PDF only
pdf:
	$(RUN) $(QUARTO) render $(BOOK_DIR) --to pdf

## slides       : generate slide .qmd files then render them as revealjs
slides:
	$(RUN) python3 gen_slides.py
	cd $(SLIDES_DIR) && $(RUN) $(QUARTO) render

## slides-execute : re-run notebooks, regenerate slides, then render
slides-execute:
	$(RUN) python3 gen_slides.py
	cd $(SLIDES_DIR) && $(RUN) $(QUARTO) render --execute

## execute      : re-run all notebooks, ignoring the freeze cache (all formats)
execute:
	$(RUN) $(QUARTO) render $(BOOK_DIR) --execute

## html-execute : re-run all notebooks and render HTML
html-execute:
	$(RUN) $(QUARTO) render $(BOOK_DIR) --to html --execute

## pdf-execute  : re-run all notebooks and render PDF
pdf-execute:
	$(RUN) $(QUARTO) render $(BOOK_DIR) --to pdf --execute

## clean        : remove rendered output and slides
clean:
	rm -rf $(BOOK_DIR)/_output $(SLIDES_DIR)

## purge        : remove output, slides AND the freeze cache (full rebuild)
purge:
	rm -rf $(BOOK_DIR)/_output $(SLIDES_DIR) $(BOOK_DIR)/_freeze
