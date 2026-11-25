.DEFAULT_GOAL := all
.DELETE_ON_ERROR:

LATEXMK ?= latexmk
PDFVIEWER ?= pdfterm

TEX_SRC := src/tex
TEX_BUILD := build/tex
BIN_DIR := build/bin
NOTES_DIR := docs/notes
PROOF_DIR := docs/pf
PROOF_PDF := $(PROOF_DIR)/main.pdf

NOTE_NAMES := prime_power_notes prime_power_notes_concise prime_power_roadmap
TEX_PDFS := $(addsuffix .pdf,$(addprefix $(TEX_BUILD)/,$(NOTE_NAMES) main))
NOTE_PDFS := $(addsuffix .pdf,$(addprefix $(NOTES_DIR)/,$(NOTE_NAMES)))
TOOL_NAMES := scan_q_native multiplier_orders explore_2adic_lifts
TOOLS := $(addprefix $(BIN_DIR)/,$(TOOL_NAMES))

.PHONY: all default notes native main display clean $(NOTE_NAMES) $(TOOL_NAMES)
all: notes main native
default: all
notes: $(NOTE_PDFS)
native: $(TOOLS)
main: $(PROOF_PDF)

$(NOTE_NAMES): %: $(NOTES_DIR)/%.pdf
$(TOOL_NAMES): %: $(BIN_DIR)/%

$(NOTE_PDFS): $(NOTES_DIR)/%.pdf: $(TEX_BUILD)/%.pdf | $(NOTES_DIR)
	cp $< $@

$(PROOF_PDF): $(TEX_BUILD)/main.pdf | $(PROOF_DIR)
	cp $< $@

$(TEX_PDFS): $(TEX_BUILD)/%.pdf: $(TEX_SRC)/%.tex $(TEX_SRC)/prime_power_notes.bib .latexmkrc makefile | $(TEX_BUILD)
	TEXINPUTS="$(abspath $(TEX_SRC)):$$TEXINPUTS" \
	BIBINPUTS="$(abspath $(TEX_SRC)):$$BIBINPUTS" \
	$(LATEXMK) -pdf -g -outdir="$(TEX_BUILD)" "$<"

$(TEX_BUILD) $(BIN_DIR) $(NOTES_DIR) $(PROOF_DIR) data:
	mkdir -p $@

display: prime_power_notes
	$(PDFVIEWER) $(NOTES_DIR)/prime_power_notes.pdf

$(BIN_DIR)/scan_q_native: scripts/scan_q_native.cpp scripts/q_prefix_scan.hpp makefile | $(BIN_DIR) data
	$(CXX) $(CPPFLAGS) -O3 -std=c++17 $(CXXFLAGS) -pthread $< $(LDFLAGS) -lflint -lgmpxx -lgmp $(LDLIBS) -o $@

$(BIN_DIR)/multiplier_orders: scripts/exp-gcd/multiplier_orders.c makefile | $(BIN_DIR) data
	$(CC) $(CPPFLAGS) -O3 -std=c17 -Wall -Wextra $(CFLAGS) -pthread $< $(LDFLAGS) $(LDLIBS) -o $@

$(BIN_DIR)/explore_2adic_lifts: scripts/explore_2adic_lifts.cpp makefile | $(BIN_DIR) data
	$(CXX) $(CPPFLAGS) -O3 -std=c++17 $(CXXFLAGS) -fopenmp $< $(LDFLAGS) -lgmpxx -lgmp $(LDLIBS) -o $@

clean:
	$(RM) -r build
	$(RM) $(NOTE_PDFS) $(PROOF_PDF)
