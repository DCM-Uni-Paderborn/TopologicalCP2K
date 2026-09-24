LATEXMK ?= latexmk
BUILD_DIR ?= .build

.PHONY: all check
all:
	mkdir -p $(BUILD_DIR)/main $(BUILD_DIR)/si
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error -outdir=$(BUILD_DIR)/main main.tex
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error -outdir=$(BUILD_DIR)/si supporting_information.tex
	cp $(BUILD_DIR)/main/main.pdf main.pdf
	cp $(BUILD_DIR)/si/supporting_information.pdf supporting_information.pdf

check: all
	grep -E 'Overfull|There were undefined|multiply defined|LaTeX Error' $(BUILD_DIR)/main/main.log $(BUILD_DIR)/si/supporting_information.log; test $$? -eq 1
