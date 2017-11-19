# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = gcpy
SOURCEDIR     = source
BUILDDIR      = build

# Rule to copy/update the tutorials from their directory under the main gcpy
# file tree
TUTORIALS_DIR = ../tutorials
TUTORIALS_SRC = $(wildcard $(TUTORIALS_DIR)/*)
tutorials:
	@mkdir -p $(SOURCEDIR)/tutorials
	@cp $(TUTORIALS_SRC) $(SOURCEDIR)/tutorials

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile tutorials

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile tutorials
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

html-noplot:
        $(SPHINXBUILD) -D plot_gallery=0 -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
        @echo
        @echo "Build finished. The HTML pages are in $(BUILDDIR)/html."