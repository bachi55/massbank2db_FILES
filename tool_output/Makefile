####
#
# The MIT License (MIT)
#
# Copyright 2021 Eric Bach <eric.bach@aalto.fi>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is furnished
# to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
####

# Example call: make scoring INPUT_DIR=test_metfrag/ NUMBER_OF_THREADS=1 -j 4 -Otarget

# Path to the MetFrag jar-file
METFRAG_JAR:=~/.local/bin/MetFrag2.4.5-CL.jar
METFRAG_NICE_VALUE:=10

USAGE_STRING:=USAGE: make target INPUT_DIR=DIRECTORY_WITH_INPUT_DATA [NUMBER_OF_THREADS=..., TIMEOUT_IN_SECS=...]

# Input and output directories
ifndef INPUT_DIR
$(error $(USAGE_STRING))
endif
CONF_EXT:=conf
CAND_EXT:=cands.gz
PEAK_EXT:=peaks
RES_DIR:=$(INPUT_DIR)

# Time out for the MetFrag scoring
ifndef TIMEOUT_IN_SECS
# Defaults to 30min per spectrum
TIMEOUT_IN_SECS:=1800
endif

# Number threads used for each MetFrag instance
ifndef NUMBER_OF_THREADS
NUMBER_OF_THREADS:=4
endif

# Add search paths for the recipes
VPATH:=$(INPUT_DIR)

INPUTS=$(subst $(INPUT_DIR),,$(wildcard $(INPUT_DIR)*/*.conf))
# metfrag/AU_001/AU039212.conf -> AU_001/AU039212.conf
OUTPUTS_CANDIDATES=$(INPUTS:.conf=.csv.gz)
# AU_001/AU039212.conf -> AU_001/AU039212.csv.gz

.PHONY: all scoring

default: all

all: scoring

scoring: $(addprefix $(RES_DIR), $(OUTPUTS_CANDIDATES))

$(RES_DIR)%.csv.gz: %.conf %.peaks %.cands.gz
	@echo "========"
	@echo "Scoring:" $<
	gunzip -fk $(word 3, $^)
	timeout $(TIMEOUT_IN_SECS) \
		nice -$(METFRAG_NICE_VALUE) java -jar $(METFRAG_JAR) \
			ParameterFile=$(word 1, $^) \
			PeakListPath=$(word 2, $^) \
			LocalDatabasePath=$(subst .gz,,$(word 3, $^)) \
			ResultsPath=$(@D) \
			NumberThreads=$(NUMBER_OF_THREADS) \
		|| echo "MetFrag scoring failed. Possibly due to timeout. Error code = $${?}."
	rm -f $(subst .gz,, $(word 3, $^))
	@echo "---------"
	@echo "Simplify:" $(basename $@)
	python simplify_metfrag_output.py $(basename $@) \
		&& gzip $(basename $@) \
		|| touch $@.failed
