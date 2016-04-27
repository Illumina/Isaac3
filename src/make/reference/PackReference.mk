################################################################################
##
## Isaac Genome Alignment Software
## Copyright (c) 2010-2014 Illumina, Inc.
## All rights reserved.
##
## This software is provided under the terms and conditions of the
## GNU GENERAL PUBLIC LICENSE Version 3
##
## You should have received a copy of the GNU GENERAL PUBLIC LICENSE Version 3
## along with this program. If not, see
## <https://github.com/illumina/licenses/>.
##
################################################################################
##
## file PackReference.mk
##
## brief Defines appropriate rules
##
## author Roman Petrovski
##
################################################################################

# first target needs to be defined in the beginning. Ohterwise includes such as
# Log.mk cause unexpected behavior
firsttarget: all

MAKEFILES_DIR:=@iSAAC_HOME@@iSAAC_FULL_DATADIR@/makefiles

# Import the global configuration
include $(MAKEFILES_DIR)/common/Config.mk

include $(MAKEFILES_DIR)/common/Sentinel.mk

# Import the logging functionalities
include $(MAKEFILES_DIR)/common/Log.mk

# Import the debug functionalities
include $(MAKEFILES_DIR)/common/Debug.mk

include $(MAKEFILES_DIR)/reference/Config.mk

ifeq (,$(REFERENCE_GENOME))
$(error "REFERENCE_GENOME is not defined")
endif

ifeq (,$(OUTPUT_FILE))
$(error "OUTPUT_FILE is not defined")
endif

ORIGINAL_GENOME_FASTA:=$(notdir $(shell $(XSLTPROC) $(GET_ANY_FASTA_PATH_XSL) $(REFERENCE_GENOME)))
GENOME_FASTA:=$(TEMP_DIR)/$(ORIGINAL_GENOME_FASTA)

# in some cases the reference can be spread over multiple .fa files. 
# CONSOLIDATED_REFERENCE_XML ensures there is only one .fa file to pack
CONSOLIDATED_REFERENCE_XML:=$(TEMP_DIR)/sorted-reference.xml

ALL_GENOME_ANNOTATION_PATHS:=\
$(shell $(XSLTPROC) --stringparam ANNOTATION_TYPE_PARAM k-uniqueness $(GET_ANNOTATION_FILE_PATHS_XSL) $(REFERENCE_GENOME))\
$(shell $(XSLTPROC) --stringparam ANNOTATION_TYPE_PARAM k-repeatness $(GET_ANNOTATION_FILE_PATHS_XSL) $(REFERENCE_GENOME))

$(CONSOLIDATED_REFERENCE_XML): $(REFERENCE_GENOME) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(REORDER_REFERENCE) --reference-genome $< --output-directory $(TEMP_DIR) --output-xml $(SAFEPIPETARGET) 

$(OUTPUT_FILE): $(ALL_GENOME_ANNOTATION_PATHS) $(CONSOLIDATED_REFERENCE_XML)
	$(CMDPREFIX) $(TAR) -czvO \
		$(foreach ga, $(ALL_GENOME_ANNOTATION_PATHS), -C $(dir $(ga)) $(notdir $(ga))) \
		-C $(CURDIR)/$(TEMP_DIR) \
		$(notdir $(GENOME_FASTA)) \
		$(notdir $(CONSOLIDATED_REFERENCE_XML)) \
		>$(SAFEPIPETARGET) 

all: $(OUTPUT_FILE)

