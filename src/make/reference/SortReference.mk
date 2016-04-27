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
## file SortReference.mk
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

include $(MAKEFILES_DIR)/reference/FindNeighbors.mk

ifeq (,$(GENOME_FILE))
$(error "GENOME_FILE is not defined")
endif

SORTED_REFERENCE_XML:=sorted-reference.xml

$(CONTIGS_XML): $(GENOME_FILE) $(TEMP_DIR)/.sentinel
	$(CMDPREFIX) $(PRINT_CONTIGS) -g $(GENOME_FILE) >$(SAFEPIPETARGET)

ANNOTATION_XML:=$(TEMP_DIR)/annotation.xml

$(ANNOTATION_XML) : $(ANNOTATION_FILE) $(ANNOTATION_REPEAT_FILE)
	$(CMDPREFIX) echo "\
	<SortedReference>\
		<FormatVersion>$(CURRENT_REFERENCE_FORMAT_VERSION)</FormatVersion>\
		<Annotations>\
			<Annotation Type='k-uniqueness' K='$(NEIGHBORHOOD_DISTANCE)'>\
				<Format>$(ANNOTATION_FORMAT)</Format>\
				<File>$(CURDIR)/$(ANNOTATION_FILE)</File>\
			</Annotation>\
			<Annotation Type='k-repeatness' K='$(NEIGHBORHOOD_DISTANCE)'>\
				<Format>$(ANNOTATION_FORMAT)</Format>\
				<File>$(CURDIR)/$(ANNOTATION_REPEAT_FILE)</File>\
			</Annotation>\
		</Annotations>\
	</SortedReference>\
	" >$(SAFEPIPETARGET)

$(SORTED_REFERENCE_XML): $(CONTIGS_XML) $(ANNOTATION_XML)
	$(CMDPREFIX) $(MERGE_REFERENCES) $(foreach part, $^, -i '$(part)') -o $(SAFEPIPETARGET)

all: $(SORTED_REFERENCE_XML)
	$(CMDPREFIX) $(LOG_INFO) "All done!"


