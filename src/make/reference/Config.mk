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
## file Config.mk
##
## brief Common configuration file for all makefiles.
##
## Defines paths, file names, extensions, etc.
##
## author Roman Petrovski
##
################################################################################

MASK_FILE_XML_SUFFIX=.xml
MASK_FILE_SUFFIX=.dat

ANNOTATION_BITS:=8
ANNOTATION_MASK_FILE_SUFFIX:=.$(ANNOTATION_BITS)bpb.gz
ANNOTATION_FORMAT:=$(ANNOTATION_BITS)bpb

NEIGHBORS_BITS:=16
NEIGHBORS_FILE_SUFFIX:=.$(NEIGHBORS_BITS)bpb.gz
NEIGHBORS_FILE_FORMAT:=$(NEIGHBORS_BITS)bpb

CURRENT_REFERENCE_FORMAT_VERSION:=8

GENOME_NEIGHBORS_PREFIX:=neighbors-1or2-
GENOME_NEIGHBORS_SUFFIX:=.1bpb
GENOME_NEIGHBORS_FORMAT:=1bpb

CONTIGS_XML:=$(TEMP_DIR)/contigs.xml
