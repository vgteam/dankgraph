.. _odgi degree:

#########
odgi degree
#########

Describe the graph in terms of node degree.

SYNOPSIS
========

**odgi degree** [**-i, --idx**\ =\ *FILE*] [*OPTION*]…

DESCRIPTION
===========

The odgi degree command describes the graph in terms of node degree.
In summarization mode, it shows the *node.count*, *edge.count*, *avg.degree*,
*min.degree*, and *max.degree*. One can also specify degree ranges streaming these into
a BED file.

OPTIONS
=======

MANDATORY OPTIONS
--------------

| **-i, --idx**\ =\ *FILE*
| Load the succinct variation graph in ODGI format from this *FILE*. The file name usually ends with *.og*. It also accepts GFAv1, but the on-the-fly conversion to the ODGI format requires additional time!

Summary Options
---------------

| **-S, --summarize**
| Summarize the graph properties and dimensions. Print to stdout the
  node.id and the node.degree.

| **-w, --windows-in**\ =\ *LEN:MIN:MAX*
| Print to stdout a BED file of path intervals where the degree is between *MIN* and *MAX*, merging the ranges not separated by more then *LEN* bp.

| **-W, --windows-out**\ =\ *LEN:MIN:MAX*
| Print to stdout a BED file of path intervals where the degree is outside *MIN* and *MAX*, merging the ranges not separated by more then *LEN* bp.

Threading
---------

| **-t, --threads**\ =\ *N*
| Number of threads to use for parallel operations.

Processing Information
----------------------

| **-P, --progress**
| Print information about the operations and the progress to stderr.

Program Information
-------------------

| **-h, --help**
| Print a help message for **odgi degree**.

..
	EXIT STATUS
	===========
	
	| **0**
	| Success.
	
	| **1**
	| Failure (syntax or usage error; parameter error; file processing
	  failure; unexpected error).
	
	BUGS
	====
	
	Refer to the **odgi** issue tracker at
	https://github.com/pangenome/odgi/issues.
