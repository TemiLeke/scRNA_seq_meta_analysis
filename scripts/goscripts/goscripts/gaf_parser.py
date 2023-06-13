#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@author: Pieter Moris
'''

import os


def importGAF(path, geneSet):
    """
    Imports a GAF file (gene association format) and generates a dictionary mapping the gene Uniprot AC to the GO ID.
    Only imports genes which are present in the provided (background) gene set, if one is provided.

    Information on the GAF 2.1 format can be found at
        http://geneontology.org/page/go-annotation-file-gaf-format-21

    Parameters
    ----------
    path : str
        The path to the file.
    geneSet : set
        A set containing the Uniprot AC's of all the genes under consideration (background).

    Returns
    -------
    dict of str mapping to set
        A dictionary that maps Uniprot ACs (str) to a set GO IDs.

    Possible improvements:
        Check for `is_obsolete` and `replaced_by`, although the replacement term should be in OBO file as an entry.

        Check for inclusion in provided gene set afterwards using:
            gafDict = {key: value for key, value in gafDict.items() if key in geneSet }
    """

    gafPath = os.path.abspath(path)

    # If no gene set was provided
    if not geneSet:
        with open(gafPath, 'r') as gafFile:
            gafDict = {}
            for line in gafFile:
                if line.startswith('UniProtKB'):
                    splitLine = line.split('\t')  # split column-wise
                    uniprotAC = splitLine[1]
                    goTerm = splitLine[4]
                    goQualifier = splitLine[3]
                    if 'NOT' not in goQualifier:  # ignore annotations with "NOT"
                        # Create new key if Uniprot AC does not already appear in dictionary
                        if uniprotAC not in gafDict:
                            gafDict[uniprotAC] = {goTerm}   # make dictionary containing Uniprot AC as key and
                        else:                               # set of GO's as value
                            gafDict[uniprotAC].add(goTerm)

        print('Retrieved', len(gafDict),
              'annotated Uniprot AC\'s from', gafPath + '\n')

        return gafDict

    # if background gene list was provided, limit the import to these genes
    else:
        with open(gafPath, 'r') as gafFile:
            gafDict = {}
            for line in gafFile:
                if line.startswith('UniProtKB'):
                    splitLine = line.split('\t')  # split column-wise
                    uniprotAC = splitLine[1]
                    goTerm = splitLine[4]
                    goQualifier = splitLine[3]
                    if 'NOT' not in goQualifier:  # ignore annotations with "NOT"
                        if uniprotAC in geneSet:  # only keep genes present in background gene set
                            # Create new key if Uniprot AC does not already appear in dictionary
                            if uniprotAC not in gafDict:
                                gafDict[uniprotAC] = {goTerm}   # make dictionary containing Uniprot AC as key and
                            else:                               # set of GO terms as value
                                gafDict[uniprotAC].add(goTerm)

        print('Retrieved', len(gafDict),
              'annotated (background filtered) Uniprot AC\'s from', gafPath + '\n')

        return gafDict


def createSubsetGafDict(subset, gafDict):
    """
    Generates a dictionary that maps the subset's Uniprot ACs to the GO IDs,
    based on the provided gene subset and the gaf dictionary.

    Parameters
    ----------
    subset : set of str
        A subset of Uniprot ACs of interest.
    gafDict : dict of str mapping to set
        A dictionary that maps Uniprot ACs (str) to a set GO IDs.
        Generated by importGAF().

    Returns
    -------
    dict of str mapping to set
        A dictionary that maps the subset's Uniprot ACs to GO IDs.
    """

    gafSubsetDict = {gene: gafDict[gene] for gene in subset if gene in gafDict}

    return gafSubsetDict


def cleanGafTerms(gafDict, filteredGOdict):
    """
    Remove GO terms that do not belong to the chosen namespace, from the gaf dictionary.
    Also removes genes entirely if none of their associated terms belong to the namespace.

    Parameters
    ----------
    gafDict : dict of str mapping to set
        A dictionary that maps gene Uniprot ACs (str) to a set GO term IDs.
        Generated by importGAF().
    filteredGOdict
        A filtered dictionary of GO objects all belonging to the same namespace.
        Generated by obo_tools.importOBO() followed by obo_tools.filterOnNamespace().

    Returns
    -------
    filteredGafDict
        The gaf dictionary after removal of GO terms belonging to different namespaces.
    """

    filteredGafDict = {}

    for gene, goids in gafDict.items():
        # for the associated GO terms of each gene in the gaf dictionary
        # find its intersection with the namespace-filtered GO terms.
        intersection = goids.intersection(filteredGOdict.keys())
        # if the intersection is not empty
        if intersection:
            # add gene keys and intersected reduced GO term sets to a dictionary
            filteredGafDict[gene] = intersection

    return filteredGafDict
