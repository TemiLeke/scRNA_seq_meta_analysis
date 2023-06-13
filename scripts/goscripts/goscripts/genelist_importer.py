#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@author: Pieter Moris
'''

import os


def importGeneList(path):
    """
    Imports the interest/background set (Uniprot AC).

    Parameters
    ----------
    path : str
        The path to the file.

    Returns
    -------
    set of str
        A set of background Uniprot AC's.

    Notes: Gene lists should not contain a header. One gene per line.

    Possible improvements: check for file structure and allow headers, comma separated lists, etc.
    """

    listPath = os.path.abspath(path)

    with open(listPath, 'r') as inGenes:
        geneSet = set([line.rstrip() for line in inGenes])

    print('Retrieved', len(geneSet), 'Uniprot AC\'s from', listPath)

    return geneSet


def isValidSubset(subset, background):
    """
    Checks if the gene subset of interest contains genes not present in the background set.
    If there are additional genes they are removed.

    Parameters
    ----------
    subset : set of str
        A subset of Uniprot ACs of interest.
    background : set of str
        A set of Uniprot ACs to be used as the background.

    Returns
    -------
    set of str
        A cleaned subset of Uniprot ACs of interest.
    """

    if subset.issubset(background):
        return subset
    else:
        missing = [AC for AC in subset if AC not in background]
        print('WARNING! The subset of interest contained genes not present in the background list.')
        print(missing)
        print('Removing these genes from the set of interest...\n')
        return subset.difference(missing)


def reportMissingGenes(geneSet, gafDict, indicator):
    """
    Finds and reports Uniprot AC's in the provided background/interest gene sets which
    are not present in the gene association file (most likely obsolete entries).
    Also returns a new set where these missing genes are removed.

    Parameters
    ----------
    geneSet : set of str
        A set of Uniprot ACs.
        Generated by importSubset() or importBackground().
    gafDict : dict of str mapping to set
        A dictionary mapping gene Uniprot AC's (str) to a set GO IDs.
        Generated by importGAF().
    indicator : str
        A string signifying whether the set is the background or interest set of genes.

    Returns
    -------
    geneSet : set of str
        The set after removal of Uniprot AC's not present in provided gene lists.
    """

    if len(geneSet) != len(gafDict):
        obsoleteGenes = [AC for AC in geneSet if AC not in gafDict]
        print('WARNING! The', indicator, 'gene set contained genes not present in the gene association file:')
        print(obsoleteGenes)
        print('Removing these genes from the', indicator, 'set...\n')
        return geneSet.difference(obsoleteGenes)
    else:
        return geneSet
