#!/usr/bin/python

"""
Unit-tests for chebi2gene
"""

import unittest
import json

from chebi2gene import *


class Chebi2GeneTestCase(unittest.TestCase):
    """ Unit-test class for chebi2gene.

    Beware, the tests are dependant on the content of the 
    database/sparql endpoint. Changing the data in the sparql endpoint
    might results in failed tests.
    """

    def test_convert_to_uniprot_id(self):
        """ Test the convert_to_uniprot_id function ."""
        data = {'key': ['http://url/to/test:1234']}
        output = convert_to_uniprot_id(data)
        self.assertEqual(output, {'key': ['1234']})

    def test_get_exact_chebi_from_search(self):
        """ Test the get_exact_chebi_from_search function ."""
        expected = {'': {'name': '',
                        'syn': ['']
                        }
                    }
        output = get_exact_chebi_from_search('-beta-carotene')
        self.assertEqual(output, expected)

    def test_get_extended_chebi_from_search(self):
        """ Test the get_extended_chebi_from_search function ."""
        expected = {'': {'name': '',
                        'syn': ['']
                        }
                    }
        output = get_extended_chebi_from_search('-beta-carotene')
        self.assertEqual(output, expected)

    def test_get_genes_of_proteins(self):
        """ Test the get_genes_of_proteins function ."""
        data = {'key': ['1234']}
        expected = {'', ['', '']}
        output = get_genes_of_proteins(data)
        self.assertEqual(output, expected)

    def test_get_organism_of_proteins(self):
        """ Test the get_organism_of_proteins function ."""
        data = {'key': ['1234']}
        expected = {'', ['', '']}
        output = get_organism_of_proteins(data)
        self.assertEqual(output, expected)

    def test_get_pathways_of_proteins(self):
        """ Test the get_pathways_of_proteins function ."""
        data = {'key': ['1234']}
        expected = {'', ['', '']}
        output = get_pathways_of_proteins(data)
        self.assertEqual(output, expected)

    def test_get_protein_of_chebi(self):
        """ Test the get_protein_of_chebi function ."""
        data = ''
        expected = {'', ['', '']}
        output = get_protein_of_chebi(data)
        self.assertEqual(output, expected)


if __name__ == '__main__':
    unittest.main()
