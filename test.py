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
        expected = {'35309': {'name': ['(5S,6R)-beta-carotene 5,6-epoxide'],
                        'syn': ['C40H56O', '(5S,6R)-beta-carotene 5,6-epoxide',
                               '(5S,6R)-5,6-epoxy-5,6-dihydro-beta,beta-carotene']
                        }
                    }
        output = get_exact_chebi_from_search('-beta-carotene')
        self.assertEqual(output, expected)

    def test_get_extended_chebi_from_search(self):
        """ Test the get_extended_chebi_from_search function ."""
        expected = {'17579': {'name': ['beta-carotene'],
                         'syn': ['all-trans-beta-carotene']
                        }
                    }
        output = get_extended_chebi_from_search('trans-beta-carotene')
        self.assertEqual(output, expected)

    def test_get_genes_of_proteins(self):
        """ Test the get_genes_of_proteins function ."""
        data = {'key': ['Q38933']}
        expected = {'Q38933': [
                      {'stop': '31104596', 'start': '31103094',
                       'sca': u'SL2.31ch04', 'name': u'Solyc04g040190.1.1',
                       'desc': u'Beta-lycopene cyclase (AHRD V1 ***- A6YS01_SOLLC)%3B'\
                       ' contains Interpro domain(s)  IPR010108  Lycopene cyclase%2C'\
                       ' beta and epsilon '},
                      {'stop': u'42291459', 'start': u'42289963',
                       'sca': u'SL2.31ch06', 'name': u'Solyc06g074240.1.1',
                       'desc': u'Lycopene beta-cyclase (AHRD V1 ***- B7U386_ACTCH)%3B'\
                       ' contains Interpro domain(s)  IPR010108  Lycopene cyclase%2C'\
                       ' beta and epsilon '},
                      {'stop': u'60349655', 'start': u'60348153',
                       'sca': u'SL2.31ch10', 'name': u'Solyc10g079480.1.1',
                       'desc': u'Beta-lycopene cyclase (AHRD V1 ***- A6YS01_SOLLC)%3B'\
                       ' contains Interpro domain(s)  IPR010108  Lycopene cyclase%2C'\
                       ' beta and epsilon '},
                      {'stop': u'2291525', 'start': u'2286570',
                       'sca': u'SL2.31ch12', 'name': u'Solyc12g008980.1.1',
                       'desc': u'Lycopene beta cyclase (AHRD V1 **** C1N7E6_MICPS)%3B'\
                       ' contains Interpro domain(s)  IPR010108  Lycopene cyclase%2C'\
                       ' beta and epsilon '}
                      ]}
        output = get_genes_of_proteins(data)
        self.assertEqual(output, expected)

    def test_get_organism_of_proteins(self):
        """ Test the get_organism_of_proteins function ."""
        data = {'key': ['Q38933']}
        expected = {'Q38933': ['Arabidopsis thaliana']}
        output = get_organism_of_proteins(data)
        self.assertEqual(output, expected)

    def test_get_pathways_of_proteins(self):
        """ Test the get_pathways_of_proteins function ."""
        data = {'key': ['Q38933']}
        expected = {'Q38933':
                    ['Carotenoid biosynthesis; beta-carotene biosynthesis.',
                     'Carotenoid biosynthesis; beta-zeacarotene biosynthesis.']
                   }
        output = get_pathways_of_proteins(data)
        self.assertEqual(output, expected)

    def test_get_protein_of_chebi(self):
        """ Test the get_protein_of_chebi function ."""
        data = '17578'
        expected = {'16740': [
                     'http://www.ebi.ac.uk/rhea#rel/controller/UNIPROT:A5W4F2',
                     'http://www.ebi.ac.uk/rhea#rel/controller/UNIPROT:P0C619',
                     'http://www.ebi.ac.uk/rhea#rel/controller/UNIPROT:P0C618',
                     'http://www.ebi.ac.uk/rhea#rel/controller/UNIPROT:A5W4F1',
                             ]}
        output = get_protein_of_chebi(data)
        self.assertEqual(output, expected)


if __name__ == '__main__':
    unittest.main()
