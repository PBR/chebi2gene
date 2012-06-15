#!/usr/bin/python

"""
Small web application to retrieve information from uniprot and itag for
a given compound.

The idea is that for one compound we are able to find out in which
reactions it is involved and what are the proteins involved in these
reactions. For each of these proteins we can find if there are genes and
genes from tomato associated with them.
"""

from flask import Flask, Response, render_template, request, redirect, url_for
from flaskext.wtf import Form, TextField

import datetime
import rdflib
import urllib
import json


# Address of the sparql server to query.
SERVER = 'http://sparql.plantbreeding.nl:8080/sparql/'
# Turn on or off the debugging mode (turn on only for development).
DEBUG = True
# Create the application.
APP = Flask(__name__)
APP.secret_key = 'df;lkhad;fkl234jbcl90-=davjnk.djbgf-*iqgfb.vkjb34hrt' \
'q2lkhflkjdhflkdjhbfakljgipfurp923243nmrlenr;k3jbt;kt'

# Stores in which graphs are the different source of information.
GRAPHS = {
    'uniprot': 'http://uniprot.pbr.wur.nl/',
    'itag': 'http://itag2.pbr.wur.nl/',
    'chebi': 'http://chebi.pbr.wur.nl/',
    'rhea': 'http://rhea.pbr.wur.nl/',
}


class ChebiIDForm(Form):
    """ Simple text field form to input the chebi identifier or the
    name of the protein.
    """
    chebi_id = TextField('Chebi ID or molecule name')


def convert_to_uniprot_id(data):
    """ Converts from RHEA Uniprot URI to Uniprot ID.

    @param data, a dictionary of String: [String] where the keys are
    reaction ID and the values are protein URI.
    @return, a dictionary of String: [String] where the keys are
    reaction ID and the values are protein ID.
    """
    for key in data:
        proteins = data[key]
        proteins2 = []
        for protein in proteins:
            prot_id = protein.rsplit(':', 1)[1]
            proteins2.append(prot_id.strip())
        data[key] = proteins2
    return data


def get_exact_chebi_from_search(name):
    """ Search the chebi database for molecule having the given string
    in their name. The data returned contains the chebi identifier, the
    name and synonyms of the molecule in chebi.

    @param name, a string, name of the molecule to search in chebi.
    @return, a dictionary containing all the molecule found for having
    the input string in their name. The data structure returned is like:
    {string: {'name': string, 'syn': [String]}}, where the keys are the
    chebi identifier and the values are dictionaries containing the
    name of the molecules and a list of its synonym.
    """
    query = '''
    PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
    PREFIX obo:<http://purl.obolibrary.org/obo#>
    SELECT DISTINCT ?id ?name ?syn
    FROM <%(chebi)s>
    WHERE {
      {
        ?id rdfs:label ?name .
        ?id obo:Synonym ?syn .
        FILTER (
            regex(?name, "%(search)s", "i")
        )
      }
    } ORDER BY ?id
    ''' % {'search': name, 'chebi': GRAPHS['chebi']}
    data_js = sparql_query(query, 'http://localhost:8890/sparql')
    if not data_js:
        return
    molecules = {}
    for entry in data_js['results']['bindings']:
        chebi_id = entry['id']['value'].rsplit('/', 1)[1].split('_')[1]
        if chebi_id in molecules:
            molecules[chebi_id]['syn'].append(entry['syn']['value'])
        else:
            molecules[chebi_id] = {
                                    'name': [entry['name']['value']],
                                    'syn': [entry['syn']['value']]
                                  }
    return molecules


def get_extended_chebi_from_search(name):
    """ Search the chebi database for molecule having the given string
    in their name or in their synonyms. The data returned contains the
    chebi identifier, the name and synonyms of the molecule in chebi.

    @param name, a string, name of the molecule to search in chebi.
    @return, a dictionary containing all the molecule found for having
    the input string in their name or in their synonyms.
    The data structure returned is like:
    {string: {'name': string, 'syn': [String]}}, where the keys are the
    chebi identifier and the values are dictionaries containing the
    name of the molecules and a list of its synonym.
    """
    query = '''
    PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
    PREFIX obo:<http://purl.obolibrary.org/obo#>
    SELECT DISTINCT ?id ?name ?syn
    FROM <%(chebi)s>
    WHERE {
      {
        ?id rdfs:label ?name .
        ?id obo:Synonym ?syn .
        FILTER (
            regex(?name, "%(search)s", "i")
            || regex(?syn, "%(search)s", "i")
        )
      }
    } ORDER BY ?id
    ''' % {'search': name, 'chebi': GRAPHS['chebi']}
    data_js = sparql_query(query, 'http://localhost:8890/sparql')
    if not data_js:
        return
    molecules = {}
    for entry in data_js['results']['bindings']:
        chebi_id = entry['id']['value'].rsplit('/', 1)[1].split('_')[1]
        if chebi_id in molecules:
            molecules[chebi_id]['syn'].append(entry['syn']['value'])
        else:
            molecules[chebi_id] = {
                                    'name': [entry['name']['value']],
                                    'syn': [entry['syn']['value']]
                                  }
    return molecules


def get_genes_of_proteins(data):
    """ Returns the genes associated with proteins.

    @param name, a dictionary where the keys are reactions identifier
    and the values lists of proteins.
    @return, a dictionary containing all the genes related with the
    proteins specified.
    The data structure returned is like:
    {string: [String]}, where the keys are the uniprot identifier and
    the values are list of gene identifier associated with the protein.
    """
    genes = {}
    for key in data:
        proteins = data[key]
        # Let's make sure the identifiers are unique
        proteins = list(set(proteins))
        query = '''
        PREFIX gene:<http://pbr.wur.nl/GENE#>
        PREFIX pos:<http://pbr.wur.nl/POSITION#>
        SELECT DISTINCT ?prot ?name ?sca ?start ?stop ?desc
        FROM <%(itag)s>
        WHERE{
            ?gene gene:Protein ?prot .
                FILTER (
                ?prot IN (
<http://purl.uniprot.org/uniprot/%(prot)s>
                )
            )
            ?gene gene:Position ?pos .
            ?pos pos:Scaffold ?sca .
            ?gene gene:Description ?desc .
            ?gene gene:FeatureName ?name .
            ?pos pos:Start ?start .
            ?pos pos:Stop ?stop .
        } ORDER BY ?name
        ''' % {'prot': '>,\n<http://purl.uniprot.org/uniprot/'.join(
            proteins), 'itag': GRAPHS['itag']}
        data_js = sparql_query(query, SERVER)
        for entry in data_js['results']['bindings']:
            prot_id = entry['prot']['value'].rsplit('/', 1)[1]
            gene = {}
            for var in ['name', 'sca', 'start', 'stop', 'desc']:
                gene[var] = entry[var]['value']

            gene['sca'] = gene['sca'].rsplit('#', 1)[1]
            if prot_id in genes:
                genes[prot_id].append(gene)
            else:
                genes[prot_id] = [gene]

    return genes


def get_pathways_of_proteins(data):
    """ Returns the pathways associated with proteins.

    @param name, a dictionary where the keys are reactions identifier
    and the values lists of proteins.
    @return, a dictionary containing all the pathways related with the
    proteins specified.
    The data structure returned is like:
    {string: [String]}, where the keys are the uniprot identifier and
    the values are list of pathways associated with the protein.
    """
    pathways = {}
    for key in data:
        proteins = data[key]
        # Let's make sure the identifiers are unique
        proteins = list(set(proteins))
        query = '''
        PREFIX gene:<http://pbr.wur.nl/GENE#>
        PREFIX uniprot:<http://purl.uniprot.org/core/>
        PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
        SELECT DISTINCT ?prot ?desc
        FROM <%(uniprot)s>
        WHERE {
            ?prot uniprot:annotation ?annot .
            ?annot rdfs:seeAlso ?url .
            ?annot rdfs:comment ?desc .
            FILTER (
                ?prot IN (
<http://purl.uniprot.org/uniprot/%(prot)s>
                )
            )
        }
        ''' % {'prot':
        '>,\n<http://purl.uniprot.org/uniprot/'.join(proteins),
                'uniprot': GRAPHS['uniprot']}
        data_js = sparql_query(query, SERVER)
        for entry in data_js['results']['bindings']:
            prot_id = entry['prot']['value'].rsplit('/', 1)[1]
            path = entry['desc']['value']
            if prot_id in pathways and path not in pathways[prot_id]:
                pathways[prot_id].append(path)
            else:
                pathways[prot_id] = [path]
    return pathways


def get_organism_of_proteins(data):
    """ Returns the all organism associated with the proteins.

    @param name, a dictionary where the keys are reactions identifier
    and the values lists of proteins.
    @return, a dictionary containing all the organism related with the
    proteins specified.
    The data structure returned is like:
    {string: [String]}, where the keys are the uniprot identifier and
    the values are list of organisms associated with the protein.
    """
    organism = {}
    for key in data:
        proteins = data[key]
        # Let's make sure the identifiers are unique
        proteins = list(set(proteins))
        query = '''
        PREFIX uniprot:<http://purl.uniprot.org/core/>
        SELECT DISTINCT ?prot ?name
        FROM <%(uniprot)s>
        WHERE {
            ?prot uniprot:organism ?orga .
            ?orga uniprot:scientificName ?name .
            FILTER (
                ?prot IN (
<http://purl.uniprot.org/uniprot/%(prot)s>
                )
            )
        }
        ''' % {'prot':
            '>,\n<http://purl.uniprot.org/uniprot/'.join(proteins),
                'uniprot': GRAPHS['uniprot']}
        data_js = sparql_query(query, SERVER)
        for entry in data_js['results']['bindings']:
            prot_id = entry['prot']['value'].rsplit('/', 1)[1]
            orga = entry['name']['value']
            if prot_id in organism and orga not in organism[prot_id]:
                organism[prot_id].append(orga)
            else:
                organism[prot_id] = [orga]
    return organism


def get_protein_of_chebi(chebi_id):
    """ Returns the all protein associated with a compound.

    @param name, a string, identifier of a compound on chebi.
    @return, a dictionary containing all the proteins related with the
    compound specified.
    The data structure returned is like:
    {string: [String]}, where the keys are reaction identifiers and the
    values are list of proteins associated with the reaction.
    """
    query = '''
    prefix bp: <http://www.biopax.org/release/biopax-level2.owl#>
    SELECT DISTINCT ?react ?xref
    FROM <%(rhea)s>
    WHERE {
      ?cmp bp:XREF <http://www.ebi.ac.uk/rhea#CHEBI:%(chebi_id)s> .
      ?dir ?p ?cmp .
      ?react ?p2 ?dir .
      ?react bp:XREF ?xref .
      FILTER (
        regex(?xref, 'UNIPROT')
      )
    }
    ''' % {'chebi_id': chebi_id, 'rhea': GRAPHS['rhea']}
    data = sparql_query(query, SERVER)
    if not data:
        return
    output = {}
    for entry in data['results']['bindings']:
        key = entry['react']['value'].split('#')[1]
        if key in output:
            output[key].append(entry['xref']['value'])
        else:
            output[key] = [entry['xref']['value']]
    return output


def sparql_query(query, server, output_format='application/json'):
    """ Runs the given SPARQL query against the desired sparql endpoint
    and return the output in the format asked (default being rdf/xml).

    @param query, the string of the sparql query that should be ran.
    @param server, a string, the url of the sparql endpoint that we want
    to run query against.
    @param format, specifies in which format we want to have the output.
    Defaults to `application/json` but can also be `application/rdf+xml`.
    @return, a JSON object, representing the output of the provided
    sparql query.
    """
    params = {
        'default-graph': '',
        'should-sponge': 'soft',
        'query': query,
        'debug': 'off',
        'timeout': '',
        'format': output_format,
        'save': 'display',
        'fname': ''
    }
    querypart = urllib.urlencode(params)
    response = urllib.urlopen(server, querypart).read()
    try:
        output = json.loads(response)
    except ValueError:
        output = {}
    return output


def run_query_via_rdflib(query, server):
    """ Runs the given query of the given server, loads the results
    rdf/xml into a rdflib.Graph and return a rdf/xml representation of
    this graph.
    This is a bit of a hack to return a nicer rdf/xml representation of
    the knowledge retrieve than older version of virtuoso offers.
    From version 6.1.5 at least, this trick should not be needed
    anymore.

    @param query, the string of the sparql query that should be ran.
    @param server, a string, the url of the sparql endpoint that we want
    to run query against.
    @return, a string, representing the rdf output of the provided query.
    """
    graph = rdflib.Graph()
    graph.parse(data=sparql_query(query, server),
        output_format="application/rdf+xml")
    return graph.serialize(format='xml')


##  Web-app


@APP.route('/', methods=['GET', 'POST'])
def index():
    """ Shows the front page.
    All the content of this page is in the index.html file under the
    templates directory. The file is full html and has no templating
    logic within.
    """
    print 'Chebi2gene %s -- %s -- %s' % (datetime.datetime.now(),
        request.remote_addr, request.url)
    form = ChebiIDForm(csrf_enabled=False)
    if form.validate_on_submit():
        try:
            int(form.chebi_id.data)
            return redirect(url_for('show_chebi',
                chebi_id=form.chebi_id.data))
        except ValueError:
            return redirect(url_for('search_chebi',
                name=form.chebi_id.data))
    return render_template('index.html', form=form)


@APP.route('/search/<name>')
def search_chebi(name):
    """ Search the CHEBI database for the name given.
    """
    print 'Chebi2gene %s -- %s -- %s' % (datetime.datetime.now(),
        request.remote_addr, request.url)
    molecules = get_exact_chebi_from_search(name)
    if molecules and len(molecules) == 1:
        return redirect(url_for('show_chebi',
                chebi_id=molecules.keys()[0]))
    return render_template('search.html', data=molecules, search=name,
        extended=False)


@APP.route('/fullsearch/<name>')
def search_chebi_extended(name):
    """ Search the CHEBI database for the name given including the
    synonyms.
    """
    print 'Chebi2gene %s -- %s -- %s' % (datetime.datetime.now(),
        request.remote_addr, request.url)
    molecules = get_extended_chebi_from_search(name)
    return render_template('search.html', data=molecules, search=name,
        extended=True)


@APP.route('/chebi/<chebi_id>')
def show_chebi(chebi_id):
    """ Shows the front page.
    All the content of this page is in the index.html file under the
    templates directory. The file is full html and has no templating
    logic within.
    """
    print 'Chebi2gene %s -- %s -- %s' % (datetime.datetime.now(),
        request.remote_addr, request.url)
    proteins = get_protein_of_chebi(chebi_id)
    if not proteins:
        return render_template('output.html', proteins=[],
        pathways=None, genes=None, organisms=None, chebi=chebi_id)
    proteins = convert_to_uniprot_id(proteins)
    pathways = get_pathways_of_proteins(proteins)
    genes = get_genes_of_proteins(proteins)
    organisms = get_organism_of_proteins(proteins)
    return render_template('output.html', proteins=proteins,
        pathways=pathways, genes=genes, organisms=organisms,
        chebi=chebi_id)


@APP.route('/csv/<chebi_id>')
def generate_csv(chebi_id):
    """ Generate a comma separated value file containing all the
    information.
    """
    print 'Chebi2gene %s -- %s -- %s' % (datetime.datetime.now(),
        request.remote_addr, request.url)
    # Regenerate the informations
    proteins = get_protein_of_chebi(chebi_id)
    proteins = convert_to_uniprot_id(proteins)
    pathways = get_pathways_of_proteins(proteins)
    genes = get_genes_of_proteins(proteins)
    organisms = get_organism_of_proteins(proteins)

    string = 'Chebi ID, Chebi URL, Rhea ID, Rhea URL, UniProt, \
Organism, Type, Name, Scaffold, Start, Stop, Description\n'
    chebi_url = 'http://www.ebi.ac.uk/chebi/searchId.do?chebiId=%s' % \
        chebi_id
    for reaction in proteins:
        react_url = 'http://www.ebi.ac.uk/rhea/reaction.xhtml?id=RHEA:%s' % \
            reaction
        for protein in proteins[reaction]:
            if protein in pathways:
                for pathway in pathways[protein]:
                    string = string + '%s,%s,%s,%s,%s,%s,Pathway,%s\n' % (
                        chebi_id, chebi_url, reaction, react_url, protein,
                        " - ".join(organisms[protein]),
                        pathway)
            if protein in genes:
                for gene in genes[protein]:
                    string = string + \
                    '%s,%s,%s,%s,%s,%s,Gene,%s,%s,%s,%s,%s\n' % (
                        chebi_id, chebi_url, reaction, react_url, protein,
                        " - ".join(organisms[protein]),
                        gene['name'], gene['sca'],
                        gene['start'], gene['stop'], gene['desc'])
    return Response(string, mimetype='application/excel')


if __name__ == '__main__':
    APP.debug = DEBUG
    APP.run()
