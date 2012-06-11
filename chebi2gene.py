#!/usr/bin/python

"""
Small RDF API for the ITAG annotation.

We keep the logic simple by querying the sparql endpoint for the desired
information and return its RDF representation.
"""

from flask import Flask, Response, render_template, request, redirect, url_for
from flaskext.wtf import Form, TextField, validators, Required

import datetime
import rdflib
import urllib
import json


"""
Address of the sparql server to query.
"""
SERVER = 'http://sparql.plantbreeding.nl:8080/sparql/'
"""
Turn on or off the debugging mode (turn on only for development).
"""
DEBUG = True
"""
Create the application.
"""
APP = Flask(__name__)
APP.secret_key = 'df;lkhad;fkl234jbcl90-=dasg789flin1234hc kh kkavjnk.djbgf-*iqgfb.vkjb34hrtq2'

"""
Stores in which graphs are the different source of information.
"""
graphs = {
    'uniprot' : 'http://uniprot.pbr.wur.nl/',
    'itag' : 'http://itag2.pbr.wur.nl/',
}


class ChebiIDForm(Form):

    chebi_id = TextField('Chebi ID or molecule name')


def convert_to_uniprot_uri(data):
    """
    
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
    """
    """
    query = '''
    PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
    PREFIX obo:<http://purl.obolibrary.org/obo#>
    SELECT DISTINCT ?id ?name ?syn
    FROM <http://chebi.pbr.wur.nl/>
    WHERE {
      {
        ?id rdfs:label ?name .
        ?id obo:Synonym ?syn .
        FILTER (
            regex(?name, "%(search)s", "i")
        )
      }
    } ORDER BY ?id
    ''' % {'search': name}
    data_js = sparqlQuery(query, 'http://localhost:8890/sparql')
    molecules = {}
    for entry in data_js['results']['bindings']:
        chebi_id = entry['id']['value'].rsplit('/', 1)[1].split('_')[1]
        if chebi_id in molecules:
            molecules[chebi_id]['syn'].append(entry['syn']['value'])
        else:
            molecules[chebi_id] = { 'name' : [entry['name']['value']],
                                    'syn' : [entry['syn']['value']]
                                  }
    return molecules


def get_extended_chebi_from_search(name):
    """
    """
    query = '''
    PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
    PREFIX obo:<http://purl.obolibrary.org/obo#>
    SELECT DISTINCT ?id ?name ?syn
    FROM <http://chebi.pbr.wur.nl/>
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
    ''' % {'search': name}
    data_js = sparqlQuery(query, 'http://localhost:8890/sparql')
    molecules = {}
    for entry in data_js['results']['bindings']:
        chebi_id = entry['id']['value'].rsplit('/', 1)[1].split('_')[1]
        if chebi_id in molecules:
            molecules[chebi_id]['syn'].append(entry['syn']['value'])
        else:
            molecules[chebi_id] = { 'name' : [entry['name']['value']],
                                    'syn' : [entry['syn']['value']]
                                  }
    return molecules
    

def get_genes_of_proteins(data):
    """
    """
    for key in data:
        proteins = data[key]['pathway']
        # Let's make sure the identifiers are unique
        proteins = list(set(proteins))
        query = '''
        PREFIX gene:<http://pbr.wur.nl/GENE#> 
        PREFIX pos:<http://pbr.wur.nl/POSITION#> 
        SELECT DISTINCT ?prot ?name ?sca ?start ?stop ?desc
        FROM <http://itag2.pbr.wur.nl/>
        WHERE{
            ?gene gene:Protein ?prot .
                FILTER (
                ?prot IN (
<http://purl.uniprot.org/uniprot/%s>
                )
            )
            ?gene gene:Position ?pos .
            ?pos pos:Scaffold ?sca .
            ?gene gene:Description ?desc .
            ?gene gene:FeatureName ?name .
            ?pos pos:Start ?start .
            ?pos pos:Stop ?stop .
        } ORDER BY ?name
        ''' % ('>,\n<http://purl.uniprot.org/uniprot/'.join(proteins))
        data_js = sparqlQuery(query, SERVER)
        genes = {}
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

        for protein in proteins:
            if not protein in genes:
                genes[protein] = []
        data[key]['genes'] = genes
    return data


def get_pathways_of_proteins(data):
    """
    """
    output = {}
    for key in data:
        proteins = data[key]
        # Let's make sure the identifiers are unique
        proteins = list(set(proteins))
        query = '''
        PREFIX gene:<http://pbr.wur.nl/GENE#>
        PREFIX uniprot:<http://purl.uniprot.org/core/>
        PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
        SELECT DISTINCT ?prot ?desc
        FROM <http://uniprot.pbr.wur.nl/>
        WHERE {
            ?prot uniprot:annotation ?annot .
            ?annot rdfs:seeAlso ?url .
            ?annot rdfs:comment ?desc .
            FILTER (
                ?prot IN (
<http://purl.uniprot.org/uniprot/%s>
                )
            )
        }
        ''' % ('>,\n<http://purl.uniprot.org/uniprot/'.join(proteins))
        data_js = sparqlQuery(query, SERVER)
        prot = {}
        for entry in data_js['results']['bindings']:
            prot_id = entry['prot']['value'].rsplit('/', 1)[1]
            if prot_id in prot:
                prot[prot_id].append(entry['desc']['value'])
            else:
                prot[prot_id] = [entry['desc']['value']]
        for protein in proteins:
            if not protein in prot:
                prot[protein] = []
        output[key] = {'pathway': prot}
    return output


def get_protein_of_chebi(chebi_id):
    """
    """
    query = '''
    prefix bp: <http://www.biopax.org/release/biopax-level2.owl#>
    SELECT DISTINCT ?react ?xref
    FROM <http://rhea.pbr.wur.nl/>
    WHERE {
      ?cmp bp:XREF <http://www.ebi.ac.uk/rhea#CHEBI:%s> .
      ?dir ?p ?cmp .
      ?react ?p2 ?dir .
      ?react bp:XREF ?xref .
      FILTER (
        regex(?xref, 'UNIPROT')
      )
    }
    ''' % chebi_id
    data = sparqlQuery(query, SERVER)
    output = {}
    for entry in data['results']['bindings']:
        key = entry['react']['value'].split('#')[1]
        if key in output:
            output[key].append(entry['xref']['value'])
        else:
            output[key] = [entry['xref']['value']]
    return output


def sparqlQuery(query, server, format = 'application/json'):
    """ Runs the given SPARQL query against the desired sparql endpoint
    and return the output in the format asked (default being rdf/xml).

    :param query, the string of the sparql query that should be ran.
    :param server, a string, the url of the sparql endpoint that we want
    to run query against.
    :param format, specifies in which format we want to have the output.
    Defaults to `application/json` but can also be `application/rdf+xml`.
    """
    params = {
        'default-graph': '',
        'should-sponge': 'soft',
        'query': query,
        'debug': 'off',
        'timeout': '',
        'format': format,
        'save': 'display',
        'fname': ''
    }
    querypart = urllib.urlencode(params)
    response = urllib.urlopen(server, querypart).read()
    return json.loads(response)


def run_query_via_rdflib(query, server):
    """ Runs the given query of the given server, loads the results
    rdf/xml into a rdflib.Graph and return a rdf/xml representation of
    this graph.
    This is a bit of a hack to return a nicer rdf/xml representation of
    the knowledge retrieve than older version of virtuoso offers.
    From version 6.1.5 at least, this trick should not be needed
    anymore.

    :param query, the string of the sparql query that should be ran.
    :param server, a string, the url of the sparql endpoint that we want
    to run query against.
    """
    graph = rdflib.Graph()
    graph.parse(data=sparqlQuery(query, server), format="application/rdf+xml")
    return graph.serialize(format='xml')


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
        except ValueError, er:
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
    if len(molecules) == 1:
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
    proteins = convert_to_uniprot_uri(proteins)
    pathways = get_pathways_of_proteins(proteins)
    genes = get_genes_of_proteins(pathways)
    return render_template('output.html', data=genes, chebi=chebi_id)

@APP.route('/csv/<chebi_id>')
def generate_csv(chebi_id):
    """ Generate a comma separated value file containing all the
    information.
    """
    print 'Chebi2gene %s -- %s -- %s' % (datetime.datetime.now(),
        request.remote_addr, request.url)
    # Regenerate the informations
    proteins = get_protein_of_chebi(chebi_id)
    proteins = convert_to_uniprot_uri(proteins)
    pathways = get_pathways_of_proteins(proteins)
    data = get_genes_of_proteins(pathways)

    string = 'Chebi ID, Chebi URL, Rhea ID, Rhea URL, UniProt \
    URL, Type, Name, Scaffold, Start, Stop, Description\n'
    chebi_url = 'http://www.ebi.ac.uk/chebi/searchId.do?chebiId=%s' % \
        chebi_id
    for reactions in data:
        react_url = 'http://www.ebi.ac.uk/rhea/reaction.xhtml?id=RHEA:%s' % \
            reactions
        if 'pathway' in data[reactions]:
            for proteins in data[reactions]['pathway']:
                for pathways in data[reactions]['pathway'][proteins]:
                    string = string + '%s,%s,%s,%s,%s,Pathway,%s\n' % (
                        chebi_id, chebi_url, reactions, react_url, proteins,
                        pathways)
        if 'genes' in data[reactions]:
            for proteins in data[reactions]['genes']:
                for gene in data[reactions]['genes'][proteins]:
                    string = string + '%s,%s,%s,%s,%s,Gene,%s,%s,%s,%s,%s\n' % (
                        chebi_id, chebi_url, reactions, react_url, proteins,
                        gene['name'], gene['sca'],
                        gene['start'], gene['stop'], gene['desc'])
    return Response(string, mimetype='application/excel')


if __name__ == '__main__':
    APP.debug = DEBUG
    APP.run()
