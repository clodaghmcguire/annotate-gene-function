#import libraries
import requests
import xmltodict
import json
import urllib.request
from urllib.request import urlopen, HTTPError

#read in csv of genes of interest
#must include a column 'Gene stable ID' that contains ensembl ID
#Ensembl IDs can be downloaded from BioMart
ensembl_ids = pd.read_csv('ensemble_geneIDs.csv')
all_gene_ids = ensembl_ids['Gene stable ID'].tolist()

gene_annotations = []
for gene_id in all_gene_ids:
    gene_dict = {}
    gene_dict['gene_id'] = gene_id
    gene_dict['gene_name'] = ""
    try:
        url = "https://mygene.info/v3/gene/"+gene_id
        response = urllib.request.urlopen(url)
        source = response.read()
        get_function = json.loads(source)
        molecular_function = []
        try:
            for term in get_function['go']['MF']:
                go_term = term['term']
                go_id = term['id']
                if go_id +": " + go_term not in molecular_function:
                    molecular_function.append(go_id +": " + go_term)
            gene_dict['molecular_function'] = ', \n'.join(molecular_function)
        except:
            gene_dict['molecular_function'] = ""
        biological_process = []
        try:
            for term in get_function['go']['BP']:
                go_term = term['term']
                go_id = term['id']
                if go_id +": " + go_term not in biological_process:
                    biological_process.append(go_id +": " + go_term)
            gene_dict['biological_process'] = ', \n'.join(biological_process)
        except:
            gene_dict['biological_process'] = ""
        molecular_pathways = []
        try:
            for pathway in get_function['pathway']['reactome']:
                if pathway['name'] not in molecular_pathways:
                    molecular_pathways.append(pathway['name'])
            gene_dict['reactome_pathways'] = ', \n'.join(molecular_pathways)
        except:
            gene_dict['reactome_pathways'] = ""
        wikipathways = []
        try:
            for pathway in get_function['pathway']['wikipathways']:
                if pathway['name'] not in wikipathways:
                    wikipathways.append(pathway['name'])
            gene_dict['wikipathways'] = ', \n'.join(wikipathways)
        except:
            gene_dict['wikipathways'] = ""
    except HTTPError:
        gene_dict['molecular_function'] = ""
        gene_dict['biological_process'] = ""
        gene_dict['reactome_pathways'] = ""
        gene_dict['wikipathways'] = ""
    
    try:
        url = "https://www.proteinatlas.org/"+gene_id+".xml"
        response = urlopen(url)
        content = response.read()
        xml_to_dict = xmltodict.parse(content)
        get_expression = json.dumps(xml_to_dict, indent=4)
        data = json.loads(get_expression)
        gene_dict['gene_name'] = data['proteinAtlas']['entry']['name']
        expression_data = data['proteinAtlas']['entry']['tissueExpression']['data']
        tissue_expression = []
        for data in expression_data:
            #print(data)
            tissue = data['tissue']['#text']
            if data['level']['#text'] == "not detected":
                #print(tissue)
                continue
            else:
                level = data['level']['#text']
                #print(tissue)
                #print(level)
                tissue_expression.append(tissue +": "+ level)
        gene_dict['expression'] = ', \n'.join(tissue_expression)
    except:
        gene_dict['expression'] = ""
    
    gene_annotations.append(gene_dict)

df = pd.DataFrame.from_dict(gene_annotations)
df.to_csv("gene_annotations.csv", index=False)
