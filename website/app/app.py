#################################################################
#################################################################
############### Dubois RNA-seq Analysis #########################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#######################################################
#######################################################
########## 1. App Configuration
#######################################################
#######################################################

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
### Flask
from flask import Flask, render_template, url_for, request, Response
import os, json, h5py
import pandas as pd

### Nbformat
import nbformat as nbf
from nbconvert.preprocessors import ExecutePreprocessor
from nbconvert import HTMLExporter
from traitlets.config import Config
c = Config()
c.HTMLExporter.preprocessors = ['nbconvert.preprocessors.ExtractOutputPreprocessor']
html_exporter_with_figs = HTMLExporter(config=c)

##### 2. Flask App #####
entry_point = '/l1000ng'
app = Flask(__name__, static_url_path='/app/static')

##### 3. Prefix middleware #####
class PrefixMiddleware(object):

	def __init__(self, app, prefix=''):
		self.app = app
		self.prefix = prefix

	def __call__(self, environ, start_response):
		if environ['PATH_INFO'].startswith(self.prefix):
			environ['PATH_INFO'] = environ['PATH_INFO'][len(self.prefix):]
			environ['SCRIPT_NAME'] = self.prefix
			return self.app(environ, start_response)
		else:
			start_response('404', [('Content-Type', 'text/plain')])
			return ["This url does not belong to the app.".encode()]
app.wsgi_app = PrefixMiddleware(app.wsgi_app, prefix=entry_point)

#######################################################
#######################################################
########## 2. Routes
#######################################################
#######################################################
##### Handles routes used to generate notebooks.

##################################################
########## 2.1 Webpages
##################################################

#############################################
########## 1. Home
#############################################
### Landing page for the website. Links to analyze() endpoint.
### Links to: analyze().

@app.route('/')
def index():

	# Read samples
	return render_template('index.html')

#############################################
########## 2. Analyze
#############################################

@app.route('/analyze')
def analyze():

	# Return JSON
	return render_template('analyze.html')

#############################################
########## 3. Results
#############################################

@app.route('/results', methods=['GET', 'POST'])
def results():

	# Get default
	example_list = ['AML001_CD34_6H:BRD-K43389675:10', 'AML001_PC3_6H:BRD-A19037878:0.37037', 'AML001_PC3_6H:BRD-A19037878:1.11111']

	# Get RIDs
	rids = request.form.get('rid_list', ','.join(example_list))

	# Return JSON
	return render_template('results.html', rids=rids)

#############################################
########## 4. Reusing Notebooks
#############################################

@app.route('/reusing_notebooks')
def reusing_notebooks():

	# RIDs
	rids = ['AML001_CD34_6H:BRD-K43389675:10,AML001_PC3_6H:BRD-A19037878:0.37037,AML001_PC3_6H:BRD-A19037878:1.11111,AML001_PC3_6H:BRD-A19037878:10,AML001_PC3_6H:BRD-A19037878:3.33333']

	# Requirements
	requirements = ['jupyter', 'requests', 'random', 'string', 'os', 'pandas', 'seaborn', 'sklearn', 'plotly']

	return render_template('reusing_notebooks.html', rids=rids, requirements=requirements)

##################################################
########## 2.2 APIs
##################################################

#############################################
########## 1. Metadata API
#############################################

@app.route('/api/metadata')
def metadata_api():
	
	# Read data
	dataframe = pd.read_table('app/static/data/signature_metadata.txt')#.head(50000)
	dataframe['checkbox'] = ''

	# Return JSON
	return json.dumps({'data': dataframe.to_dict(orient='records')})

#############################################
########## 2. Notebook API
#############################################

@app.route('/api/notebook', methods=['GET', 'POST'])
def notebook_api():
	
	# Read data
	if request.method == 'POST':
		request_data = request.form
	elif request.method == 'GET':
		request_data = {'request_type': 'html', 'rids': 'AML001_CD34_6H:BRD-K43389675:10,AML001_PC3_6H:BRD-A19037878:0.37037,AML001_PC3_6H:BRD-A19037878:1.11111,AML001_PC3_6H:BRD-A19037878:10,AML001_PC3_6H:BRD-A19037878:3.33333'}

	# Read notebook
	notebook = nbf.read('app/static/notebooks/l1000ng.ipynb', as_version=4)

	# Replace ids
	notebook['cells'][0]['source'] = notebook['cells'][0]['source'].replace('{rid_list}', str(request_data['rids'].split(',')))

	# Return file
	if request_data['request_type'] == 'ipynb':

		# Convert to string
		notebook_string = nbf.writes(notebook)

		# Return
		return Response(notebook_string, mimetype="text", headers={"Content-disposition": "attachment; filename=l1000ng_notebook.ipynb"})

	# Return notebook
	elif request_data['request_type'] == 'html':
		
		# Initialize preprocess
		print('executing...')
		if os.environ.get('SECRET_KEY'):
			ep = ExecutePreprocessor(timeout=600)#, kernel_name='venv')
		else:
			ep = ExecutePreprocessor(timeout=600, kernel_name='venv')

		# Execute
		ep.preprocess(notebook, {'metadata': {'path': '.'}})

		# Convert
		notebook_html = html_exporter_with_figs.from_notebook_node(notebook)[0]

		# Return HTML
		return notebook_html

#############################################
########## 3. Data API
#############################################

@app.route('/api/data', methods=['POST'])
def data_api():

	# Get RIDs
	rid_list = request.json['rid_list']
	
	# Read h5
	f = h5py.File('app/static/data/signatures.h5', 'r')

	# Get indices
	indices = [index for index, rid in enumerate(f['meta']['sample']['rid']) if rid in rid_list]

	# Get dataframes
	dataframes = {
		'signatures': pd.DataFrame(f['data']['cd'][:, indices], columns=f['meta']['sample']['rid'][indices], index=f['meta']['gene']['pr_gene_symbol']).rename_axis('gene_symbol').astype(float).reset_index(),
		'metadata': pd.DataFrame({key: value.value[indices] for key, value in f['meta']['sample'].items()})
	}

	# Convert to JSON
	result_json = json.dumps({key: value.to_dict(orient='records') for key, value in dataframes.items()})

	# Return
	return result_json
