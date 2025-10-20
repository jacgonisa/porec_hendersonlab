import sys
import matplotlib.pyplot as plt
import re 
import plotly.express as px
from dash import Dash, dcc, html, Input, Output
import numpy as np
#from plotly.subplots import make_subplots
#import plotly.graph_objects as go

plt.rcParams["figure.figsize"] = [16,9]

#enzyme  = sys.argv[1]
fasta_f  = sys.argv[1]
bed_file = sys.argv[2]

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def reverse(dna):
    r_dna = "".join(complement.get(base, base) for base in reversed(dna))
    return r_dna

def read_sequence(fasta_f):
    raw_data = {}
    with open(fasta_f,"r") as f:
        line = f.readline()
        while line:
            if line.count(">") > 0:
                this_name = line.strip().split(">")[1]
                line = f.readline()
                raw_data[this_name] = line.strip()
            line = f.readline()
    return raw_data

def read_bed(bed_file):
    bed_data = {}
    with open(bed_file,"r") as f:
        line = f.readline()
        while line:
            tmp = line.strip().split("\t")
            if tmp[0] in bed_data.keys():
                bed_data[tmp[0]].append(int(tmp[1]))
                bed_data[tmp[0]].append(int(tmp[2]))
            else:
                bed_data[tmp[0]] = [int(tmp[1])]
                bed_data[tmp[0]].append(int(tmp[2]))
            line = f.readline()
    return bed_data

raw_data = read_sequence(fasta_f)
bed_data = read_bed(bed_file)

def find_sites(enzyme_name):
    data = {}
    r_enzyme_name = reverse(enzyme_name)
    for k in raw_data.keys():
        data[k] = []
        for m in re.finditer(enzyme_name,raw_data[k]):
            data[k].append(m.start())
        for m in re.finditer(r_enzyme_name,raw_data[k]):
            data[k].append(m.start())
    return data

bins=1000
std=10
mean = 10
chrom="Chr1"
enzyme = "GATC"
demodropdown = "Chr1"

app = Dash(__name__)


app.layout = html.Div([
    html.H4('Restriction enzyme cut site locations'),
    dcc.Graph(id="graph"),
    dcc.Dropdown(list(raw_data.keys()), "Chr1", id='demodropdown'),
    html.P("Bins:"),
    dcc.Slider(id="bins",min=10, max=168539, value=1000, 
               marks={10: '10', 10000: '10k', 1000:"1k", 5000:"5k", #6000:"6", 7000:"7", 8000:"8", 9000:"9",
               20000: '20k', 30000: '30k', 40000: '40k', 50000: '50k', 60000: '60k', 70000: '70k', 80000: '80k', 90000: '90k', 100000:"100k"}),
    html.P("Enzyme:"),
    dcc.Textarea(
        id='enzy_text',
        value=enzyme,
        style={'width': '50%', 'height': 30},
    ),
])


@app.callback(
    Output("graph", "figure"), 
    Input("bins", "value"), 
    #Input("std", "value"),
    Input('demodropdown', 'value'),
    Input('enzy_text', 'value'))
def display_color(bins, demodropdown, enzy_text):
    #data = np.random.normal(mean, std, size=500) # replace with your own data source
    data = find_sites(enzy_text)
    #r_enzyme = reverse(enzyme)
    #fig = make_subplots(rows=2, cols=1, row_heights=[0.7, 0.3],shared_xaxes=True,)

    fig = px.histogram(data[demodropdown],nbins=bins, title=demodropdown+" ("+str(len(data[demodropdown]))+" total sites)") #, row=1,col=1) #, range_x=[-10, 10])
    if demodropdown in bed_data.keys():
        for b in bed_data[demodropdown]:
            fig.add_vline(x=b,line_dash="dash")
    fig.update_layout(xaxis_title_text='Position in chromosome', yaxis_title_text='# of cut sites')
    return fig

app.run_server(debug=True)

