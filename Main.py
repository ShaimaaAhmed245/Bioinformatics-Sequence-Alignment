import os, io, random
import string
import numpy as np
import pandas as pd

import subprocess

from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from io import StringIO, BytesIO
from Bio.Align.Applications import MuscleCommandline

import panel as pn
import panel.widgets as pnw
pn.extension()

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models import CustomJS, TextInput, Paragraph
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot
from bokeh.io import curdoc

from bokeh.io import *
from bokeh.io import show
from bokeh.plotting import figure, output_file, show
from bokeh.models import *
from bokeh.events import *
from bokeh.models.widgets import FileInput
from bokeh.palettes import YlGn9

from pybase64 import b64decode

from GlobalAlignment import *
from LocalAlignment import *
from MSAMetrics import *

def get_colors(seqs):
    # make colors for bases in sequence
    # TAGC- --> green,red,orange,blue,white
    text = [i for s in list(seqs) for i in s]
    clrs =  {'A':'red','T':'green','G':'orange','C':'blue','-':'white'}
    colors = [clrs[i] for i in text]
    return colors
 
# Create inital figure 
p = figure(x_range=(0, 100),
           plot_height=250,
           plot_width=1000,
           title="Sequence Alignment",
           tools="xpan,reset",
           y_range=(0, 100) 
           )
# for letters of sequence to be appeared on the figure (A,C,G,T,-) 
glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
            text_font="monospace",text_font_size="9pt")
# for colored rectangle of sequence
rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
            line_color=None, fill_alpha=0.4)
# data for the changed values of the figure to be updated when changed
source = ColumnDataSource(dict(x=[], y=[], recty=[], text=[], colors=[]))
p.add_glyph(source, glyph)
p.add_glyph(source, rects)

# Create inital figure for matrix
short_labels = ["Seq1", "Seq2", "Seq3", "Seq4", "Seq5"]
cm = np.array([[0., 0., 0., 0., 0.],
               [0., 0., 0., 0., 0.],
               [0., 0., 0., 0., 0.],
               [0., 0., 0., 0., 0.],
               [0., 0., 0., 0., 0.]])

cm_df = pd.DataFrame(cm, index=short_labels, columns=short_labels)
cm_df = pd.DataFrame(cm_df.stack(),columns=['val']).reset_index()
cm_source = ColumnDataSource(cm_df)
mapper = LinearColorMapper(palette="Viridis256")

met = figure(
    title=f"Metrics",
    toolbar_location=None,
    x_range=short_labels,
    y_range=list(reversed(short_labels)),
    )

met.rect(
    source=cm_source,
    x='level_1',
    y='level_0',
    width=1, 
    height=1, 
    fill_color={'field':'val','transform': mapper},
    )

met.text(
    source=cm_source,
    x='level_1',
    y='level_0',
    text='val',
    text_font_size='10pt',
    x_offset=-30,
    y_offset=10
    )

# Bokeh sequence alignment view
def view_alignment(aln, fontsize="9pt", plot_width=800):
    # make sequence and id lists  and text from the aln object
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]    
    text = [i for s in list(seqs) for i in s]
    colors = get_colors(seqs)    
    N = len(seqs[0])
    S = len(seqs)    
    width = .4

    x = np.arange(1,N+1)
    y = np.arange(0,S,1)
    #creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    #flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    # use recty for rect coords with an offset
    recty = gy+.5
    # now we can update the ColumnDataSource with all the arrays
    source.data.update([('x', gx),('y', gy),('recty', recty),('text', text),('colors', colors)])
    # id of the sequence 
    Dict = {}
    for i in range(0,len(ids)):
        Dict[i+1] = ids[i]

    # plot height update
    plot_height = len(seqs)*15+50
    p.plot_height= plot_height

    if N>100:
        viewlen=100
    else:
        viewlen=N
    #view_range is for the close up view
    # view_range = (0,viewlen)
    p.x_range.end = viewlen
    p.y_range.end = S

    p.grid.visible = False
    p.xaxis.major_label_text_font_style = "bold"
    p.yaxis.minor_tick_line_width = 0
    p.yaxis.major_tick_line_width = 0
    p.yaxis.ticker = list(range(1,len(ids)+1))
    p.yaxis.major_label_overrides = Dict


def upload_fit_data(attr, old, new):
    print("fit data upload succeeded")
    # decode data to be read
    decoded = b64decode(new)
    # not changed binary another decode will be applied
    decoded = decoded.decode("utf-8")
    records = SeqIO.parse(StringIO(decoded), "fasta")
    count = 0
    for rec in records:
        if count == 0:
            text_input_seq1.value = str(rec.seq)
            count += 1
        elif count == 1:
            text_input_seq2.value = str(rec.seq)
            count += 1
        elif count == 2:
            text_input_seq3.value = str(rec.seq)
            count += 1
        elif count == 3:
            text_input_seq4.value = str(rec.seq)
            count += 1
# browse button -->  only fasta files
file_input = FileInput(accept=".fasta")
# call the callback function (upload_fit_data) when choose file
file_input.on_change('value', upload_fit_data)

def update_metrics():
    mat1 = np.zeros((5,5)) 
    mat_df = pd.DataFrame(mat1, index=short_labels, columns=short_labels)
    mat_df = pd.DataFrame(mat_df.stack(),columns=['val']).reset_index()
    cm_source.data.update(mat_df)

#get record of sequence from text input 
def get_record(act):
    if act == 0:
        x = text_input_seq1.value
    elif act == 1:
        x = text_input_seq2.value
    elif act == 2:
        x = text_input_seq3.value
    elif act == 3:
        x = text_input_seq4.value
    elif act == 4:
        x = text_input_seq5.value
    return x

def is_empty():
    count = 0
    if text_input_seq1.value == '':
        count += 1
    if text_input_seq2.value == '':
        count += 1
    if text_input_seq3.value == '':
        count += 1
    if text_input_seq4.value == '':
        count += 1
    if text_input_seq5.value == '':
        count += 1
    return count
    # if count == 5:
    #     return True
    # else:
    #     return False

# checking if the seq contains any letters rather than (A,T,C,G)
def check_letters(seq):
    count = 0
    for i in range(0,len(seq)):
        if seq[i] != 'A' and seq[i] != 'G' and seq[i] != 'T' and seq[i] != 'C':
            count += 1
    
    if count == 0:
        return True
    else:
        return False

# check gap, match, mismatch digits and not empty 
def check_MGM():
    flag = False
    if text_input_gap.value == " " or len(text_input_gap.value) == 0 or text_input_gap.value.lstrip('-+').isdigit() == False:
        flag = True
    if text_input_match.value == " " or len(text_input_match.value) == 0 or text_input_match.value.lstrip('-+').isdigit() == False:
        flag = True
    if text_input_mismatch.value == " " or len(text_input_mismatch.value) == 0 or text_input_mismatch.value.lstrip('-+').isdigit() == False:
        flag = True
    return flag

def global_alignment():
    if (len(checkbox_button_group.active) == 2 and is_empty() <= 3 and check_MGM() == False):
        callback_holder.visible = False
        callback_holder.text = ""
        # pass the match, mismatch and gap value to the global aligment class
        p1 = GlobalAlignment(int(text_input_gap.value), int(text_input_match.value), int(text_input_mismatch.value))
        seq1 = get_record(checkbox_button_group.active[0]).upper()
        seq2 = get_record(checkbox_button_group.active[1]).upper()
        # check letter (A,G,C,T) only # len not equal zero (not empty)
        if check_letters(seq1) == True and check_letters(seq2) == True and len(seq1) != 0 and len(seq2) != 0:
            # call function to align the 2 sequences in the global alignment class
            # return the score and the 2 sequnces aligned
            score, xSeq, ySeq= p1.getAlignedSequences(seq1,seq2)
            # set score value
            text_ouput_seqScore.value = str(score)
            # A A --> AA
            xSeq = ''.join(xSeq)
            ySeq = ''.join(ySeq)
            # create new list ofsequene records with the aligned sequences to be represented
            seq_records = [SeqRecord(Seq(xSeq),id="Seq" + str(checkbox_button_group.active[0]+1)), SeqRecord(Seq(ySeq),id="Seq"+ str(checkbox_button_group.active[1]+1))]
            # representation of the seuences
            view_alignment(seq_records, plot_width=900)
            update_metrics()
        else:
            callback_holder.visible = True
            callback_holder.text = "Your letters should be A G T C"       
    else:
        if check_MGM() == True:
            callback_holder.visible = True
            callback_holder.text = "Please write Gap, Match, Mismatch for global alignment"
        elif is_empty() > 3:
            callback_holder.visible = True
            callback_holder.text = "Please write yor sequence for global alignment"
        else:
            callback_holder.visible = True
            callback_holder.text = "Choose 2 sequences only for global alignment"

def local_alignment():
    if (len(checkbox_button_group.active) == 2 and is_empty() <= 3 and check_MGM() == False):
        callback_holder.visible = False
        callback_holder.text = ""
        # pass the match, mismatch and gap value to the local aligment class
        p1 = LocalAlignment(int(text_input_gap.value), int(text_input_match.value), int(text_input_mismatch.value))
        x = get_record(checkbox_button_group.active[0]).upper()
        y = get_record(checkbox_button_group.active[1]).upper()
        # check letter (A,G,C,T) only # len not equal zero (not empty)
        if check_letters(x) == True and check_letters(y) == True and len(x) != 0 and len(y) != 0:
            # call function to align the 2 sequences in the local alignment class
            # return the score and the 2 sequnces aligned
            bestScore,xSeq, ySeq = p1.getSequence(x,y)
            text_ouput_seqScore.value = str(bestScore)
            xSeq = ''.join(xSeq)
            ySeq = ''.join(ySeq)
            # create new list ofsequene records with the aligned sequences to be represented
            seq_records = [SeqRecord(Seq(xSeq),id="Seq" + str(checkbox_button_group.active[0]+1)), SeqRecord(Seq(ySeq),id="Seq"+ str(checkbox_button_group.active[1]+1))]
            view_alignment(seq_records, plot_width=900)
            update_metrics()
        else:
            callback_holder.visible = True
            callback_holder.text = "Your letters should be A G T C"
    else:
        if check_MGM() == True:
            callback_holder.visible = True
            callback_holder.text = "Please write Gap, Match, Mismatch for local alignment"
        elif is_empty() > 3:
            callback_holder.visible = True
            callback_holder.text = "Please write yor sequence for local alignment"
        else:
            callback_holder.visible = True
            callback_holder.text = "Choose 2 sequences only for local alignment"

# check on the clicked buttons sequence should be greater than 2 and sequence not empty and MGM not empty (match,mismatch and gap)
def multiple_alignment():
    if (len(checkbox_button_group.active) > 2 and len(checkbox_button_group.active) <= 5 and is_empty() <=2 and check_MGM() == False):
        # error message disapper and empty
        callback_holder.visible = False
        callback_holder.text = ""
        seqs = []
        flag = 0
        for i in range (0,len(checkbox_button_group.active)):
            seq = get_record(checkbox_button_group.active[i]).upper()
            # if there is an issue in letters
            if check_letters(seq) == False:
                flag = 1
            # sequence is empty
            if len(seq) == 0:
                flag = 1
            seq = SeqRecord(Seq(seq),id="Seq"+str((checkbox_button_group.active[i]+1)))
            seqs.append(seq)
        if flag == 0:
            # write new fasta file with the input sequences by user
            filename = 'MSAseq.fasta'
            SeqIO.write(seqs, filename, "fasta")
            # pass the new fasta file created to muscle (align) to be aligned
            # the output will be the MSA of the sequences (output)
            output = subprocess.check_output(["C:\\Users\\ADMIN\Downloads\\Bio_Project\\muscle5.1.win64.exe",
            "-align", r"C:\Users\ADMIN\Downloads\Bio_Project\Bio_Project\MSAseq.fasta",
            "-output", r"C:\Users\ADMIN\Downloads\Bio_Project\Bio_Project\MSAseqOut.fasta"],
            text=True)
            # read the MSA output file to show the results of the alignment in the figure
            seq_records = list(SeqIO.parse(r"C:\Users\ADMIN\Downloads\Bio_Project\Bio_Project\MSAseqOut.fasta", "fasta"))
            view_alignment(seq_records, plot_width=900)
            if dropdown.value == "Sum of Pairs":
                similarity,matrix = sum_of_pairs(seq_records, int(text_input_gap.value), int(text_input_match.value), int(text_input_mismatch.value))
            elif dropdown.value == "Percent Identity":
                similarity,matrix = percent_identity(seq_records)
            elif dropdown.value == "Mutual Identity":
                similarity,matrix = Mutual_Identity(seq_records)
            else:
                similarity,matrix = sum_of_pairs(seq_records, int(text_input_gap.value), int(text_input_match.value), int(text_input_mismatch.value))

            text_ouput_seqScore.value = str(similarity)

            mat1 = np.zeros((5,5))
            mat1[0:len(seq_records), 0:len(seq_records)] = matrix  
            mat_df = pd.DataFrame(mat1, index=short_labels, columns=short_labels)
            mat_df = pd.DataFrame(mat_df.stack(),columns=['val']).reset_index()
            cm_source.data.update(mat_df)

        else:
            callback_holder.visible = True
            callback_holder.text = "Your letters should be A G T C"
    else:
        if check_MGM() == True:
            callback_holder.visible = True
            callback_holder.text = "Please write Gap, Match, Mismatch for multiple alignment"
        elif is_empty() > 2:
            callback_holder.visible = True
            callback_holder.text = "Please write yor sequence for multiple alignment"
        else:
            callback_holder.visible = True
            callback_holder.text = "Choose more than 2 sequences for multiple alignment"

# error message
callback_holder = PreText(text='', css_classes=['hidden'], visible=False)

# update function should take attribute , old and new 
# update for the matrix score of MSA based on the chosen method 
def update(attr,old,new):
    multiple_alignment()

dropdown = Select(value="Sum of Pairs", options=["Sum of Pairs","Percent Identity","Mutual Identity"],width=200)
dropdown.on_change('value', update)

GA_button = Button(label="Global Alignment", button_type="primary", width=200)
GA_button.on_click(global_alignment)

LA_button = Button(label="Local Alignment", button_type="primary", width=200)
LA_button.on_click(local_alignment)

MSA_button = Button(label="Multiple Alignment", button_type="primary", width=200)
MSA_button.on_click(multiple_alignment)


# USER INTERACTIONS
# input text to change values based on user input sequences 
text_input_seq1 = TextInput(value="", title="Sequence 1")
text_input_seq2 = TextInput(value="", title="Sequence 2")
text_input_seq3 = TextInput(value="", title="Sequence 3")
text_input_seq4 = TextInput(value="", title="Sequence 4")
text_input_seq5 = TextInput(value="", title="Sequence 5")

text_ouput_seqScore = TextInput(value="", title="Sequence Score")

text_input_gap = TextInput(value="", title="Sequence Gap")
text_input_match = TextInput(value="", title="Sequence Match")
text_input_mismatch = TextInput(value="", title="Sequence Mismatch")


LABELS = ["Sequence 1", "Sequence 2", "Sequence 3", "Sequence 4", "Sequence 5"]
checkbox_button_group = CheckboxButtonGroup(labels=LABELS)

title = Div(text='<h1 style="text-align: center;align-items: center;align-self: center;align-content: center;place-items: center;position: center;"> Sequence Alignment </h1>')

input_MGM = Row(text_input_gap,text_input_match,text_input_mismatch)
input_seq1 = Row(text_input_seq1,text_input_seq2)
input_seq2 = Row(text_input_seq3,text_input_seq4)

row = Row(GA_button, LA_button, MSA_button,dropdown,file_input,text_ouput_seqScore)
app = Column(title,input_MGM,row,checkbox_button_group,callback_holder,input_seq1,input_seq2,text_input_seq5,p,met)
curdoc().add_root(app)

