#!/usr/bin/env python
# -*- coding: utf-8 -*-
#***************************************************************************
#*                       ______   ____    __°   ______
#*                      / ____/  /  _/   /_/   / ____/
#*                     / /_      / /    /_/   / / __
#*                    / __/    _/ /   _/_/   / /_/ /
#*                   /_/      /___/  /___/   \____/
#*
#*    FUNCTIONAL IMAGING AND INSTRUMENTATION GROUP - UNIVERSITA' DI PISA
#*
#***************************************************************************
#*
#*  Project     : Laboratorio di Fisica Medica
#!  @file         utility.py
#!  @brief        Utility package
#*
#*  Author(s)   : Giancarlo Sportelli (GK)
#*                see AUTHORS for complete info
#*  License     : see LICENSE for info
#*
#***************************************************************************
#*
#*                             R e v i s i o n s
#*
#*--------------------------------------------------------------------------
#*  Timestamp             Author    Version    Description
#*--------------------------------------------------------------------------
#*  22:48 08/02/2016      GK         0.1       Initial design
#*  further revisions are tagged in the git repository
#***************************************************************************

from __future__ import division
from __future__ import print_function
import sys
import numpy
from scipy.interpolate import UnivariateSpline
from fractions import gcd
from ipywidgets import IntProgress, HTML, VBox
from IPython.display import display

def log_progress(sequence, every=None, size=None, show=True):
    is_iterator = False
    if size is None:
        try:
            size = len(sequence)
        except TypeError:
            is_iterator = True
    if size is not None:
        if every is None:
            if size <= 200:
                every = 1
            else:
                every = size / 200 # every 0.5%
    else:
        assert every is not None, 'sequence is iterator, set every'

    if is_iterator:
        progress = IntProgress(min=0, max=1, value=1)
        progress.bar_style = 'info'
    else:
        progress = IntProgress(min=0, max=size, value=0)
    label = HTML()
    box = VBox(children=[label, progress])
    fixed_label = ''
    if show:
        display(box)
        if isinstance(show,str):
            fixed_label = show

    index = 0
    try:
        for index, record in enumerate(sequence, 1):
            if index == 1 or index % every == 0:
                if is_iterator:
                    label.value = '{fixed_label}{index} / ?'.format(fixed_label=fixed_label,index=index)
                else:
                    progress.value = index
                    #label.value = u'{fixed_label}{index} / {size}'.format(
                    label.value = u'{fixed_label}{index}'.format(
                        fixed_label=fixed_label,
                        index='%.1f %%'%(100.*index/size),
                        #index=index,
                        #size=size
                    )
            yield record
    except:
        progress.bar_style = 'danger'
        raise
    else:
        progress.bar_style = 'success'
        progress.value = index
        label.value = '{fixed_label}{result}'.format(fixed_label=fixed_label,result=str('%.1f %%'%(100.*index/size) or '?'))

def lcm(numbers):
    return reduce(lambda x, y: (x*y)/gcd(x,y), numbers, 1)

def fwhm(x,y):
    return numpy.diff(UnivariateSpline(x, y-y.max()/2).roots())[0]

def offset(x,y):
    return numpy.sum(UnivariateSpline(x, y-y.max()/2).roots())/2

def p(stream,events=None):
    if events == None:
        events = int(len(stream)/10)
    for i in range(events):
        print ('-'*30)
        print ('Event',i)
        print ('-'*30)
        print(format_uint16_event_1(stream[i*10:i*10+5]   ))
        print(format_uint16_event_1(stream[i*10+5:i*10+10]))

def print_err(*args,**kwargs):
    kwargs['file']=sys.stderr
    print(*args,**kwargs)

def format_uint16_event_1(event):
    return format_str_event_1([bin(j)[2:].rjust(16,'0') for j in event])

def format_str_event_1(event):
    w = event[0]
    s = ''
    try:
        s += ' '.join([str(i) for i in [w[:3],w[3],w[4:10],w[10:16],' | DAQ: {} ({}), EVT: {} ({})'.format(w[4:10],int(w[4:10],2),w[10:16],int(w[10:16],2))]])
    except:
        s += 'Bad word: ' + str(w)
    s += '\n'
    for i,w in enumerate(event[1:]):
        try:
            s += ' '.join([str(i) for i in [w[:3],w[3],w[4:16],'  | TBD{}: {}, {}: {} ({})'.format(3-i,w[3],['XA','XB','YA','YB'][i],w[4:16],int(w[4:16],2))]])
        except:
            s += 'Bad word: ' + str(w)
        s += '\n'
    s = s.replace('TBD2',' SNG')
    s = s.replace('TBD1',' PAD')
    s = s.replace('TBD0',' DLY')
    return s
