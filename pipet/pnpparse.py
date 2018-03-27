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
#!  @file         pnpparse.py
#!  @brief        Raw data parsing package
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

from __future__ import print_function
import numpy
import sys

raw_type = numpy.uint16
raw_type_size = 2
single_type_size = 10
channels = ['xa','xb','ya','yb']
extended_channels = channels + ['sum']
single_type = numpy.dtype([
    ('mrk',numpy.uint16), ('dip',numpy.uint16), ('dco',numpy.uint16),
    ('tb0',numpy.uint16), ('tb1',numpy.uint16), ('tb2',numpy.uint16), ('tb3',numpy.uint16),
    ('ck0',numpy.uint16), ('ck1',numpy.uint16), ('ck2',numpy.uint16), ('ck3',numpy.uint16), ('ck4',numpy.uint16),
    ('xa',numpy.uint16), ('xb',numpy.uint16), ('ya',numpy.uint16), ('yb',numpy.uint16),
    ('sum',numpy.uint32)
    ])
event_type = numpy.dtype([('a',single_type),('b',single_type),('pos',numpy.uint32),('status',numpy.uint16)])
event_type_size = 2 * single_type_size
chk_vector = numpy.array([[4,0,1,2,3]],dtype=numpy.uint16)
convmap_field_type = numpy.dtype([('name','a3'),('word',numpy.uint8),('len',numpy.uint8),('ofs',numpy.uint8)])
convmap = numpy.array([('mrk',0,6,0), ('dip',0,6,6), ('dco',0,1,12), ('tb0',4,1,12), ('tb1',3,1,12), ('tb2',2,1,12), ('tb3',1,1,12),
                       ('ck0',0,3,13), ('ck1',1,3,13), ('ck2',2,3,13), ('ck3',3,3,13), ('ck4',4,3,13),
                       ('xa',1,12,0), ('xb',2,12,0), ('ya',3,12,0), ('yb',4,12,0),
                      ],dtype=convmap_field_type)
dummy_dip = 15
dummy_chk = 5

def and_reduce(*l):
    if len(l) == 1:
        return l[0].astype(bool)
    elif len(l) == 2:
        return numpy.logical_and(l[0],l[1])
    elif len(l) > 2:
        return and_reduce(and_reduce(*l[:2]),and_reduce(*l[2:]))

def or_reduce(*l):
    if len(l) == 1:
        return l[0].astype(bool)
    elif len(l) == 2:
        return numpy.logical_or(l[0],l[1])
    elif len(l) > 2:
        return or_reduce(or_reduce(*l[:2]),or_reduce(*l[2:]))

def header():
    s  =        '-' * (11+9+15+len(extended_channels)*2*6) + '\n'
    s += '           FLAGS    DD  MRK   DAQ  ' + ''.join([''.join([(j+i).upper().rjust(6) for i in extended_channels]) for j in ['a','b']])
    s += '\n' + '-' * (11+9+15+len(extended_channels)*2*6)
    return s

def evt2str(e):
    flags_a = ''.join([['·','●'][e['a']['tb%d'%i]] for i in range(4)])
    flags_b = ''.join([['·','●'][e['b']['tb%d'%i]] for i in range(4)])
    type_a = ' ' if e['a']['dco'] == 0 else '*'
    type_b = ' ' if e['b']['dco'] == 0 else '*'
    s  = ('%s%s '%(flags_a,flags_b))
    s += ('%s%s '%(type_a,type_b))
    s += str(e['a']['mrk']).rjust(2) + '/' + str(e['b']['mrk']).ljust(3)
    s += str(e['a']['dip']).rjust(2) + '/' + str(e['b']['dip']).ljust(3)
    s += ''.join([''.join([str(e[j][i]).rjust(6) for i in extended_channels]) for j in ['a','b']])
    return s

def raw2evt(raw):
    sng = raw2sng(raw)
    evt = numpy.zeros(sng.size/2,dtype=event_type)
    evt['a'] = sng[::2]
    evt['b'] = sng[1::2]
    marker_match_check = and_reduce(
        evt['a']['mrk'] != evt['b']['mrk'],
        evt['a']['dip'] != dummy_dip,
        evt['b']['dip'] != dummy_dip
        )
    if marker_match_check.any():
        t = numpy.nonzero(marker_match_check)[0]
        print ('Warning! Marker match check failed (%d times, first occurrences at event %s).'%(t.size,', '.join(map(str,t[:5]+1))),file=sys.stderr)

    evt_mrk = numpy.choose(evt['a']['dip'] == dummy_dip, [evt['a']['mrk'], evt['b']['mrk']])
    marker_continuity_check = numpy.concatenate([
            [0],and_reduce(
                (evt_mrk[1:] - evt_mrk[:-1]) != 1,
                (evt_mrk[:-1] - evt_mrk[1:]) != 63,
                numpy.logical_not(and_reduce(evt['a'][1:]['tb0'],evt['a'][:-1]['tb0'],evt_mrk[1:]-evt_mrk[:-1]==0)),
                numpy.logical_not(and_reduce(evt['b'][1:]['tb0'],evt['b'][:-1]['tb0'],evt_mrk[1:]-evt_mrk[:-1]==0)),
                )
        ])
    if marker_continuity_check.any():
        t = numpy.nonzero(marker_continuity_check)[0]
        print ('Warning! Marker continuity check failed (%d times, first occurrences at event %s).'%(t.size,', '.join(map(str,t[:5]+1))),file=sys.stderr)
    return evt

def raw2sng(raw):
    step = single_type_size/raw_type_size
    sng = numpy.zeros(raw.size/step,dtype=single_type)
    for i in convmap:
        sng[i['name']] = numpy.bitwise_and(numpy.right_shift(raw[i['word']::step],i['ofs']),int('1'*i['len'],2))
    ck_array = numpy.array([sng['ck%d'%i] for i in range(chk_vector.size)]).transpose().ravel()
    signature_check = and_reduce(
        ck_array != numpy.repeat(chk_vector.transpose(),sng.size,axis=1).transpose().ravel(),
        ck_array != dummy_chk)
    if signature_check.any():
        t = numpy.nonzero(signature_check)[0]
        print ('Warning! Signature check failed (%d times, first occurrences at word %s).'%(t.size,', '.join(map(str,t[:5]+1))),file=sys.stderr)
    for i in channels:
        sng['sum'] += sng[i]
    return sng

def evt2raw(evt):
    sng = numpy.zeros(evt.size*2,dtype=single_type)
    sng[::2]  = evt['a']
    sng[1::2] = evt['b']
    raw = sng2raw(sng)
    return raw

def sng2raw(sng):
    step = single_type_size/raw_type_size
    raw = numpy.zeros(sng.size*step,dtype=raw_type)
    raw = numpy.bitwise_or(raw,numpy.left_shift(numpy.bitwise_and(numpy.repeat(chk_vector,sng.size,axis=0).ravel(),0x7),13))
    for i in convmap:
        raw[i['word']::step] = numpy.bitwise_or(raw[i['word']::step],numpy.left_shift(numpy.bitwise_and(sng[i['name']],int('1'*i['len'],2)),i['ofs']))
    return raw

def crop(evt,beg,end,step=1):
    if end != None:
        evt = evt[:end*step]
    if beg != None:
        evt = evt[beg*step:]
    return evt

def filter_events(evt,args):
    if args.d == 0:
        evt = evt[and_reduce(evt['a']['dco']==0,evt['b']['dco']==0)]
    elif args.d == 1:
        evt = evt[or_reduce(evt['a']['dco']==1,evt['b']['dco']==1)]
    #print args.fl
    if args.fl == 0:
        evt = evt[
            and_reduce(
                evt['a']['tb0']==0,evt['a']['tb1']==0,
                evt['a']['tb2']==0,evt['a']['tb3']==0,
                evt['b']['tb0']==0,evt['b']['tb1']==0,
                evt['b']['tb2']==0,evt['b']['tb3']==0
                )
        ]
    elif args.fl == 1:
        evt = evt[
            or_reduce(
                evt['a']['tb0']==1,evt['a']['tb1']==1,
                evt['a']['tb2']==1,evt['a']['tb3']==1,
                evt['b']['tb0']==1,evt['b']['tb1']==1,
                evt['b']['tb2']==1,evt['b']['tb3']==1
                )
        ]
    filter_indices = numpy.ones(evt.size,dtype=numpy.bool)
    for i in extended_channels:
        if vars(args)['l'+i] != None:
            filter_indices = and_reduce(filter_indices,evt['a'][i] >= vars(args)['l'+i])
            filter_indices = and_reduce(filter_indices,evt['b'][i] >= vars(args)['l'+i])
        if vars(args)['u'+i] != None:
            filter_indices = and_reduce(filter_indices,evt['a'][i] <= vars(args)['u'+i])
            filter_indices = and_reduce(filter_indices,evt['b'][i] <= vars(args)['u'+i])
    evt = evt[filter_indices]
    return evt

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Parse a plugnpet acquisition file',
    epilog='Example: ./pnpparse.py -i pet_20120225_190516.evt.dat -o pet_20120225_190516_lt4000.evt.dat -usum 4000')
    parser.add_argument("-i", dest="ipath", help="Input binary file", metavar="binary filename", action="store", required=True)
    parser.add_argument("-o", dest="opath", help="Output binary file", metavar="binary filename", action="store")
    parser.add_argument("-p", help="Print events", action="store_true")
    parser.add_argument("-d", help="Filter delayed [0=No delayed, 1=Only delayed, 2=Both]", metavar="option", action="store",type=int,default=2)
    parser.add_argument("-fl", help="Filter flagged [0=Not flagged, 1=Only flagged, 2=Both]", metavar="option", action="store",type=int,default=2)
    parser.add_argument("-fe", help="Maximum bytes to read from file", metavar="bytes", action="store",type=int,default=-1)
    parser.add_argument("-bb", help="Start from (before filtering)", metavar="event number", action="store",type=int)
    parser.add_argument("-be", help="End to (before filtering)", metavar="event number", action="store",type=int)
    parser.add_argument("-ab", help="Start from (after filtering)", metavar="event number", action="store",type=int)
    parser.add_argument("-ae", help="End to (after filtering)", metavar="event number", action="store",type=int)
    for i in extended_channels:
        parser.add_argument("-l"+i, help="Lower threshold for ADC "+i, metavar="threshold", action="store",type=int)
    for i in extended_channels:
        parser.add_argument("-u"+i, help="Upper threshold for ADC "+i, metavar="threshold", action="store",type=int)
    args = parser.parse_args()
    if args.fe == -1:
        max_bytes = -1
    else:
        max_bytes = args.fe / raw_type_size
    raw = numpy.fromfile(args.ipath,dtype=raw_type,count=args.fe)
    raw = crop(raw,args.bb,args.be,event_type_size/raw_type_size)
    evt = raw2evt(raw)
    filter_events(evt,args)
    evt = crop(evt,args.ab,args.ae)
    if args.opath != None:
        raw = evt2raw(evt)
        raw.tofile(args.opath)
    if args.p:
        print (header())
        for i,e in enumerate(evt):
            print (('%d:'%(i+1)).rjust(10),evt2str(e))
