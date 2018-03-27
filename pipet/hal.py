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
#!  @file         hal.py
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
import os
import sys
import time
import datetime
import json
import ok
import re
import numpy
import itertools
import inspect
import traceback
from scipy.interpolate import UnivariateSpline
from fractions import gcd
from .utility import *

C_OK_PIPE_ERRORS = ['InvalidEndpoint','InvalidBlockSize','Failed', 'Timeout']
C_BTPIPE_READY_DETPH = 1024
C_DAQ_EVENT_BYTES = 20
C_READ_BUF_MAX_SIZE = 128*1024*1024 # 128 MB
C_READ_BUF_MAX_EVENTS = ((C_READ_BUF_MAX_SIZE//C_DAQ_EVENT_BYTES)//C_BTPIPE_READY_DETPH)*C_BTPIPE_READY_DETPH
C_PEDESTAL_AUTO_EVENTS = 4*1024*1024
C_WORD_SIZE = 2
C_BUS_MAP = {
    'x' : {
        'read'  : [(0x30,i) for i in range(16)]
                + [(0x31,i) for i in range(16)]
                + [(0x32,i) for i in range(16)]
                + [(0x33,i) for i in range(16)],
        'write' : [(0x10,i) for i in range(16)]
                + [(0x11,i) for i in range(16)]
                + [(0x12,i) for i in range(16)]
                + [(0x13,i) for i in range(16)],
        'oe'    : [(0x18,i) for i in range(16)]
                + [(0x19,i) for i in range(16)]
                + [(0x1A,i) for i in range(16)]
                + [(0x1B,i) for i in range(16)]
        },
    'y' : {
        'read'  : [(0x34,i) for i in range(16)]
                + [(0x35,i) for i in range(16)]
                + [(0x36,i) for i in range(16)]
                + [(0x37,i) for i in range(16)],
        'write' : [(0x14,i) for i in range(16)]
                + [(0x15,i) for i in range(16)]
                + [(0x16,i) for i in range(16)]
                + [(0x17,i) for i in range(16)],
        'oe'    : [(0x1C,i) for i in range(16)]
                + [(0x1D,i) for i in range(16)]
                + [(0x1E,i) for i in range(16)]
                + [(0x1F,i) for i in range(16)]
    }
}

C_DAQ_MAP = {
    '1': {
        # signal : (bus, bit, mode)
        'nclr' : ('x',0,'w'),
        'dco'  : ('x',1,'w'),
        'trg'  : ('x',2,'w'),
        'dav'  : ('x',3,'r'),
        'req'  : ('x',4,'w'),
        'ack'  : ('x',5,'r'),
        'aux'  : ('x',6,'w'),
        'rej'  : ('x',7,'w'),
        'bsy'  : ('x',8,'r'),
        'tbd0' : ('x',10,'w'),
        'tbd1' : ('x',11,'w'),
        'tbd2' : ('x',12,'w'),
        'tbd3' : ('x',13,'w'),
        'rd0'  : ('x',30,'r'),
        'rd1'  : ('x',31,'r'),
        'rd2'  : ('x',32,'r'),
        'rd3'  : ('x',33,'r'),
        'rd4'  : ('x',34,'r'),
        'rd5'  : ('x',35,'r'),
        'rd6'  : ('x',36,'r'),
        'rd7'  : ('x',37,'r'),
        'rd8'  : ('x',38,'r'),
        'rd9'  : ('x',39,'r'),
        'rd10' : ('x',40,'r'),
        'rd11' : ('x',41,'r'),
        'rd12' : ('x',42,'r'),
        'rd13' : ('x',43,'r'),
        'rd14' : ('x',44,'r'),
        'rd15' : ('x',45,'r'),
        'mrk0' : ('x',50,'w'),
        'mrk1' : ('x',51,'w'),
        'mrk2' : ('x',52,'w'),
        'mrk3' : ('x',53,'w'),
        'mrk4' : ('x',54,'w'),
        'mrk5' : ('x',55,'w'),
    },
    '2': {
        'nclr' : ('y',0,'w'),
        'dco'  : ('y',1,'w'),
        'trg'  : ('y',2,'w'),
        'dav'  : ('y',3,'r'),
        'req'  : ('y',4,'w'),
        'ack'  : ('y',5,'r'),
        'aux'  : ('y',6,'w'),
        'rej'  : ('y',7,'w'),
        'bsy'  : ('y',8,'r'),
        'tbd0' : ('y',10,'w'),
        'tbd1' : ('y',11,'w'),
        'tbd2' : ('y',12,'w'),
        'tbd3' : ('y',13,'w'),
        'rd0'  : ('y',30,'r'),
        'rd1'  : ('y',31,'r'),
        'rd2'  : ('y',32,'r'),
        'rd3'  : ('y',33,'r'),
        'rd4'  : ('y',34,'r'),
        'rd5'  : ('y',35,'r'),
        'rd6'  : ('y',36,'r'),
        'rd7'  : ('y',37,'r'),
        'rd8'  : ('y',38,'r'),
        'rd9'  : ('y',39,'r'),
        'rd10' : ('y',40,'r'),
        'rd11' : ('y',41,'r'),
        'rd12' : ('y',42,'r'),
        'rd13' : ('y',43,'r'),
        'rd14' : ('y',44,'r'),
        'rd15' : ('y',45,'r'),
        'mrk0' : ('y',50,'w'),
        'mrk1' : ('y',51,'w'),
        'mrk2' : ('y',52,'w'),
        'mrk3' : ('y',53,'w'),
        'mrk4' : ('y',54,'w'),
        'mrk5' : ('y',55,'w'),
    },
}

C_ACQUISITION_MODE_MAP = {
        'auto'    : {'ep': 0x00, 'value': 0x0<<3, 'mask': 0x3<<3},
        'single_a': {'ep': 0x00, 'value': 0x1<<3, 'mask': 0x3<<3},
        'single_b': {'ep': 0x00, 'value': 0x2<<3, 'mask': 0x3<<3},
        'coinc'   : {'ep': 0x00, 'value': 0x3<<3, 'mask': 0x3<<3}
}

def make_frame_list(tot_events,frame_events):
    tot_events, frame_events = int(tot_events), int(frame_events)
    return [frame_events]*(tot_events//frame_events)+[tot_events%frame_events]

C_ACQUISITION_MODE_FRAME_SCHEME = {
        'auto'    : lambda events, rates: make_frame_list(events,C_PEDESTAL_AUTO_EVENTS),
        'single_a': lambda events, rates: make_frame_list(events,max(lcm((C_BTPIPE_READY_DETPH,C_DAQ_EVENT_BYTES)),(rates['cfd_a']//C_BTPIPE_READY_DETPH)*C_BTPIPE_READY_DETPH)),
        'single_b': lambda events, rates: make_frame_list(events,max(lcm((C_BTPIPE_READY_DETPH,C_DAQ_EVENT_BYTES)),(rates['cfd_b']//C_BTPIPE_READY_DETPH)*C_BTPIPE_READY_DETPH)),
        'coinc'   : lambda events, rates: make_frame_list(events,max(lcm((C_BTPIPE_READY_DETPH,C_DAQ_EVENT_BYTES)),(rates['cnc_a']//C_BTPIPE_READY_DETPH)*C_BTPIPE_READY_DETPH)),
}

C_CONFIGURATION_MAP = {
    'oscillator1_on' : lambda x: {'ep': 0x00, 'value': int(bool(x))<<0, 'mask': 1<<0},
    'oscillator2_on' : lambda x: {'ep': 0x00, 'value': int(bool(x))<<1, 'mask': 1<<1},
    'delay_a' : lambda x: {'ep': 0x0E, 'value': int(x), 'mask': 0xffff},
    'delay_b' : lambda x: {'ep': 0x0F, 'value': int(x), 'mask': 0xffff},
    'acquisition_on' : lambda x: {'ep': 0x00, 'value': int(bool(x))<<2, 'mask': 1<<2},
    'acquisition_mode'  : lambda x: C_ACQUISITION_MODE_MAP[x],
    'block_size'     : lambda x: {'ep': 0x04, 'value': int(x)         , 'mask': 0x07ff},
}

class Bus():
    def __init__(self,fpga,name):
        self.fpga = fpga
        self.name = name
        self.verbose = True
    def read_update(self):
        self.fpga.UpdateWireOuts()
    def write_update(self):
        self.fpga.UpdateWireIns()
    def oe(self,i,val,update=True):
        if self.verbose:
            print ('Bus',self.name,'oe',i,val,end=' ')
        ep, bit = C_BUS_MAP[self.name]['oe'][i]
        if self.verbose:
            print ('-> ep',hex(ep),'bit',bit)
        self.fpga.SetWireIn(ep,int(bool(val))<<bit,mask=1<<bit,update=update)
    def write(self,i,val,update=True):
        if self.verbose:
            print ('Bus',self.name,'write',i,val,end=' ')
        ep, bit = C_BUS_MAP[self.name]['write'][i]
        if self.verbose:
            print ('-> ep',hex(ep),'bit',bit)
        self.fpga.SetWireIn(ep,int(bool(val))<<bit,mask=1<<bit,update=update)
    def read(self,i,update=True):
        if self.verbose:
            print ('Bus',self.name,'read',i,end=' ')
        ep, bit = C_BUS_MAP[self.name]['read'][i]
        if self.verbose:
            print ('-> ep',hex(ep),'bit',bit)
        return int((self.fpga.GetWireOut(ep,update=update)>>bit) & 1)

class Daq():
    def __init__(self,bus_array,daq_id):
        self.read_timeout_s = 1
        self.read_poll_interval_s = 0.01
        self.daq_id = daq_id
        self.bus_array = bus_array
        wide_signals = {}
        is_wide_re = re.compile(r'(?P<name>[a-z]+)(?P<number>\d+)')
        for i in C_DAQ_MAP[self.daq_id]:
            bus_n, bit, mode = C_DAQ_MAP[self.daq_id][i]
            setattr(self,i,self.make_signal_handler(bus_n,bit,mode))
            r = is_wide_re.match(i)
            if r:
                rgd = r.groupdict()
                name = rgd['name']
                number = int(rgd['number'])
                if name in wide_signals:
                    if wide_signals[name]['mode'] != mode:
                        raise RuntimeError('Mode mismatch in wide_signal array: '+name)
                    wide_signals[name]['array'][number] = (bus_n,bit)
                else:
                    wide_signals[name] = {
                        'array': {number : (bus_n,bit)},
                        'mode' : mode
                    }
        for i in wide_signals:
            setattr(self,i,self.make_wide_signal_handler(wide_signals[i]))

    def oe_all(self):
        for i in C_DAQ_MAP[self.daq_id]:
            bus_n, bit, mode = C_DAQ_MAP[self.daq_id][i]
            if mode == 'w':
                getattr(self,i)(oe=1)

    def oe_none(self):
        for i in C_DAQ_MAP[self.daq_id]:
            bus_n, bit, mode = C_DAQ_MAP[self.daq_id][i]
            if mode == 'w':
                getattr(self,i)(oe=0)

    def print_all(self):
        print ('NCLR:',self.nclr())
        print ('DCO :',self.dco())
        print ('TRG :',self.trg())
        print ('DAV :',self.dav())
        print ('REQ :',self.req())
        print ('ACK :',self.ack())
        print ('AUX :',self.aux())
        print ('BSY :',self.bsy())
        print ('TBD :',self.tbd())
        print ('RD  :',self.rd())
        print ('MRK :',self.mrk())

    def trig(self):
        self.trg(o=0,oe=1,read=False)
        self.trg(o=1,read=False)
        self.trg(o=0,read=False)

    def read_word(self):
        if not self.dav():
            print ('No data available.')
            return
        self.req(o=0,oe=1,read=False)
        if self.ack():
            raise RuntimeError('ack is high while req is low.')
        self.req(o=1,read=False)
        start = time.time()
        while True:
            if time.time() - start > self.read_timeout_s:
                raise RuntimeError('Timeout: Daq didn\'t ack after req.')
            if self.ack():
                break
            time.sleep(self.read_poll_interval_s)
        rd = self.rd()
        self.req(o=0,read=False)
        return rd

    def print_event_format(self):
        print ('| WORD | D15 | D14 | D13 |  D12 | D11 | D10 | D09 | D08 | D07 | D06 | D05 | D04 | D03 | D02 | D01 | D00 |')
        print ('---------------------------------------------------------------------------------------------------------')
        print ('|  0   |  1  |  0  |  0  |  DCO |              DAQ ID               |                MRK                |')
        print ('|  1   |  0  |  0  |  0  | TBD3 |                                  XA                                   |')
        print ('|  2   |  0  |  0  |  1  | TBD2 |                                  XB                                   |')
        print ('|  3   |  0  |  1  |  0  | TBD1 |                                  YA                                   |')
        print ('|  4   |  0  |  1  |  1  | TBD0 |                                  YB                                   |')

    def read_event(self,print_event=True):
        event = [self.read_word() for i in range(5)]
        if print_event:
            print (format_str_event_1(event))
        return event

    def make_signal_handler(self,bus_n,bit,mode):
        if mode == 'w':
            def signal_setter(o=None,oe=None,read=True):
                if o != None:
                    self.bus_array[bus_n].write(bit,o)
                if oe != None:
                    self.bus_array[bus_n].oe(bit,oe)
                if read:
                    return self.bus_array[bus_n].read(bit)
            return signal_setter
        elif mode == 'r':
            def signal_getter():
                return self.bus_array[bus_n].read(bit)
            return signal_getter

    def make_wide_signal_handler(self,signals_dict):
        if signals_dict['mode'] == 'w':
            signal_width = max([n for n in signals_dict['array']])+1
            def wide_signal_setter(o=None,oe=None,read=True):
                if o != None:
                    bitarray = map(int,(bin(o)[2:].rjust(len(signals_dict['array']),'0'))[::-1])
                    for n in signals_dict['array']:
                        bus_n, bit = signals_dict['array'][n]
                        self.bus_array[bus_n].write(bit,bitarray[n],update=False)
                        self.bus_array[bus_n].write_update()
                if oe != None:
                    for n in signals_dict['array']:
                        bus_n, bit = signals_dict['array'][n]
                        self.bus_array[bus_n].oe(bit,oe,update=False)
                        self.bus_array[bus_n].write_update()
                if read:
                    bitarray = ['0']*signal_width
                    for n in signals_dict['array']:
                        bus_n, bit = signals_dict['array'][n]
                        bitarray[n] = self.bus_array[bus_n].read(bit)
                    return ''.join(map(str,bitarray[::-1]))
            return wide_signal_setter
        elif signals_dict['mode'] == 'r':
            signal_width = max([n for n in signals_dict['array']])+1
            def wide_signal_getter():
                bitarray = ['0']*signal_width
                for n in signals_dict['array']:
                    bus_n, bit = signals_dict['array'][n]
                    bitarray[n] = self.bus_array[bus_n].read(bit)
                return ''.join(map(str,bitarray[::-1]))
            return wide_signal_getter

class okDevice():
    def __init__(self,verbose=False):
        self.verbose = verbose

    def fsm_decode(self,name,code):
        if self.fsm_map == None:
            return 'NODEFS'
        elif name not in self.fsm_map:
            return 'UNREG'
        else:
            for i in self.fsm_map[name]:
                if self.fsm_map[name][i] == code:
                    return i
            return 'NOTFOUND'

    def UpdateWireOuts(self):
        if self.verbose:
            print('UpdateWireOuts()')
        self.xem.UpdateWireOuts()

    def UpdateWireIns(self):
        if self.verbose:
            print('UpdateWireIns()')
        self.xem.UpdateWireIns()

    def GetWireOut(self, ep, update=True):
        if update:
            self.UpdateWireOuts()
        ret = self.xem.GetWireOutValue(ep)
        if self.verbose:
            print('GetWireOutValue({}) -- = "{}" ({})'.format(hex(ep),bin(ret)[2:].rjust(16,'0'),ret))
        return ret

    def SetWireIn(self, ep, value, mask=0xffff, update=True):
        if self.verbose:
            print('SetWireInValue({},"{}","{}")'.format(hex(ep),bin(value)[2:].rjust(16,'0'),bin(mask)[2:].rjust(16,'0')))
        self.xem.SetWireInValue(ep,value,mask)
        if update:
            self.UpdateWireIns()

    def Trigger(self, ep, n):
        if self.verbose:
            print('ActivateTriggerIn({},{})'.format(hex(ep),n))
        self.xem.ActivateTriggerIn(ep, n)

    def IsTriggered(self, ep, n):
        if self.verbose:
            print('UpdateTriggerOuts()')
        self.xem.UpdateTriggerOuts()
        ret = self.xem.IsTriggered(ep, n)
        if self.verbose:
            print('IsTriggered({},{}) -- = {}'.format(hex(ep),n,str(ret)))
        return ret

    def WaitForTrigger(self, ep, n, timeout = None, polling_time = 1e-3):
        start = time.time()
        while True:
            self.IsTriggered(ep, n)
            time.sleep(polling_time)
            if timeout and time.time() - start > timeout:
                break

    def WriteToPipeIn(self, ep, buf):
        if self.verbose:
            print('WriteToPipeIn({},buf)'.format(hex(ep)))
        return self.xem.WriteToPipeIn(ep, buf)

    def ReadFromPipeOut(self, ep, buf, bsize = None):
        if bsize == None:
            if self.verbose:
                print('ReadFromPipeOut({},buf)'.format(hex(ep)))
            return self.xem.ReadFromPipeOut(ep, buf)
        else:
            if self.verbose:
                print('ReadFromBlockPipeOut({},{},buf)'.format(hex(ep),bsize))
            return self.xem.ReadFromBlockPipeOut(ep, bsize, buf)

    def InitializeDevice(self, bitfile):
        self.xem = ok.okCFrontPanel()
        if (self.xem.NoError != self.xem.OpenBySerial('')):
            raise RuntimeError("A device could not be opened. Is one connected?")

        self.devInfo = ok.okTDeviceInfo()
        if (self.xem.NoError != self.xem.GetDeviceInfo(self.devInfo)):
            raise RuntimeError("Unable to retrieve device information.")

        print("-"*60)
        print("Pipet learning system - University of Pisa")
        print("-"*60)
        print(datetime.datetime.now().strftime('Board initialized: %H:%M:%S %d-%m-%Y'))
        print("    Firmware path: %s (%s)"%(bitfile,time.strftime("%d/%m/%Y %H:%M",time.localtime(os.path.getmtime(bitfile)))))
        print("-"*60)

        self.xem.LoadDefaultPLLConfiguration()

        # Download the configuration file.
        if (self.xem.NoError != self.xem.ConfigureFPGA(bitfile)):
            raise RuntimeError("FPGA configuration failed.")

        # Check for FrontPanel support in the FPGA configuration.
        if (False == self.xem.IsFrontPanelEnabled()):
            raise RuntimeError("FrontPanel support is not available.")

        print ("Device ready.") #print ("FrontPanel support is available.")

        # Used for debug:
        firmware_path = os.path.split(bitfile)[0]
        fsm_map_file = os.path.join(firmware_path,'fsm-map.json')
        self.fsm_map = None
        if os.path.isfile(fsm_map_file):
            try:
                self.fsm_map = json.load(open(fsm_map_file))
            except:
                print('Warning! Bad fsm map definitions!',file=sys.stderr)
                print(traceback.format_exc(),file=sys.stderr)

class pipet():
    def __init__(self):
        '''Main instance constructor'''
        self.fpga = okDevice(verbose=False)
        self.bus_array = {'x': Bus(self.fpga,'x'), 'y': Bus(self.fpga,'y')}
        self.daq1 = Daq(self.bus_array,'1')
        self.daq2 = Daq(self.bus_array,'2')
    def set_bus_verbosity(self,verbosity):
        '''Changes system verbosity (for debug only)'''
        for b in self.bus_array:
            self.bus_array[b].verbose = verbosity
    def init(self, bitfile = '../firmware/default.bit'):
        '''Initializes the device and uploads the firmware'''
        self.fpga.InitializeDevice(bitfile)
        init_time = time.time()
        self.set_bus_verbosity(False)
        self.daq1.oe_all()
        self.daq2.oe_all()
        self.reset_daqs()
        self.oe_init()
        time.sleep(2-(time.time()-init_time)) # Needs to wait 2 seconds after the initialization to make the counters steady
    def reset_daqs(self):
        '''Reset DAQ boards'''
        self.daq1.nclr(0)
        self.daq2.nclr(0)
        self.daq1.nclr(1)
        self.daq2.nclr(1)
        time.sleep(0.5)
    def config(self,name,value,update=True):
        '''Configures FPGA registers (not intended for the user)'''
        self.fpga.SetWireIn(update=update,**(C_CONFIGURATION_MAP[name](value)))
    def delay(self,steps=0):
        '''Sets the delay steps between DCFD_A and DCFD_B'''
        C_MAX_CFD_DELAY_STEPS = 256
        assert(abs(steps)<C_MAX_CFD_DELAY_STEPS)
        if steps > 0:
            self.config('delay_a',steps,update=False)
            self.config('delay_b',0,update=True)
        elif steps < 0:
            self.config('delay_a',0,update=False)
            self.config('delay_b',-steps,update=True)
        else:
            self.config('delay_a',0,update=False)
            self.config('delay_b',0,update=True)
        time.sleep(2)
    def oe_init(self):
        '''Enable output buses (not intended for the user)'''
        self.bus_array['x'].oe(16,1,update=False)
        self.bus_array['y'].oe(16,1,update=False)
        self.bus_array['x'].oe(18,1,update=False)
        self.bus_array['y'].oe(18,1,update=False)
        self.bus_array['x'].oe(20,1,update=False)
        self.bus_array['x'].oe(20,1,update=False)
        self.bus_array['x'].oe(21,1,update=False)
        self.bus_array['y'].oe(21,1,update=True)
    def osc(self,enable):
        '''Enables the internal oscillator to the DCFD output for CW measurements'''
        self.bus_array['x'].write(16,0,update=False)
        self.bus_array['y'].write(16,0,update=False)
        self.bus_array['x'].write(18,0,update=False)
        self.bus_array['y'].write(18,0,update=True)
        self.config('oscillator1_on',enable,update=False)
        self.config('oscillator2_on',enable,update=True)
    def read_acq_pipe(self,length):
        '''Read data from internal FIFO (not intended for the user)'''
        assert(length >= C_DAQ_EVENT_BYTES)
        assert(length % C_WORD_SIZE == 0)
        assert(length <= C_READ_BUF_MAX_SIZE)
        block_size = int((numpy.nonzero(numpy.logical_not(length%numpy.arange(2,1025,2)))[0][-1]+1)*2)
        buf = bytearray(length)
        #print ('reading {} B frames with {} B blocks'.format(len(buf),block_size))
        self.config('block_size',block_size,update=True)
        ret = self.fpga.ReadFromPipeOut(0xA1, buf, bsize = block_size)
        if ret < 0:
            print ('Error:',[i for i in C_OK_PIPE_ERRORS if ret == getattr(self.fpga.xem,i)])
            return None
        return numpy.frombuffer(buf,dtype=numpy.uint16)
    def acquire(self,mode='auto',events=1000,frames=None,show=False):
        '''Acquire data either in auto, single_a, single_b or coinc mode'''
        if mode not in C_ACQUISITION_MODE_MAP:
            raise RuntimeError('Unknown acquisition mode')
        self.reset_daqs()
        self.config('acquisition_on',0,update=True)
        self.config('acquisition_on',1,update=False)
        self.config('acquisition_mode',mode,update=True)
        ret = []
        if frames == None:
            time.sleep(1)
            r = self.rates()
            frames = C_ACQUISITION_MODE_FRAME_SCHEME[mode](events,r)
            for cur_events in log_progress(frames,show=show):
                ret.append(self.read_acq_pipe(cur_events*C_DAQ_EVENT_BYTES))
        else:
            for i in log_progress(range(frames),show=show):
                ret.append(self.read_acq_pipe(events*C_DAQ_EVENT_BYTES))
        self.config('acquisition_on',0,update=True)
        return numpy.concatenate(ret)
    def rates(self,print_rates=False):
        '''Return trigger rates'''
        cfd_a_lsb = self.fpga.GetWireOut(0x24,update=True)
        cfd_a_msb = self.fpga.GetWireOut(0x25,update=False)
        cfd_b_lsb = self.fpga.GetWireOut(0x26,update=False)
        cfd_b_msb = self.fpga.GetWireOut(0x27,update=False)
        cnc_a_lsb = self.fpga.GetWireOut(0x28,update=False)
        cnc_a_msb = self.fpga.GetWireOut(0x29,update=False)
        cnc_b_lsb = self.fpga.GetWireOut(0x2A,update=False)
        cnc_b_msb = self.fpga.GetWireOut(0x2B,update=False)
        dly_a_lsb = self.fpga.GetWireOut(0x2C,update=False)
        dly_a_msb = self.fpga.GetWireOut(0x2D,update=False)
        dly_b_lsb = self.fpga.GetWireOut(0x2E,update=False)
        dly_b_msb = self.fpga.GetWireOut(0x2F,update=False)
        cfd_a = cfd_a_lsb | (cfd_a_msb << 16)
        cfd_b = cfd_b_lsb | (cfd_b_msb << 16)
        cnc_a = cnc_a_lsb | (cnc_a_msb << 16)
        cnc_b = cnc_b_lsb | (cnc_b_msb << 16)
        dly_a = dly_a_lsb | (dly_a_msb << 16)
        dly_b = dly_b_lsb | (dly_b_msb << 16)
        r = {
                'cfd_a': cfd_a,
                'cfd_b': cfd_b,
                'cnc_a': cnc_a,
                'cnc_b': cnc_b,
                'dly_a': dly_a,
                'dly_b': dly_b
            }
        if print_rates:
            print ('cfd_a: ',r['cfd_a'],'Hz')
            print ('cfd_b: ',r['cfd_b'],'Hz')
            print ('cnc_a: ',r['cnc_a'],'Hz')
            print ('cnc_b: ',r['cnc_b'],'Hz')
            print ('dly_a: ',r['dly_a'],'Hz')
            print ('dly_b: ',r['dly_b'],'Hz')
        return r

pet = pipet()
dq1 = pet.daq1
dq2 = pet.daq2