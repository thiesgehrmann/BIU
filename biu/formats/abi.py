#!/usr/bin/env python
#
# abifpy.py
# python module for reading abi trace files
# http://github.com/bow/abifpy
# Modified by Thies Gehrmann
# ABI File reference: http://www6.appliedbiosystems.com/support/software_community/ABIF_File_Format.pdf
# Relevant table: Table 5

"""Python module for reading .ab1 trace files."""

import datetime
import struct
from os.path import splitext, basename

from sys import version_info

from .. import utils

plt         = utils.py.loadExternalModule('matplotlib.pylab')
MaxNLocator = utils.py.loadExternalModule('matplotlib.ticker', 'MaxNLocator')
patches     = utils.py.loadExternalModule('matplotlib.patches') 

RELEASE = False
__version_info__ = ('1', '0', )
__version__ = '.'.join(__version_info__)
__version__ += '-dev' if not RELEASE else ''


__all__ = ['ABI']

# dictionary for deciding which values to extract and contain in self.data
EXTRACT = {
    'APFN2': 'Sequencing Analysis parameters file name',
    'APXV1': 'Analysis Protocol XML schema version',
    'APrN1': 'Analysis Protocol settings name',
    'APrV1': 'Analysis Protocol settings version',
    'APrX1': 'Analysis Protocol XML string',
    'CMNT1': 'Sample Comment',
    'CTID1': 'Container Identifier, a.k.a. plate barcode',
    'CTNM1': 'Container name, usually identical to CTID, but not necessarily so',
    'CTTL1': 'Comment Title',
    'CpEP1': 'Capillary type electrophoresis. 1 for a capillary based machine. 0 for a slab gel based machine.',
    'DATA1': 'Channel 1 raw data',
    'DATA2': 'Channel 2 raw data',
    'DATA3': 'Channel 3 raw data',
    'DATA4': 'Channel 4 raw data',
    'DATA5': 'Short Array holding measured volts/10 (EP voltage) during run',
    'DATA6': 'Short Array holding measured milliAmps trace (EP current) during run',
    'DATA7': 'Short Array holding measured milliWatts trace (Laser EP Power) during run',
    'DATA8': 'Short Array holding measured oven Temperature (polymer temperature) trace during run',
    'DATA9': 'Channel 9 processed data',
    'DATA10': 'Channel 10 processed data',
    'DATA11': 'Channel 11 processed data',
    'DATA12': 'Channel 12 processed data',
    # Prism 3100/3100-Avant may provide DATA105
    #          3130/3130-XL may provide DATA105
    # 3530/3530-XL may provide DATA105-199, 9-12, 205-299
    'DSam1': 'Downsampling factor',
    'DySN1': 'Dye set name',
    'Dye#1': 'Number of dyes',
    'DyeN1': 'Dye 1 name',
    'DyeN2': 'Dye 2 name',
    'DyeN3': 'Dye 3 name',
    'DyeN4': 'Dye 4 name',
    'DyeW1': 'Dye 1 wavelength',
    'DyeW2': 'Dye 2 wavelength',
    'DyeW3': 'Dye 3 wavelength',
    'DyeW4': 'Dye 4 wavelength',
    # 'DyeN5-N': 'Dye 5-N Name',
    # 'DyeW5-N': 'Dye 5-N Wavelength',
    'EPVt1': 'Electrophoresis voltage setting (volts)',
    'EVNT1': 'Start Run event',
    'EVNT2': 'Stop Run event',
    'EVNT3': 'Start Collection event',
    'EVNT4': 'Stop Collection event',
    'FWO_1': 'Base Order. Sequencing Analysis Filter wheel order. Fixed for 3500 at "GATC"',
    'GTyp1': 'Gel or polymer Type',
    'InSc1': 'Injection time (seconds)',
    'InVt1': 'Injection voltage (volts)',
    'LANE1': 'Lane/Capillary',
    'LIMS1': 'Sample tracking ID',
    'LNTD1': 'Length to detector',
    'LsrP1': 'Laser Power setting (micro Watts)',
    'MCHN1': 'Instrument name and serial number',
    'MODF1': 'Data collection module file',
    'MODL1': 'Model number',
    'NAVG1': 'Pixels averaged per lane',
    'NLNE1': 'Number of capillaries',
    'OfSc1': 'List of scans that are marked off scale in Collection. (optional)',
    # OvrI and OrvV are listed as "1-N", and "One for each dye (unanalyzed
    # and/or analyzed data)"
    'OvrI1': 'List of scan number indexes that have values greater than 32767 but did not '
             'saturate the camera. In Genemapper samples, this can have indexes with '
             'values greater than 32000. In sequencing samples, this cannot have '
             'indexes with values greater than 32000.',
    'OvrI2': 'List of scan number indexes that have values greater than 32767 but did not '
             'saturate the camera. In Genemapper samples, this can have indexes with '
             'values greater than 32000. In sequencing samples, this cannot have '
             'indexes with values greater than 32000.',
    'OvrI3': 'List of scan number indexes that have values greater than 32767 but did not '
             'saturate the camera. In Genemapper samples, this can have indexes with '
             'values greater than 32000. In sequencing samples, this cannot have '
             'indexes with values greater than 32000.',
    'OvrI4': 'List of scan number indexes that have values greater than 32767 but did not '
             'saturate the camera. In Genemapper samples, this can have indexes with '
             'values greater than 32000. In sequencing samples, this cannot have '
             'indexes with values greater than 32000.',
    'OvrV1': 'List of color data values found at the locations listed in the OvrI tag. '
             'There must be exactly as many numbers in this array as in the OvrI array.',
    'OvrV2': 'List of color data values found at the locations listed in the OvrI tag. '
             'There must be exactly as many numbers in this array as in the OvrI array.',
    'OvrV3': 'List of color data values found at the locations listed in the OvrI tag. '
             'There must be exactly as many numbers in this array as in the OvrI array.',
    'OvrV4': 'List of color data values found at the locations listed in the OvrI tag. '
             'There must be exactly as many numbers in this array as in the OvrI array.',
    'PDMF1': 'Sequencing Analysis Mobility file name chosen in collection',
    'RMXV1': 'Run Module XML schema version',
    'RMdN1': 'Run Module name (same as MODF)',
    'RMdX1': 'Run Module XML string',
    'RPrN1': 'Run Protocol name',
    'RPrV1': 'Run Protocol version',
    'RUND1': 'Run Started Date',
    'RUND2': 'Run Stopped Date',
    'RUND3': 'Data Collection Started Date',
    'RUND4': 'Data Collection Stopped date',
    'RUNT1': 'Run Started Time',
    'RUNT2': 'Run Stopped Time',
    'RUNT3': 'Data Collection Started Time',
    'RUNT4': 'Data Collection Stopped Time',
    'Rate1': 'Scanning Rate. Milliseconds per frame.',
    'RunN1': 'Run Name',
    'SCAN1': 'Number of scans',
    'SMED1': 'Polymer lot expiration date',
    'SMLt1': 'Polymer lot number',
    'SMPL1': 'Sample name',
    'SVER1': 'Data collection software version',
    'SVER3': 'Data collection firmware version',
    'Satd1': 'Array of longs representing the scan numbers of data points, which are flagged as saturated by data collection (optional)',
    'Scal1': 'Rescaling divisor for color data',
    'Scan1': 'Number of scans (legacy - use SCAN)',
    'TUBE1': 'Well ID',
    'Tmpr1': 'Run temperature setting',
    'User1': 'Name of user who created the plate (optional)',
    'PLOC1': 'Array of peak locations edited by user',
    'PLOC2': 'Array of peak locations as called by Basecaller',
    'PRJT1': 'SeqScape 2.0 project template name',
    'PROJ4': 'SeqScape 2.0 project name',
    'PSZE1': 'Plate size. The number of sample positions in the container. Current allowed values: 96, 384.',
    'PTYP1': 'Plate type. Current allowed values: 96-Well, 384-Well.',
    'PuSc1': 'Median pupscore',
    'QV201': 'QV20+ value',
    'QV202': 'One of "Pass", "Fail", or "Check"',
    'QcPa1': 'QC parameters',
    'QcRn1': 'Trimming and QC code',
    'QcRs1': 'QC warnings, a concatenated comma separated string',
    'QcRs2': 'QC errors, a concatenated comma separated string',
    'RGOw1': 'The name entered as the Owner of a Results Group, in the Results Group Editor. Implemented as the user name from the results group.',
    'RInj1': 'Reinjection number. The reinjection number that this sample belongs to. Not present if there was no reinjection.',
    'RNmF1': 'Raman normalization factor',
    'RevC1': 'for whether the sequence has been complemented',
    'RunN1': 'Run name (which, for 3500, is different from injection name)',
    'S/N%1': 'Signal strength for each dye',
    'SMID1': 'Polymer first installed date',
    'SMRn1': 'Number of runs (injections) processed with the current polymer (runs allowed - runs remaining)',
    'SPAC1': 'Average peak spacing used in last analysis',
    'SPAC2': 'Basecaller name - corresponds to name of bcp file.',
    'SPAC3': 'Average peak spacing last calculated by the Basecaller.',
    'SPEC1': 'Sequencing Analysis Specimen Name',
    'SVER2': 'Basecaller version number',
    'SVER4': 'Sample File Format Version String',
    'ScPa1': 'The parameter string of size caller',
    'ScSt1': 'Raw data start point. Set to 0 for 3500 data collection.',
    'SpeN1': 'Active spectral calibration name',
    'TrPa1': 'Timming parameters',
    'TrSc1': 'Trace score.',
    'TrSc2': 'One of "Pass", "Fail", or "Check"',
    'phAR1': 'Trace peak aria ratio',
    'phCH1': 'Chemistry type ("term", "prim", "unknown"), based on DYE_1 information',
    'phDY1': 'Dye ("big", "d-rhod", "unknown"), based on mob file information',
    'phQL1': 'Maximum Quality Value',
    'phTR1': 'Set Trim region',
    'phTR2': 'Trim probability',
          }     

# dictionary for unpacking tag values
_BYTEFMT = {
            1: 'b',     # byte
            2: 's',     # char
            3: 'H',     # word
            4: 'h',     # short
            5: 'i',     # long
            6: '2i',    # rational, legacy unsupported
            7: 'f',     # float
            8: 'd',     # double
            10: 'h2B',  # date
            11: '4B',   # time
            12: '2i2b', # thumb
            13: 'B',    # bool
            14: '2h',   # point, legacy unsupported
            15: '4h',   # rect, legacy unsupported
            16: '2i',   # vPoint, legacy unsupported
            17: '4i',   # vRect, legacy unsupported
            18: 's',    # pString
            19: 's',    # cString
            20: '2i',   # Tag, legacy unsupported
           }

# header structure
_HEADFMT = '>4sH4sI2H3I'

# directory data structure
_DIRFMT = '>4sI2H4I'

# to handle py3 IO
def py3_get_string(byte):
    if version_info[0] < 3:
        return byte
    else:
        return byte.decode()

def py3_get_byte(string):
    if version_info[0] < 3:
        return string
    else:
        return string.encode()

class ABI(object):
    """Class representing trace file."""
    def __init__(self, in_file, trimming=False):        
        self._handle = open(in_file, 'rb')
        try:
            self._handle.seek(0)
            if not self._handle.read(4) == py3_get_byte('ABIF'):
                raise IOError('Input is not a valid trace file')
        except IOError:
            self._handle = None
            raise
        else:
            # header data structure:
            # file type, file, version, tag name, tag number, element type code,
            # element size, number of elements, data size, data offset, handle,
            # file type, file version
            # dictionary for containing file metadata
            self.data = {}
            # dictionary for containing extracted directory data
            self.tags = {}
            self.trimming = trimming
            # values contained in file header
            self._handle.seek(0)
            header = struct.unpack(_HEADFMT, 
                     self._handle.read(struct.calcsize(_HEADFMT)))
            # file format version
            self.version = header[1]

            # build dictionary of data tags and metadata
            for entry in self._parse_header(header):
                key = entry.tag_name + str(entry.tag_num)
                self.tags[key] = entry
                # only extract data from tags we care about
                if key in EXTRACT:
                    # e.g. self.data['well'] = 'B6'
                    self.data[key] = self.get_data(key)

            self.id = self._get_file_id(in_file)
            self.name = self.get_data('SMPL1')
            self.seq = self.get_data('PBAS2')
            self.qual = ''.join([chr(ord(value) + 33) for value in self.get_data('PCON2')])
            self.qual_val = [ord(value) for value in self.get_data('PCON2')]

            if trimming:
                self.seq, self.qual, self.qual_val = map(self.trim, 
                                                        [self.seq, self.qual,
                                                        self.qual_val])

    def __repr__(self):
        """Represents data associated with the file."""
        if len(self.seq) > 10:
            seq = "{0}...{1}".format(self.seq[:5], self.seq[-5:])
            qual_val = "[{0}, ..., {1}]".format(
                      repr(self.qual_val[:5])[1:-1], 
                      repr(self.qual_val[-5:])[1:-1])
        else:
            seq = self.seq
            qual_val = self.qual_val

        return "{0}({1}, qual_val:{2}, id:{3}, name:{4})".format(
                self.__class__.__name__, repr(seq), qual_val,
                repr(self.id), repr(self.name))
    
    def _parse_header(self, header):
        """Generator for directory contents."""
        # header structure:
        # file signature, file version, tag name, tag number, 
        # element type code, element size, number of elements
        # data size, data offset, handle
        head_elem_size = header[5]
        head_elem_num = header[6]
        head_offset = header[8]
        index = 0
        
        while index < head_elem_num:
            start = head_offset + index * head_elem_size
            # added directory offset to tuple
            # to handle directories with data size <= 4 bytes
            self._handle.seek(start)
            dir_entry =  struct.unpack(_DIRFMT, 
                        self._handle.read(struct.calcsize(_DIRFMT))) + (start,)
            index += 1
            yield _TraceDir(dir_entry, self._handle)

    def _get_file_id(self, in_file):
        """Returns filename without extension."""
        return splitext(basename(in_file))[0]

    def close(sel):
        """Closes the Trace file object."""
        self._handle.close()
    

    def get_data(self, key):
        """Returns data stored in a tag."""
        return self.tags[key].tag_data

    def seq_remove_ambig(self, seq):
        """Replaces extra ambiguous bases with 'N'."""
        import re
        seq = self.seq
        return re.sub("K|Y|W|M|R|S", 'N', seq)

    def export(self, out_file="", fmt='fasta'):
        """Writes the trace file sequence to a fasta file.
        
        Keyword argument:
        out_file -- output file name (detault 'tracefile'.fa)
        fmt -- 'fasta': write fasta file, 'qual': write qual file, 'fastq': write fastq file

        """
        if out_file == "":
            file_name = self.id
            if fmt == 'fasta':
                file_name += '.fa'
            elif fmt == 'qual':
                file_name += '.qual'
            elif fmt == 'fastq':
                file_name += '.fq'
            else:
                raise ValueError('Invalid file format: {0}.'.format(fmt))
        else:
            file_name = out_file
        
        if fmt == 'fasta':
            contents = '>{0} {1}\n{2}\n'.format(
                        self.id, 
                        self.name, 
                        self.seq)
        elif fmt == 'qual':
            contents = '>{0} {1}\n{2}\n'.format(
                        self.id, 
                        self.name, 
                        ' '.join(map(str, self.qual_val)))
        elif fmt == 'fastq':
            contents = '@{0} {1}\n{2}\n+{0} {1}\n{3}\n'.format(
                        self.id, 
                        self.name, 
                        self.seq, ''.join(self.qual))

        with open(file_name, 'w') as out_file:
            out_file.writelines(contents)

    def trim(self, seq, cutoff=0.05):
        """Trims the sequence using Richard Mott's modified trimming algorithm.
        
        Keyword argument:
        seq -- sequence to be trimmed
        cutoff -- probability cutoff value

        Trimmed bases are determined from their segment score, ultimately
        determined from each base's quality values. 
        
        More on:
        http://www.phrap.org/phredphrap/phred.html
        http://www.clcbio.com/manual/genomics/Quality_trimming.html
        """
        # set flag for trimming
        start = False
        # set minimum segment size
        segment = 20
        trim_start = 0
        
        if len(seq) <= segment:
            raise ValueError('Sequence can not be trimmed because \
                             it is shorter than the trim segment size')
        else:
            # calculate probability back from formula used
            # to calculate phred qual values
            score_list = [cutoff - (10 ** (qual/-10.0)) for 
                         qual in self.qual_val]

            # calculate cummulative score_list
            # if cummulative value < 0, set to 0
            # first value is set to 0 (assumption: trim_start is always > 0)
            running_sum = [0]
            for i in range(1, len(score_list)):
                num = running_sum[-1] + score_list[i]
                if num < 0:
                    running_sum.append(0)
                else:
                    running_sum.append(num)
                    if not start:
                        # trim_start = value when cummulative starts to be > 0
                        trim_start = i
                        start = True

            # trim_finish = index of the highest cummulative value,
            # marking the segment with the highest cummulative score 
            trim_finish = running_sum.index(max(running_sum)) 

            return seq[trim_start:trim_finish]

    def chromatogram(self, ax=None, xlim=None, highlight=None, seq=True, trim_threshold=50, legend=True):
        """
        Draw the trace.
        Inputs:
            ax: Matplotlib axis object. Axis to use. If not defined, it is made.
            xlim: two-tuple. The range of sequence to plot. Default is full sequence
            highlight: two-tuple. Highlight the range specified in this region.
            seq: Boolean. Write the sequence on top of plot.
            trim_threshold: Integer. Cut the sequence from the right if all signals are less than this value.
            legend: Boolean. Draw the legend
        """
        
        # If we don't have an axis... Make it
        if ax is None:
            fig, axes = utils.figure.subplots(ncols=1, nrows=1, figsize=(10, 3))
            ax = axes[0]
        #fi
    
        # Re-index the x-axis to correspond to each base.
        xrange = [ 1 ]
        ploc2 = self.data['PLOC2']
        for i in range(len(ploc2)-1):
            curr_peak = ploc2[i]
            next_peak = ploc2[i+1]
            steps = next_peak - curr_peak # Number of measurements inbetween peaks
            step_size = 1 / steps
            xrange.extend([ xrange[-1] + step_size * (i+1) for i in range(steps) ])
        #efor
        maxy = 0
        
        # Plot the chromatograms
        channels = [ 'DATA%d' % c for c in [9,10,11,12] ]
        for c, color, legend in zip(channels, 'krgb', self.data['FWO_1']):
            ax.plot(xrange, self.data[c][ploc2[0]:][:len(xrange)], label=legend, c=color)
            maxy = max(maxy, max(self.data[c][:len(xrange)]))
        #efor
        
        # Plot the vertical lines to indicate position
        for x in range(10, int(xrange[-1]), 10):
            ax.plot([x, x], [0, maxy], c='grey', linestyle=':')
        #efor
        
        # Determine where we should cut thr sequence
        right_trim = len(ploc2)
        for i, peak in list(enumerate(ploc2))[::-1]:
            if max([self.data[c][peak] for c in channels]) >= trim_threshold:
                break
            #fi
            right_trim = i
        #efor
        
        # Set the range we wanted to plot.
        if xlim is None:
            xlim = (ploc2[0], int(xrange[-1]))
        #fi
        xlim = (xlim[0], min(xlim[1], right_trim))
        ax.set_xlim(xlim[0], xlim[1])
        
        # Hide the top spine
        ax.spines['top'].set_visible(False)
        
        # Write the sequence
        colors = dict(zip(self.data['FWO_1'], 'krgb'))
        if seq:
            for i, base in enumerate(self.seq[xlim[0]:xlim[1]-1]):
                ax.text(i+xlim[0]+1, maxy, base, color=colors[base], ha='center')
            #efor
        #fi
        
        # Draw the legend
        if legend:
            ax.legend()
        #fi

        if highlight is not None:
            x1, x2 = highlight
            rect = patches.Rectangle((x1,0), (x2-x1), maxy, linewidth=1, edgecolor='#a6bddb',facecolor='#a6bddb', alpha=100)
            ax.add_patch(rect)
        #fi

        # Only show integer positions
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        
        return ax
    #edef

class _TraceDir(object):
    """Class representing directory content."""
    def __init__(self, tag_entry, handle):
        self.tag_name = py3_get_string(tag_entry[0])
        self.tag_num = tag_entry[1]
        self.elem_code = tag_entry[2]
        self.elem_size = tag_entry[3]
        self.elem_num = tag_entry[4]
        self.data_size = tag_entry[5]
        self.data_offset = tag_entry[6]
        self.data_handle = tag_entry[7]
        self.tag_offset = tag_entry[8]

        # if data size is <= 4 bytes, data is stored inside the directory
        # so offset needs to be changed
        if self.data_size <= 4:
            self.data_offset = self.tag_offset + 20

        self.tag_data = self._unpack(handle)

    def __repr__(self):
        """Represents data associated with a tag."""
        summary = ['tag_name: {0}'.format(repr(self.tag_name))]
        summary.append('tag_number: {0}'.format(repr(self.tag_num)))
        summary.append('elem_code: {0}'.format(repr(self.elem_code)))
        summary.append('elem_size: {0}'.format(repr(self.elem_size)))
        summary.append('elem_num: {0}'.format(repr(self.elem_num)))
        summary.append('data_size: {0}'.format(repr(self.data_size)))
        summary.append('data_offset: {0}'.format(repr(self.data_offset)))
        summary.append('data_handle: {0}'.format(repr(self.data_handle)))
        summary.append('tag_offset: {0}'.format(repr(self.tag_offset)))
        summary.append('tag_data: {0}'.format(repr(self.tag_data)))
       
        return '\n'.join(summary)

    def _unpack(self, handle):
        """Returns tag data"""
        if self.elem_code in _BYTEFMT:
            
            # because ">1s" unpacks differently from ">s"
            num = '' if self.elem_num == 1 else str(self.elem_num)
            fmt = "{0}{1}{2}".format('>', num, _BYTEFMT[self.elem_code])
            start = self.data_offset
    
            handle.seek(start)
            data = struct.unpack(fmt, handle.read(struct.calcsize(fmt)))
            
            # no need to use tuple if len(data) == 1
            if self.elem_code not in [10, 11] and len(data) == 1:
                data = data[0]

            # account for different data types
            if self.elem_code == 2:
                return py3_get_string(data)
            elif self.elem_code == 10:
                return datetime.date(*data)
            elif self.elem_code == 11:
                return datetime.time(*data)
            elif self.elem_code == 13:
                return bool(data)
            elif self.elem_code == 18:
                return py3_get_string(data[1:])
            elif self.elem_code == 19:
                return py3_get_string(data[:-1])
            else:
                return data
        else:
            return None
