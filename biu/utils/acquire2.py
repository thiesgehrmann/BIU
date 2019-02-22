# Utilities to acquire files

from . import msgUtils as msg
from . import exeUtils as exe
from . import fsUtils as fs
from . import pyUtils as py
from .. import settings
from .. import ops

import os
from collections import namedtuple
import hashlib
import ftplib

###############################################################################

class AcquireFile(object):
    """
    An object to help manage filepaths and files before they are actually created.
    """
    def __init__(self, dirname, basename):
        """
        Specify a file in the Acquire pipeline
        """
        
        self.dirname  = dirname
        self.basename = basename
        
        if not(isinstance(self.basename, str) and self.basename):
            raise ValueError("No basename is specified.")
        #fi
    #edef
    
    def set_dirname(self, dirname):
        """
        Set the dirname for this file.
        """
        self.dirname = os.path.abspath(os.path.expanduser(dirname))
    #edef
    
    def provis_path(self, dirname):
        """
        Get the provisional path, given a specified dirname.
        SUBJECT TO CHANGE
        """
        return '%s/%s' % (dirname, self.basename)
    #edef
    
    @property
    def path(self):
        if self.dirname is None:
            raise RuntimeError("No dirname has been specified.")
        #fi
        return '%s/%s' % (self.dirname, self.basename)
    #edef
    
    @property
    def marker_path(self):
        """
        The path to the marker file that verifies that it was correctly processed.
        """
        return self.path + '.__biu.acquire_exists__'
    #edef
    
    @property
    def exists(self):
        """
        Determine if the file exists.
        """
        return os.path.exists(self.path) and os.path.exists(self.marker_path)
    #edef
    
    def set_exists(self):
        """
        Set the exists flag for this file
        """
        if os.path.exists(self.path):
            fs.touchFile(self.marker_path)
        else:
            raise RuntimeError("Cannot set exists state for file that does not exist. '%s'." % self.path)
        #fi
    #edef
        
    def __str__(self):
        return '%s/%s' % ('__UNDEFINED__' if self.dirname is None else self.dirname, self.basename)
    #edef
    
    def __repr__(self):
        return str(self)
    #edef
#eclass

###############################################################################
    
class AcquireFixedFile(AcquireFile):
    """
    For pre-existing files, we don't want:
     * to allow the dirnames to be changed
     * the marker_path to be created
     
    This class handles this case
    """
    
    def __init__(self, file=None, directory=None):
        """
        Initialize an AcquireFixedFile object.
        
        If you have a file:
          AcquirePreExistingFile('./my_file.txt')
        If you have a directory:
          AcquirePreExistingFile(directory='./my_directory/')
        """
        dirname  = None
        basename = None
        if file is not None:
            file = os.path.abspath(os.path.expanduser(file))
            dirname  = os.path.dirname(file)
            basename = os.path.basename(file)
        elif directory is not None:
            directory = os.path.abspath(os.path.expanduser(directory))
            dirname, basename = os.path.split(directory)
        #fi
        
        super(AcquireFixedFile, self).__init__(dirname, basename)
    #edef

    @property
    def marker_path(self):
        return self.path
    #edef
    
    def set_dirname(self, dirname):
        pass
    #edef
    
    def provis_path(self, dirname):
        """
        Get the provisional path, given a specified dirname.
        SUBJECT TO CHANGE
        """
        return '%s/%s' % (self.dirname, self.basename)
    #edef
#eclass

###############################################################################

class AcquireFinalFile(AcquireFile):
    """
    For final files, we don't want:
     * to allow the dirnames to be changed (The filename should be FIXED)
     
    This class handles this case
    """
    
    def __init__(self, file=None, directory=None):
        """
        Initialize an AcquireFinalFile object.
        
        If you have a file:
          AcquireFinalFile('./my_file.txt')
        If you have a directory:
          AcquireFinalFile(directory='./my_directory/')
        """
        dirname  = None
        basename = None
        if file is not None:
            file = os.path.abspath(os.path.expanduser(file))
            dirname  = os.path.dirname(file)
            basename = os.path.basename(file)
        elif directory is not None:
            directory = os.path.abspath(os.path.expanduser(directory))
            dirname, basename = os.path.split(directory)
        #fi
        
        super(AcquireFinalFile, self).__init__(dirname, basename)
    #edef
    
    def set_dirname(self, dirname):
        pass
    #edef
    
    def provis_path(self, dirname):
        """
        Get the provisional path, given a specified dirname.
        SUBJECT TO CHANGE
        """
        return '%s/%s' % (self.dirname, self.basename)
    #edef
#eclass

###############################################################################

class AcquireStep(object):
    """
    A step in an acquire pipeline.
    """
    __slots__ = [ 'name', 'inputs', 'output', 'action' ]
    def __init__(self, name, inputs, output_file, action):
        """
        Initialize an AcquireStep object
        
        parameters:
        -----------
        name: String. The name of the step
        inputs: A list of AcquireFile objects. The required files for this step.
        output: An AcquireFile object. This will be the output file.
        action: A function. format: action(inputs, output) -> Boolean
                inputs is a list of file paths for the inputs
                output is a file path for the output.
                
        """
        
        if not isinstance(name, str):
            raise ValueError("Specified name is not string.")
        elif hasattr(inputs, '__iter__') and (not all([ isinstance(a, AcquireFile) for a in inputs ])):
            raise ValueError("Specified input is not a list of AcquireFile objects.")
        elif not isinstance(output_file, AcquireFile):
            raise ValueError("Specified output_file is not an AcquireFile object.")
        elif not hasattr(action, '__call__'):
            raise ValueError("The provided action should be callable")
        #fi
        
        self.name   = name
        self.inputs = inputs
        self.output = output_file
        self.action = action
    #edef
    
    def do(self, where, redo):
        """
        Perform the action of this AcquireStep object.
        parameters:
        -----------
        where: String. Where to perform the pipeline.
        redo:  boolean. Re-do this step, even if it has already been completed?
        """
        
        for i in self.inputs:
            i.set_dirname(where)
            if not i.exists:
                raise RuntimeError("One of the inputs for step '%s' does not exist: '%s'." % (self.name, i.path))
            #fi
        #efor
        
        self.output.set_dirname(where)
        
        if (not self.output.exists) or redo:
            fs.mkdirname(self.output.path) # Make the folder for the file.
            status = self.action([ i.path for i in self.inputs], self.output.path)
            if status != 0:
                raise RuntimeError("Execution failed at step '%s'." % self.name)
            else:
                self.output.set_exists()
            #fi
        #fi
    #edef
    
    def path(self, where):
        """
        The provisional path for the output of this step.
        """
        return self.output.provis_path(where)
    #edef
#eclass 

###############################################################################

class Acquire2(object):
    """
    Define a pipeline to acquire and process a file.
    
    Current functions:
    ---------------------
    
    Initial file acquisition:
     * local: from local file
     * touch: create local file
     * curl:  from the web. 
     * wget:  from the web. Note: doesn't work on mac os...
     * ftp:   from the web via FTP
     * lftp:  from the web via LFTP
     
    Process the file:
     * unzip
     * gunzip
     * untar
     * select : From a directory (e.g. after unzipping), select a file
     
     * sort
     * tabix
     
     * bgzip
     * gzip
     * bzip
     
    Combine pipelines:
     * merge (cat|zcat)
     
    Finalize:
     * finalize: Finalize the file's location (either through a hard copy, or a symbolic link)
     
     
    Example pipeline:
    --------------------------
    
    a = Acquire2().curl("http://zzz.bwh.harvard.edu/plink/dist/plink-1.07-i686.zip")
    a = a.unzip().select("plink-1.07-i686/README.txt").sort().finalize('./test.README.md')
    
    file_path = a.acquire().path
    
    """
    
    AcquireStep = AcquireStep
    AcquireFile = AcquireFile
    AcquireFinalFile = AcquireFinalFile
    AcquireFixedFile = AcquireFixedFile
    STATUS_SUCCESS = 0
    STATUS_FAILURE = 1
    
    def __init__(self, file=None, where=settings.getDownloadDir(), redo=False, steps=None):
        """
        Initialize an Acquire2 object
        parameters:
        -----------
        file: Specify immediately a local file. ( Acquire2(myfile) is identical to Acquire2().local(myfile) )
        where: Specify a location where the pipeline should take place
        redo:  Redo the pipeline, regardless of whether it should be repeated or not.
        steps: For internal use.
        """
        self.where = os.path.abspath(os.path.expanduser(where))
        self.redo  = redo

        steps = [] if steps is None else steps
        for step in steps:
            if not isinstance(step, AcquireStep):
                raise ValueError("Steps parameter must be list of AcquireSteps")
            #fi
        #efor
        
        self.steps = steps
        
        if (file is not None):
            step = AcquireStep("Local", [], AcquireFixedFile(file), lambda i, o: self.STATUS_SUCCESS)
            self.steps.append(step)
        #fi
    #edef
    
    def copy(self):
        """
        Make a copy of the current Acquire2 pipeline
        """
        return Acquire2(where=self.where, redo=self.redo, steps=self.steps)
    #edef
    
    @property
    def output(self):
        """
        Return the output AcquireFile of this pipeline.
        """
        
        if len(self.steps) < 1:
            raise ValueError("No steps have been specified yet...")
        #fi
        
        return self.steps[-1].output
    #edef
    
    @property
    def exists(self):
        """
        Return whether or not the file exists.
        """
        return os.path.exists(self.path)
    #edef
    
    @property
    def path(self):
        """
        Return the path of the output file of this pipeline.
        """
        
        if len(self.steps) < 1:
            raise ValueError("No steps have been specified yet...")
        #fi
        
        return self.steps[-1].path(self.where)
    #edef
    
    def acquire(self):
        """
        Run the acquire pipeline.
        
        Returns the output AcquireFile object for this pipeline.
        """
        if self.output.exists and not self.redo:
            return self.output
        #fi
        
        for step in self.steps:
            step.do(self.where, self.redo)
        #efor
        
        return self.output
    #edef
    
    def __str__(self):
        dstr = "Acquire object.\n"
        dstr += ' Re-do steps: %s\n' % ('yes' if self.redo else 'no')
        dstr += " Current steps:\n"
        for step in self.steps:
            dstr += '  * %s -> %s\n' % (step.name, step.path(self.where))
        #efor
        return dstr
    #edef
    
    def __repr__(self):
        return str(self)
    #edef
    
    def add_step(self, step):
        """
        Add a step to the current pipeline.
        parameters:
        -----------
        step: AcquireStep object
        
        Returns:
        --------
        A new Acquire2 object with the added step.
        """
        
        if not isinstance(step, AcquireStep):
            raise ValueError("step must be an AcquireStep object!")
        #fi
        
        return Acquire2(where=self.where, redo=self.redo, steps=self.steps + [step])
    #edef
    
    def set_where(self, where):
        """
        Change the location of where a file should be set.
        """
        return Acquire2(where=where, redo=self.redo, steps=self.steps)
    #edef
    
    def set_redo(self, redo):
        """
        Specify that we should redo the whole pipeline.
        """
        return Acquire2(where=self.where, redo=redo, steps=self.steps)
    #edef
    
    def local(self, file):
        """
        Specify a local file.
        """
        step = AcquireStep("Local", [], AcquireFixedFile(file), lambda i,o: self.STATUS_SUCCESS)
        return self.add_step(step)
    #edef
    
    def touch(self, file):
        """
        Use a local file with a specific name. Touch the file.
        """
        step = AcquireStep("Touch", [], AcquireFile(file), lambda i,o: fs.touchFile(i[0].path))
        return self.add_step(step)
    #edef
    
    def finalize(self, finalpath, ln=False):
        """
        Move the file to a final location.
        parameters:
        -----------
        ln: Boolean. Make a symbolic link instead of a physical copy.
        """
        from shutil import copy
        from os import symlink
        
        def _rm_if_exists(f):
            if os.path.exists(f):
                if os.path.isfile(f):
                    os.remove(f)
                elif os.path.isdir(f):
                    import shutil
                    shutil.rmtree(f)
                elif os.path.islink(f):
                    import pathlib
                    pathlib.rmlink(f)
                else:
                    raise ValueError("I don't know what this file is... '%s'" % f)
                #fi
            #fi
        #edef
        
        def _copy(i, o):
            _rm_if_exists(o)
            copy(i[0],o)
            return self.STATUS_SUCCESS
        #edef
        
        def _ln(i,o):
            _rm_if_exists(o)
            symlink(i[0],o)
            return self.STATUS_SUCCESS
        #edef
        fn = _ln if ln else _copy
            
        step = AcquireStep("Finalize%s" % ('(symlink)' if ln else ''), [self.output],
                           AcquireFinalFile(finalpath), fn)
        
        return self.add_step(step)
    #edef
    
    def curl(self, url, cookieURL=None, username=None, password=None, ext=None):
        """
        curl: Download a file using curl
        Inputs: url: The URL to retrieve
                cookieURL : The URL of a login page
                username: If logging in, the username
                password: If logging in, the password
                ext: Optionally specify a file extension
        """
        
        ext         = '' if ext is None else ('.' + ext)
        curl_hash   = ops.lst.hash([ url, cookieURL, username, password ]) + ext
        output_file = AcquireFile(dirname=None, basename=curl_hash)
        
        def _curl(output, url, cookieURL=None, username=None, password=None):
            """
            curl: Download a file using curl
            Inputs: url: The URL to retrieve
                    cookieURL : The URL of a login page
                    username: If logging in, the username
                    password: If logging in, the password
                    ext: Optionally specify a file extension
            Output: Acquire object
            """

            cookieFile = None
            if (cookieURL is not None) or (username is not None) or (password is not None):
                cookieFile = self.__fileName + '.cookie'
                p = exe.runCommand('curl -b %s --user "%s:%s" "%s"' % (cookieFile, username, password, cookieURL))
                if p != 0:
                    msg.error("Could not set cookie... for site '%s'" % cookieURL)
                    return p
                #fi
            #fi
            p = exe.runCommand("curl -L %s '%s' > '%s'" % ( '-c "%s"' % cookieFile if cookieFile is not None else  '', url, output), shell=True, verbose=True)
            return p
        #edef
        
        step = AcquireStep("Curl(%s)" % url, [], output_file, lambda i, o: _curl(o, url, cookieURL, username, password))
        return self.add_step(step)
    #edef
    
    def _ftp(self, server, location, username=None, password=None, ext=None):
        """
        ftp: Download a file with FTP
        Inputs: server: The server to access
                location: The path to the file on the server
                username: The username to access the server (if None, default is used)
                password: The password to access the erver
                ext: Optionally specify a file extension
        Output: Acquire object
        """
        ext         = '' if ext is None else ('.' + ext)
        ftp_hash    = ops.lst.hash([ server, location, username, password ]) + ext
        output_file = AcquireFile(dirname=None, basename=ftp_hash)

        def _ftp(output_file, server, location, username, password):
            conn = ftplib.FTP(server)
            if (username is not None) and (password is not None):
                conn.login(username, password)
            else:
                conn.login()
            #fi
            p = conn.retrbinary(location, open(output_filename, 'wb').write)
            p = int(result.split(' ')[0])
            if p == 226:
                return 0
            #fi
            return p
        #edef
        
        step = AcquireStep("FTP(%s/%s)" % (server, location), [], output_file,
                           lambda i,o: _ftp(o, server, location, username, password))
        return self.add_step(step)
    #edef

    def lftp(self, server, location, username, password, ext=None):
        """
        lft: Download a file from an ftp server (but for example with sftp access)
        Inputs: server: The server to access
                location: The path to the file on the server
                username: The username to access the server (if None, default is used)
                password: The password to access the erver
                ext: Optionally specify a file extension
        Output: Acquire object
        """
        ext         = '' if ext is None else ('.' + ext)
        ftp_hash    = ops.lst.hash([ server, location, username, password ]) + ext
        output_file = AcquireFile(dirname=None, basename=ftp_hash)

        if not exe.exists('lftp'):
            raise RuntimeError("'lftp' is not installed. Please install in order to continue.")
            return 1
        #fi

        def _lftp(output_file, server,location, username, password):
            cmd = "echo -en  'open \"%s\"\\nuser \"%s\" \"%s\"\\ncat \"%s\"' | lftp > '%s'" % (server, username, password, location, output_file)
            p = exe.runCommand(cmd, shell=True, verbose=True)
            return p
        #edef
        
        step = AcquireStep("LFTP(%s/%s)" % (server, location), [], output_file,
                   lambda i,o: _lftp(o, server, location, username, password))
        return self.add_step(step)
    #edef


    def wget(self, url, ext=None):
        """
        wget: Download a file with wget
        Inputs: url: URL to retrieve
                ext: Optionally specify a file extension
        Output: Acquire object
        """
        ext         = '' if ext is None else ('.' + ext)
        ftp_hash    = ops.lst.hash([ url ]) + ext
        output_file = AcquireFile(dirname=None, basename=ftp_hash)

        def _wget(output_file, url):
            cmd = "wget -O '%s' '%s'" % ( output_file, url )
            p = exe.runCommand(cmd, verbose=True)
            return p
        #edef
        step = AcquireStep("WGET(%s)" % (url), [], output_file,
                           lambda i,o: _wget(o, url))
        return self.add_step(step)
    #edef
    
    def unzip(self, select=None):
        """
        Unzip a ZIP file.
        parameters:
        ------------
        select: String. Select a file from the zip file.
        """
        inputs      = [self.output]
        unzip_dir   = AcquireFile(inputs[-1].dirname, inputs[-1].basename + '.__unzipped__')
        
        def _unzip(input_zip, output_dir):
            p = exe.runCommand("unzip -o -d '%s' '%s'" % (output_dir, input_zip), verbose=True)
            return p
        #edef
        
        step = AcquireStep("Unzip", inputs, unzip_dir, lambda i,o: _unzip(i[0], o))
        added = self.add_step(step)
        
        if select is not None:
            return added.select(select)
        else:
            return added
        #fi
    #edef
    
    def untar(self, select=None):
        """
        Untar a TAR file.
        parameters:
        ------------
        select: String. Select a file from the TAR file.
        """
        inputs      = [self.output]
        output_dir   = AcquireFile(inputs[-1].dirname, inputs[-1].basename + '.__untarred__')
        
        def _unzip(input_tar, output_dir):
            p = exe.runCommand("tar -xf '%s' -C '%s'" % (input_tar, output_dir), verbose=True)
            return p
        #edef
        
        step = AcquireStep("Untar", inputs, output_dir, lambda i,o: _untar(i[0], o))
        added = self.add_step(step)
        
        if select is not None:
            return added.select(select)
        else:
            return added
        #fi
    #edef
    
    def gunzip(self, select=None):
        """
        gunzip a file
        parameters:
        ------------
        select: String. Select a file from the gzipped file.
        """
        inputs      = [self.output]
        output_dir   = AcquireFile(inputs[-1].dirname, inputs[-1].basename + '.__gunzipped__')
        
        def _gunzip(input_gzip, output_dir):
            p = exe.runCommand("gunzip < '%s' > '%s'" % (input_gzip, output_dir), verbose=True, shell=True)
            return p
        #edef
        
        step = AcquireStep("Gunzip", inputs, output_dir, lambda i,o: _gunzip(i[0], o))
        added = self.add_step(step)
        
        if select is not None:
            return added.select(select)
        else:
            return added
        #fi
    #edef
    
    def select(self, select):
        """
        From a directory, select a file
        parameters:
        ------------
        select: String. Select a file from the directory
        """
        inputs      = [self.output]
        select_file = inputs[-1].basename + '/' + select
        
        output = AcquireFile(inputs[-1].dirname, select_file)
        
        step = AcquireStep("Select", inputs, output, lambda i,o: self.STATUS_SUCCESS)
        return self.add_step(step)
    #edef


    def tabix(self, file_type=None, seq=0, start=1, end=2):
        """
        Make a tabix index of a bgzipped file.
        parameters:
        -----------
        Specify either a filetype, or specific columns
        """
        inputs   = [self.output]
        out_file = self.output.basename + '.tbi'
        output        = AcquireFile(inputs[-1].dirname, out_file)

        def _tabix(infile, fileType, seq, start, end):
            if fileType is not None:
                cmd = "tabix -p %s '%s'" % (fileType, infile)
            elif isinstance(seq, int) and isinstance(start, int) and isinstance(end, int):
                cmd = "tabix -s %d -b %d -e %d %s" % (seq, start, end, infile)
            else:
                raise ValueError("columns indexes not correctly specified.")
            #fi
            p = exe.runCommand(cmd, verbose=True)
            return p
        #edef
        
        step = AcquireStep("Tabix", inputs, output, lambda i,o: _tabix(i[0], file_type, seq, start, end))
        return self.add_step(step)
    #edef

    def _sort(self, options=None):
        """
        Sort the contents of a file using the unix sort command
        parameters:
        -----------
        options: Additional parameters for the unix sort command
        """
        inputs   = [self.output]
        out_file = self.output.basename +  '.__sorted__'
        output  = AcquireFile(inputs[-1].dirname, out_file)
        
        def _sort(input_file, output_file, options):
            p = exe.runCommand("sort %s < '%s' > '%s'" % ('' if options is None else options, input_file, output_file), shell=True)
            return p
        #edef
        
        step = AcquireStep("Sort(%s)" % options, inputs, output, lambda i,o: _sort(i[0], o, options))
        return self.add_step(step)
    #edef

    def bgzip(self):
        """
        bgzip a file.
        """
        inputs   = [self.output]
        out_file = self.output.basename + '.bgz'
        output   = AcquireFile(inputs[-1].dirname, out_file)
    
        def _bgzip(infile, outfile):
            p = exe.runCommand("bgzip < '%s' > '%s'" % (infile, outfile), verbose=True, shell=True)
            return p
        #edef
        
        step = AcquireStep("Bgzip", inputs, output, lambda i,o: _bgzip(i[0], o))
        return self.add_step(step)
    #edef

    def gzip(self):
        """
        gzip a file.
        """
        inputs   = [self.output]
        out_file = self.output.basename + '.gz'
        output   = AcquireFile(inputs[-1].dirname, out_file)
    
        def _gzip(infile, outfile):
            p = exe.runCommand("gzip < '%s' > '%s'" % (infile, outfile), verbose=True, shell=True)
            return p
        #edef
        
        step = AcquireStep("Gzip", inputs, output, lambda i,o: _gzip(i[0], o))
        return self.add_step(step)
    #edef

    def bzip(self):
        """
        bzip a file.
        """
        inputs   = [self.output]
        out_file = self.output.basename + '.bz'
        output   = AcquireFile(inputs[-1].dirname, out_file)
    
        def _bzip(infile, outfile):
            p = exe.runCommand("bzip < '%s' > '%s'" % (infile, outfile), verbose=True, shell=True)
            return p
        #edef
        
        step = AcquireStep("Bzip", inputs, output, lambda i,o: _bzip(i[0], o))
        return self.add_step(step)
    #edef
    
    def cmd(self, cmd, placeholder=None):
        """
        Call a shell function on the file.
        parameters:
        ------------
        cmd: The shell command to execute.
        placeholder: A string that will be replaced with the filename.
        
        if placeholder is None, then the following command will be executed:
        `cmd filename`
        
        otherwise:
        
        cmd("cat {} | tee './test.tee'", '{}')
            will result in:
        `cat filename | tee './test.tee'"`
        """
        inputs   = [self.output]
        out_file = self.output.basename + '.cmd'
        output   = AcquireFile(inputs[-1].dirname, out_file)

        def _cmd(infile, outfile, placeholder):
            if placeholder is not None:
                cmd = cmd.replace(placeholder, infile)
            else:
                cmd = cmd + ' ' + infile
            #fi
            p = exe.runCommand("%s > '%s'" % (cmd, outfile), shell=True, verbose=True)
            return p
        #edef
        
        step = AcquireStep("CMD", inputs, output, lambda i,o: _cmd(i[0], o, cmd))
        return self.add_step(step)
    #edef

    def func(self, function):
        """
        Call a python function on the pipeline.
        
        parameters:
        -----------
        function: The function to be run.
            Must return 0 (Or Acquire.STATUS_SUCCESS) in success.
            should have form:
            
            def myfunc(inputs, output):
                '''
                inputs: list of input file paths. Usually just one.
                output: the output file path.
                
                Your function should create this output file.
                '''
                
                return Acquire2.STATUS_SUCCESS
            #edef
        """
        inputs   = [self.output]
        out_file = self.output.basename + '.func'
        output   = AcquireFile(inputs[-1].dirname, out_file)

        if not hasattr(function, '__call__'):
            raise ValueError("The provided function should be callable.")
        
        step = AcquireStep("Func", inputs, output, function)
        return self.add_step(step)
    #edef
    
    def merge(self, acquire_objects, method="cat", pass_where=True, pass_redo=True):
        """
        Merge the output of several acquire pipelines.
        parameters:
        -----------
        acquire_objects: A list of acquire objects
        method: String. How to merge the outputs of the pipelines [cat|zcat]
        pass_where: boolean. Pass the current object's where parameter to each provided acquire object
        pass_redo:  boolean. Pass the current object's redo parameter to each provided acquire object
        
        NOTE: the where and redo parameters of the current pipeline must be set BEFORE running merge, otherwise it isn't passed into the provided acquire objects (unless pass_where and pass_redo are specified)
        """
        
        new_acquire = Acquire2()
        
        if pass_where:
            acquire_objects = [ ao.set_where(self.where) for ao in acquire_objects ]
        #fi
        if pass_redo:
            acquire_objects = [ ao.set_redo(self.redo) for ao in acquire_objects ]
        #fi
        
        inputs = [ ao.output for ao in acquire_objects ]
        
        output = AcquireFile(self.where, ops.lst.hash([ ao.path for ao in acquire_objects ]) + '.__merged__')
        
        def _merge(inputs, output, acquire_objects, method):
            
            inputs = [ ao.acquire().path for ao in acquire_objects ]
            
            if method == 'cat':
                cmd = "cat '%s' > '%s'" % ("' '".join(inputs), output)
            elif method == 'zcat':
                cmd = "cat '%s' | zcat > '%s'" % ("' '".join(inputs), output)
            else:
                raise NotImplementedError("Method '%s' is not implemented for merge" % method)
            #fi
            p = exe.runCommand(cmd, verbose=True, shell=True)
            return p
        #edef
        
        step = AcquireStep("Merge(%s %d pipelines)" % (method, len(acquire_objects)), [], output, lambda i,o: _merge(i,o,acquire_objects, method))
        return self.add_step(step)
    #edef

#eclass
        
        
        