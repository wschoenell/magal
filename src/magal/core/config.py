import ConfigParser
import ast
import os

__author__ = 'william'

class MagalConfig(ConfigParser.RawConfigParser):
    """
    This class abstracts the configuration files and checks for all the MagAl scripts.
    """

    def __init__(self, config_file, config_type):
        """
        Init a MagAl config instance.

        Parameters
        ----------
        config_file : string
            Path to the configuration file

        config_type : string
            Means which script is calling this object.

        See Also
        --------
        ConfigParser.RawConfigParser
        """
        ConfigParser.RawConfigParser.__init__(self)
        self.read(config_file)
        if type == 'library':
            self._load_library_vars()

    def _get_inputfiles(self, section, files):
        """
        Parameters
        ----------
        section : string
            config section to get the filenames.

        files : list
            Filenames list.

        Returns
        -------
        input_files : tuple
            All the asked inputfiles.

        Raises
        ------
        """
        for f_ in files:
            fname = os.path.expandvars(self.get(section, f_))

    def _load_library_vars(self):
        """
        Loads magal_library related configfile parameters.
        """
        libfile =
        library_type = os.path.expandvars(self.get('LibraryGeneral', 'library_type'))
        filter_file = os.path.expandvars(self.get('LibraryGeneral', 'filter_file'))

        filter_id = self.get('LibraryGeneral', 'filterset')  #TODO: Check if it is valid.
        z_from = self.getfloat('LibraryGeneral', 'z_from')
        z_to = self.getfloat('LibraryGeneral', 'z_to')
        z_step = self.getfloat('LibraryGeneral', 'z_step')
        try:
            cosmology = ast.literal_eval(self.get('LibraryGeneral', 'cosmology'))
            from astropy.cosmology import FlatLambdaCDM
            cosmo = FlatLambdaCDM(H0=cosmology['H0'], Om0=cosmology['Om0'])
        except ConfigParser.NoOptionError:
            print 'Adopting the WMAP9 comology:'
            from astropy.cosmology import WMAP9 as cosmo
            print cosmo.__doc__

        try:
            allow_overwrite = self.getboolean('LibraryGeneral', 'allow_overwrite')
        except ConfigParser.NoOptionError:
            allow_overwrite = False
        if os.path.exists(libfile):
            if allow_overwrite:
                try:
                    os.unlink(libfile)
                except:
                    print u'Cannot DELETE file {0:s}.'.format(libfile)
                    # return 1
            else:
                print u'File {0:s} already exists.'.format(libfile)
                # return 1


        # 1 - Two-exponential library parameters:
        if library_type == 'two_exp':
            inp_file = os.path.expandvars(self.get('LibraryParameters', 'input_file'))
            basedir = os.path.expandvars(self.get('LibraryParameters', 'bases_dir'))
            basefile = os.path.expandvars(self.get('LibraryParameters', 'base_file'))
            fraction_type = self.get('LibraryParameters', 'fraction_type')  # light or mass
            try:
                lambda_norm = self.getfloat('LibraryParameters', 'lambda_norm')  # mass normalization wl
            except ConfigParser.NoOptionError:
                lambda_norm = 4020
            LibModel = LibraryModel(library_type, inp_file, basedir, basefile, fraction_type, lambda_norm)  # Init
        # 2 - STARLIGHT library parameters:
        elif library_type == 'starlight_sdss':
            type = 'starlight_sdss'
            inp_file = os.path.expandvars(self.get('LibraryParameters', 'input_file'))
            tables_dir = os.path.expandvars(self.get('LibraryParameters', 'tables_dir'))
            input_dir = os.path.expandvars(self.get('LibraryParameters', 'input_dir'))
            output_dir = os.path.expandvars(self.get('LibraryParameters', 'output_dir'))
            LibModel = LibraryModel(type, inp_file, tables_dir, input_dir, output_dir)

        #### -- ####

