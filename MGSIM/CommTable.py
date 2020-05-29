from Utils import _table
import sys, os
import re
import numpy as np

class CommTable(_table):
    """pandas DataFrame subclass
    """
    
    def __init__(self, *args, **kwargs):
        """Subclassing pandas dataframe. CommTable.from_[csv,etc]()
        Can be used for reading in a table.
        """
        _table.__init__(self, *args, **kwargs)


    def taxonInLib(self, taxon_name, libID):
        """Checking whether taxon in the selected library.
        Parameters
        ----------
        taxon_name : str
            taxon name
        libID : str
            library ID

        Return
        ------
        b : boolean
        """
        libID = re.sub('\D+', '', libID)
        dfSub = self.df.loc[self.df['taxon_name'] == taxon_name]
        libs = [str(x) for x in dfSub['library'].tolist()]
        return libID in libs

            
    def get_taxonAbund(self, taxon_name, libID=None, abs_abund=False):
        """Getting the abundance(s) of taxon in the comm file.

        Parameters
        ----------
        taxon_name : str
            name of taxon
        libID : str
            library ID. If None, all libraries selected.
        abs_abund : bool
             return absolute abundance instead of relative abundance

        Returns
        -------
        abund_vals : iterable 
            abundance values for the taxon
        """
        retCol = 'abs_abund' if abs_abund else 'rel_abund_perc'
        assert retCol in self.df.columns, \
            '"{}" column not found'.format(retCol)
        if libID is not None:            
            df_sub =  self.df.loc[(self.df['library'] == libID) &
                                  (self.df['taxon_name'] == taxon_name)]
        else:
            df_sub =  self.df.loc[(self.df['taxon_name'] == taxon_name)]

        return df_sub[retCol].tolist()

    
    def get_unique_libIDs(self):
        """Get all unique libIDs from the community table."""
        return self.df['library'].unique().tolist()

    def get_unique_taxon_names(self):
        """Get all unique taxon names from the community table."""
        return self.df['taxon_name'].unique()        

    def get_total_abund(self, libID=None, abs_abund=False):
        """Get total abundance of each library or all libraries"""
        retCol = 'abs_abund' if abs_abund else 'rel_abund_perc'
        if libID is not None:
            return self.df.loc[self.df['library'] == libID][retCol].sum()
        else:
            return self.df[retCol].sum()
        
    @property
    def abs_abund(self):
        """The absolute abundance of each taxon based on user-value.

        Parameters
        ----------
        abs_abund : int
            total abundance of all taxa in each library
        """
        return self.df['abs_abund']

    @abs_abund.setter
    def abs_abund(self, abs_abund):
        """rel_abund_perc stored as a percent"""
        try:
            x = self.df['rel_abund_perc']
            self.df['abs_abund'] = x / 100 * int(abs_abund)
        except KeyError:
            raise KeyError('"rel_abund_perc" column not found in comm file')

        self.df['abs_abund'] = np.round(self.df['abs_abund'].tolist(), 0)
        self.df['abs_abund'] = self.df['abs_abund'].astype(int)
