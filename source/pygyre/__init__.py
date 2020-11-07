# Module  : pygyre top level module
# Purpose : Python support for GYRE
#
# Copyright 2020 Rich Townsend & The GYRE Team
#
# This file is part of PYGYRE. PYGYRE is free software: you can
# redistribute it and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation, version 3.
#
# PYGYRE is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import h5py as h5
import re
import astropy.table as tb

def read_output(filename, filetype='auto'):

    """
    Read GYRE output data from a file
    
    Parameters
    -----------
    
    filename : string
        Input file name

    filetype : string (default 'auto')
        Input file type (one of 'HDF', 'TXT' or 'auto')
    
    Returns
    -------

    table: astropy.table.Table
        Table containing file data
    """

    # If necessary, determine the filetype

    if filetype == 'auto':

        if h5.is_hdf5(filename):
            filetype = 'HDF'
        else:
            filetype = 'TXT'

    # Create the table

    if filetype == 'HDF':
        tab = _read_output_hdf(filename)
    elif filetype == 'TXT':
        tab = _read_output_txt(filename)
    else:
        raise Exception("Invalid filetype '{:s}'".format(filetype))

    return tab

##

def _read_output_hdf(filename):

    # Create the table

    tab = _read_generic_hdf(filename)

    # Convert complex columns and attributes

    complex_dtype = np.dtype([('re', '<f8'), ('im', '<f8')])

    for colname in tab.colnames:

        if tab[colname].dtype == complex_dtype:
            tab[colname] = tab[colname]['re'] + 1j*tab[colname]['im']

    for atrname in tab.meta:

        if tab.meta[atrname].dtype == complex_dtype:
            tab.meta[atrname] = tab.meta[atrname]['re'] + 1j*tab.meta[atrname]['im']

    return tab

##

def _read_output_txt(filename):

    # Create the table

    tab = _read_generic_txt(filename)

     # Merge real/imag columns

    for colname in tab.colnames:

        mo = re.fullmatch(r'Re\((.*)\)', colname)

        if mo is not None:

            colname = mo.group(1)

            colname_re = 'Re({:s})'.format(colname)
            colname_im = 'Im({:s})'.format(colname)
            
            tab[colname] = tab[colname_re] + 1j*tab[colname_im]

            tab.remove_columns([colname_re, colname_im])

    return tab

##

def read_model(filename, filetype='auto'):

    """
    Read stellar model data from a file
    
    Parameters
    -----------
    
    filename : string
        Input file name

    filetype : string (default 'auto')
        Input file type (one of 'MESA', 'GSM', 'POLY')
    
    Returns
    -------

    table: astropy.table.Table
        Table containing file data
    """

    # If necessary, determine the filetype

    if filetype == 'auto':

        if h5.is_hdf5(filename):

            with h5.File(filename, 'r') as file:
                if 'M_star' in file.attrs:
                    filetype = 'GSM'
                elif 'n_poly' in file.attrs:
                    filetype = 'POLY'

        else:

            filetype = 'MESA'

    # Create the table

    if filetype == 'GSM':
        tab = _read_model_gsm(filename)
    elif filetype == 'MESA':
        tab = _read_model_mesa(filename)
    elif filetype == 'POLY':
        tab = _read_model_poly(filename)
    else:
        raise Exception("Invalid filetype '{:s}'".format(filetype))

    return tab

##

def _read_model_gsm(filename):

    # Create the table

    tab = _read_model_gsm(filename)

    return tab

##

def _read_model_mesa(filename):

    # Create the table

    tab = _read_generic_txt(filename)

    return tab

##

def _read_model_poly(filename):

    # Create the table

    tab = read_poly(filename)

    # Add in structure coefficients

    z = tab['z']
    theta = tab['theta']
    dtheta = tab['dtheta']
    n_poly = tab['n_poly']
    B = tab['B']
    t = tab['t']
    mu = tab['mu']

    z_s = tab.meta['z_s']
    mu_s = tab.meta['mu_s']
    Gamma_1 = tab.meta['Gamma_1']

    tab['x'] = z/z_s

    with np.nditer([z, theta, dtheta, n_poly, t, mu, None, None, None, None]) as it:

        for z_, theta_, dtheta_, n_poly_, t_, mu_, V_2, As, U, c_1 in it:

            if z_ != 0:
                if theta_ != 0:
                    V_2[...] =  -(n_poly_ + 1.)*z_s**2*dtheta_/(theta_*z_)
                    As[...] = V_2*(z_/z_s)**2*(n_poly_/(n_poly_+1.) - 1./Gamma_1)
                else:
                    V_2[...] = np.Inf
                    As[...] = np.Inf
                U[...] = z_**3*t_*theta_**n_poly_/mu_
                c_1[...] = (z_/z_s)**3/(mu_/mu_s)
            else:
                V_2[...] = (n_poly_ + 1.)*z_s**2/3.
                As[...] = 0.
                U[...] = 3.
                c_1[...] = 3.*mu_s/z_s**3

        tab['V_2'] = it.operands[-4]
        tab['As'] = it.operands[-3]
        tab['U'] = it.operands[-2]
        tab['c_1'] = it.operands[-1]

    # Add in physical data

    tab['rho/rho_c'] = t*theta**n_poly
    tab['P/P_c'] = 1/B*(n_poly[0]+1.)/(n_poly+1.)*t**(1.-1./n_poly)*tab['rho/rho_c']**(1.+1./n_poly)
    tab['m/M'] = mu/mu_s

    return tab

##

def read_poly(filename):

    # Read data from the file

    with h5.File(filename, 'r') as file:

        meta = dict(zip(file.attrs.keys(),file.attrs.values()))

        cols = {}
    
        for key in file.keys():
            cols[key] = file[key][...]

    # Evaluate auxillary data

    z = cols['z']
    theta = cols['theta']
    dtheta = cols['dtheta']

    n_z = len(z)
    n_r = meta['n_r']

    n_poly = meta['n_poly']

    if n_r <= 1:
        meta['Delta_b'] = np.empty([0])
    Delta_b = meta['Delta_b']

    t = np.empty(n_r)
    B = np.empty(n_r)
    z_b = np.empty(n_r-1)
    v_i = np.empty(n_r)
    v_f = np.empty(n_r)
    mu_i = np.empty(n_r)
    mu_f = np.empty(n_r)

    t[0] = 1.
    B[0] = 1.
    v_i[0] = 0.
    mu_i[0] = 0.

    n_z = len(cols['z'])

    n_poly_z = np.empty([n_z])
    t_z = np.empty([n_z])
    B_z = np.empty([n_z])
    mu_z = np.empty([n_z])

    n_poly_z[0] = n_poly[0]
    t_z[0] = t[0]
    B_z[0] = B[0]
    mu_z[0] = 0.

    i = 0

    for k in range(1, n_z):

        if z[k] == z[k-1]:
            z_b[i] = z[k-1]
            v_f[i] = z[k-1]**2*dtheta[k-1]
            mu_f[i] = mu_i[i] - t[i]/B[i]*(v_f[i] - v_i[i])
            i += 1
            t[i] = t[i-1]*np.exp(n_poly[i-1]*np.log(theta[k-1]) - n_poly[i]*np.log(theta[k]) + Delta_b[i-1])
            B[i] = B[i-1]*(n_poly[i-1]+1)/(n_poly[i]+1)*theta[k]**(n_poly[i]+1)/theta[k-1]**(n_poly[i-1]+1)*t[i]**2/t[i-1]**2
            v_i[i] = z[k]**2*dtheta[k]
            mu_i[i] = mu_f[i-1]
 
        n_poly_z[k] = n_poly[i]
        t_z[k] = t[i]
        B_z[k] = B[i]
        mu_z[k] = mu_i[i] - t[i]/B[i]*(z[k]**2*dtheta[k] - v_i[i])

    v_f[-1] = z[-1]**2*dtheta[-1]
    mu_f[-1] = mu_i[-1] - t[-1]/B[-1]*(v_f[-1] - v_i[-1])

    # Add auxillary data to cols & meta

    cols['B'] = B_z
    cols['t'] = t_z
    cols['mu'] = mu_z
    cols['n_poly'] = n_poly_z

    meta['t'] = t
    meta['B'] = B
    meta['n_z'] = n_z
    meta['z_b'] = z_b
    meta['v_i'] = v_i
    meta['v_f'] = v_f
    meta['mu_i'] = mu_i
    meta['mu_f'] = mu_f
    meta['mu_s'] = mu_z[-1]
    meta['z_s'] = z[-1]

    # Create the table

    tab = tb.Table(cols, meta=meta)

    return tab

##

def _read_generic_hdf(filename):

    # Read data from the file

    with h5.File(filename, 'r') as file:

        meta = dict(zip(file.attrs.keys(),file.attrs.values()))

        cols = {}
    
        for key in file.keys():
            cols[key] = file[key][...]

    # Create the table

    tab = tb.Table(cols, meta=meta)

    return tab

##

def _read_generic_txt(filename):

    # Read data from the file

    with open(filename, 'r') as file:

        file.readline()

        meta_keys = file.readline().split()
        meta_vals = [_interpret(val) for val in file.readline().split()]

        meta = dict(zip(meta_keys, meta_vals))

        file.readline()
        file.readline()

        col_keys = file.readline().split()
        col_vals = []
    
        while True :
            line = file.readline()
            if not line : break
            col_vals.append([_interpret(val) for val in line.split()])

        file.close()

    cols = {}

    for i in range(0,len(col_keys)):
        cols[col_keys[i]] = np.array([col_vals[i] for col_vals in col_vals])

    # Create the table

    return tb.Table(cols, meta=meta)

##

def _interpret(s):

    try:

        return int(s)

    except ValueError:

        sr = re.sub(r'D', 'E', s)
        sr = re.sub(r'(\d)([+-])', r'\1E\2', sr)

        try:

            return float(sr)

        except ValueError:

            return s
