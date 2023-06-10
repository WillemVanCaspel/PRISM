#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 11:25:58 2020

@author: wim
"""

import numpy as np
from scipy.io import FortranFile
import pyshtools as sh
from wrf import interplevel
from netCDF4 import Dataset
from pyhwm2014 import HWM142D
from nrlmsise00 import msise_flat
from scipy import interpolate
from datetime import datetime

#%

latfac = -1 # -1 switches from 90-(-90) to (-90)-90 for shtools output

csp_1 = 1 # Condon-Shortley phase factor(s)
csp_2 = 1

"""
##############################################################################
# SIGPRIM FORTRAN-PYTHON I/O                                                 #
##############################################################################

"""

def gauleg(ymin,
           ymax,
           ydim):
    """
    Numerical recipe for cos(phi) term used in Nabla operator.
    
    INPUT:  - a: -1 (90 deg South),
            - b: 1 (90 deg North),
            - c: y-dimension of the grid (functions as K2 in SIGPRIM)
    
    OUTPUT: - mu-term in cos(phi) calculation
    
    """
    
    coefs, w = np.zeros(ydim), np.zeros(ydim)

    m=(ydim + 1) / 2
    xm = 0.5 * (ymax + ymin)
    xl = 0.5 * (ymax - ymin)
    
    for i in range(1,int(m)):
        z = np.cos(np.pi * (i - 0.25) / (ydim + 0.5))    
        p1=1.
        p2=0.
        
        for j in range(1,int(ydim)):
            p3 = p2
            p2 = p1
            p1 = ((2. * j - 1.) * z * p2-(j - 1.) * p3) / j
        
        pp = ydim * (z * p1 - p2) / (z * z - 1.)
        z1 = z
        z = z1 - p1 / pp
        coefs[i - 1] = xm - xl * z
        coefs[int(ydim) - i] = xm + xl * z
        w[i - 1] = 2. * xl / ((1. - z * z) * pp * pp)
        w[int(ydim) - i] = w[i - 1]
    
    return coefs


def curlz(u,
          v,
          dphi,
          dlam):
    """
    Calculate z-component of the curl as defined in Saravanan p. 28.
    
    - INPUT:  - Even order square [dim, dim] (DH) or [dim, 2* dim] (DH2) grid
                of (u,v) winds
              - dphi (polar) and dlam (azimuthal) grid spacing in RADIANS
              
    - OUTPUT: - Even complex spherical harmonics [N,M] of the curlz of (u,v)
                [in spherical coordinates]              
                
    """
    
    a = 6.37e6 # radius earth [m]
    
    ydim, xdim = np.shape(u)
        
    mu = gauleg(1.,-1.,ydim)
    nui = 1 / np.sqrt(1 - mu**2) # 1 / cos(phi)
    nu = np.sqrt(1 - mu**2)      # cos(phi)
    
    # multiply V_lam (u) by cos(phi)
    for k in range(xdim):
        u[:,k] = nu * u[:,k]
        
    # gradient terms
    dv_dlam = np.gradient(v,dlam,axis=1)
    du_dphi = np.gradient(u,dphi,axis=0)
    
    for k in range(xdim):
        dv_dlam[:,k] = nui / a * dv_dlam[:,k]
        du_dphi[:,k] = nui / a * du_dphi[:,k]
        
    total_field = dv_dlam - du_dphi
    
    # spherical harmonics expansion
    sh_curl = sh.expand.SHExpandDHC(total_field)
    even_sh = sh_curl[0,:,:]
    
    return even_sh


def div(u,
        v,
        dphi,
        dlam):
    
    """
    Calculate div as defined in Saravanan p. 28.
    
    - INPUT:  - Even order square [dim, dim] (DH) or [dim, 2* dim] (DH2) grid
                of (u,v) winds
              - dphi (polar) and dlam (azimuthal) grid spacing in RADIANS
              
    - OUTPUT: - Even complex spherical harmonics [N,M] of the curlz of (u,v)
                [in spherical coordinates]              
                
    """
    
    a = 6.37e6 # radius earth [m]
    
    ydim, xdim = np.shape(u)
        
    mu = gauleg(1.,-1.,ydim)
    nui = 1 / np.sqrt(1 - mu**2) # 1 / cos(phi)
    nu = np.sqrt(1 - mu**2)      # cos(phi)
    
    # multiply V_phi (u) by cos(phi)
    for k in range(xdim):
        v[:,k] = nu * v[:,k]
        
    # gradient terms
    du_dlam = np.gradient(u,dlam,axis=1)
    dv_dphi = np.gradient(v,dphi,axis=0)
    
    for k in range(xdim):
        du_dlam[:,k] = nui / a * du_dlam[:,k]
        dv_dphi[:,k] = nui / a * dv_dphi[:,k]
        
    total_field = du_dlam + dv_dphi
    
    # spherical harmonics expansion
    sh_div = sh.expand.SHExpandDHC(total_field)
    even_sh = sh_div[0,:,:]
    
    return even_sh


def uv_cross_product(v_11,
                     v_12,
                     lat_grid,
                     lon_grid):
    """ 
    Calculates k x (nabla_H phi), as Saravanan p. 17.

    Orientation sigprim input vectors: - v_11 = i unit vec (zonal)
                                       - v_12 = j unit vec (meridional)
                                       
    Output: (u,v) = A x B in spherical coordinates, with A = (0*i, 0*j, 1*k)
                               and B the intput vector (v_11*i, v_12*j, 0*k)
                               
    Function transforms input vectors from spherical coordinates to cartesian
    coordinates, calculates the cross-products, and then transforms back to 
    spherical coordinates.
    
    """
      
    ydim, xdim = np.shape(v_11)
    
    # initialize phase grid and convert to radians
    phase_phi = np.zeros((ydim,xdim))
    phase_lam = np.zeros((ydim,xdim))

    for k in range(xdim):
        phase_phi[:,k] = lat_grid # 90 to -90 deg
    for k in range(ydim):
        phase_lam[k,:] = lon_grid # 0 to 360 deg
        
    phase_phi = phase_phi / 180 * np.pi
    phase_lam = phase_lam / 180 * np.pi

    # phi from [pi / 2, -pi / 2] to [0, pi] relative to positive z-axis
    phase_phi = -(phase_phi - np.pi / 2) # 0-2pi or 2pi-0 makes no diff here

    sin_phi = np.sin(phase_phi) 
    cos_phi = np.cos(phase_phi) 
    sin_lam = np.sin(phase_lam)
    cos_lam = np.cos(phase_lam)
    
    # vectors should be in [k,j,i] order, k-vector is zero here
    sv_11 = np.zeros((ydim,xdim))
    sv_12 = v_12
    sv_13 = v_11
    
    # i and j 2nd input vector are zero, k = 1 
    sv_21 = np.ones((ydim,xdim))
    sv_22 = np.zeros((ydim,xdim))
    sv_23 = np.zeros((ydim,xdim))
    
    # spherical to cartesian transformation of k-vector
    a_1 = sin_phi*cos_lam*sv_21 + cos_phi*cos_lam*sv_22 - sin_lam*sv_23
    a_2 = sin_phi*sin_lam*sv_21 + cos_phi*sin_lam*sv_22 + cos_lam*sv_23
    a_3 = cos_phi*sv_21 - sin_phi*sv_22 + 0*sv_23
    
    # Nabla_psi vector
    b_1 = sin_phi*cos_lam*sv_11 + cos_phi*cos_lam*sv_12 - sin_lam*sv_13
    b_2 = sin_phi*sin_lam*sv_11 + cos_phi*sin_lam*sv_12 + cos_lam*sv_13
    b_3 = cos_phi*sv_11 - sin_phi*sv_12 + 0*sv_13
         
    # definition cross product for two vectors (x,y,z) in cartesian coords
    s1 = a_2 * b_3 - a_3 * b_2
    s2 = a_3 * b_1 - a_1 * b_3
    s3 = a_1 * b_2 - a_2 * b_1
    
    # transform back to spherical coordinates [k,j,i]
    sv_31 = sin_phi*cos_lam*s1 + sin_phi*sin_lam*s2 + cos_phi*s3
    sv_32 = cos_phi*cos_lam*s1 + cos_phi*sin_lam*s2 - sin_phi*s3
    sv_33 = -sin_lam*s1 + cos_lam*s2 + 0*s3
        
    # return i, j vectors for zonal and meridional winds along surface
    u = sv_33
    v = sv_32
    
    return u, v


def read_uv(file,
            folder,
            rec_n,
            grid_dim,
            levels=None,
            interp=False,
            grid_type='DH'):
    """
    Construct horizontal wind field from divergence and vorticity spectral
    coefficients.
    
    INPUT: - file: file to read in
           - folder: path to file to read in
           - record number (corresponding to some model time)
           - dimension of output grid (e.g. 64 for a lat-lon grid of 64 by 64)
           - levels: pressure levels on which output is interpolated; optional
           - grid_type: 'DH'  for [grid_dim, grid_dim] output
                        'DH2' for [grid_dim, 2 * grid_dim] output size
                        
    OUTPUT: - zonal & meridional wind fields on [alt,lat,lon] grid

    """

    a = 6.37e6 # radius earth [m]
    
    f = FortranFile(folder + file,'r')
    
    # read in grid dimensions / metadata gunk
    md = f.read_record(dtype=np.int32)
    (k1, k2, m, n, l) = (md[0], md[1], md[2], md[3], md[4])
    
    slv = f.read_record(dtype=np.float32)
    lat = f.read_record(dtype=np.float32)
    zstdlvl = f.read_record(dtype=np.float32)
    
#    np.save('/home/wim/Desktop/Projects/Model/vertical_grids/sigprim_lvls_ERAz.npy',slv)
       
    nmax_grid = int(grid_dim / 2) # max n + 1, controls output grid resolution
    
    if grid_type == 'DH':
        grid_multiplier = 1
    else:
        grid_multiplier = 2

    # y, x, z dim
    u_arr = np.zeros((nmax_grid * 2,nmax_grid * 2 * grid_multiplier, l))
    v_arr = np.zeros((nmax_grid * 2,nmax_grid * 2 * grid_multiplier, l))
    
    for k in range(rec_n):
        # read in spherical harmonics
        jvor = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
        jdiv = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
        jpot = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
        pssp = f.read_record(dtype=np.complex64).reshape(m+1,n+1,order='F')
        geob = f.read_record(dtype=np.float32).reshape(k1,k2,order='F')
    f.close()
    
    # surface pressure field
    even = pssp[:,:]      # p_{m,n} coefficients
    odd = np.copy(even)   # p_{-m,n}

    for m_i in range(m + 1):
            odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
        
    coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                      dtype=np.complex64)

    coeffs[0,:,:] = even
    coeffs[1,1:,:] = odd[1:,:]
    
    # change m,n arangement to that used by shtools
    ct = np.swapaxes(coeffs,1,2)
    
    # zero pad to square dim. nmax_grid >= m + 1
    shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
    shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
    x = sh.SHCoeffs.from_array(shs,normalization='4pi',csphase=csp_2)
    
    grid = x.expand(grid=grid_type,extend=False)
    ps = np.exp(np.real(grid.to_array())) # exponent due to coordinate system
    p_field = np.zeros((nmax_grid * 2,nmax_grid * 2 * grid_multiplier, l))
    
    for level in range(l):
        # pressure levels of (u,v)-field calculated from model sigma-lvls
        p_field[:,:,level] = ps * slv[level]
        
    # calculate u,v field from vorticity and divergence; see Saravanan manual
    psi = np.copy(jvor)
    xi = np.copy(jdiv)
    
    # Del operator x-y coeffs from numerical recipe gauleg(ymin,ymax,ydim)
    mu = gauleg(-1.,1.,nmax_grid * 2)
    nui = 1 / np.sqrt(1 - mu**2) 
        
    for lvl in range(l):
        for n_sh in range(1,n + 1):
            # n_sh runs from 1 to (n + 1); jvor(n=0) and jdiv(n=0) are (0,0j)
            eigenvalue_leg = -n_sh * (n_sh + 1) / (a ** 2)  # Saravanan p. 17
            psi[:,n_sh,lvl] = 1 / eigenvalue_leg * jvor[:,n_sh,lvl]
            xi[:,n_sh,lvl] = 1 / eigenvalue_leg * jdiv[:,n_sh,lvl]
       
    # vorticity term
    for lvl_i in range(l):
        even = psi[:,:,lvl_i] # p_{m,n} coefficients
        odd = np.copy(even)   # p_{-m,n}
        
        for m_i in range(m + 1):
            odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
        
        coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                          dtype=np.complex64)
    
        coeffs[0,:,:] = even
        coeffs[1,1:,:] = odd[1:,:]
        
        # change m,n arangement to that used by shtools
        ct = np.swapaxes(coeffs,1,2)
        
        # zero pad to square dim. nmax_grid >= m + 1
        shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
        shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
        x = sh.SHCoeffs.from_array(shs,normalization='4pi',csphase=csp_2)
        
        grid = x.expand(grid=grid_type,extend=False)
        data = np.real(grid.to_array())
        lats = grid.lats()
        lons = grid.lons()

        dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
        dlam = (lons[1] - lons[0]) / 180 * np.pi

        dim_x = np.shape(data)[1]
                   
        grad_y, grad_x = np.gradient(data,dphi,dlam)
        x_coeffs = np.zeros(np.shape(grad_x))
        y_coeffs = np.zeros(np.shape(grad_y))
            
        for j in range(dim_x):
            x_coeffs[:,j] = nui / a # nui is symmetric 
            y_coeffs[:,j] = 1 / a   # dave uses np.sqrt(1 - mu**2) / a here instead?
        
        u, v = uv_cross_product(x_coeffs * grad_x,
                                y_coeffs * grad_y,
                                lats,
                                lons)
    
        u_arr[:,:,lvl_i] = u
        v_arr[:,:,lvl_i] = v
        
        # if lvl_i == 60:
        #     plt.figure()
        #     plt.title('u vor')
        #     plt.pcolormesh(u)
        #     plt.colorbar()
            
        #     plt.figure()
        #     plt.title('v vor')
        #     plt.pcolormesh(v)
        #     plt.colorbar()
        
    # divergence term
    for lvl_i in range(l):
        even = xi[:,:,lvl_i] # p_{m,n} coefficients
        odd = np.copy(even)   # p_{-m,n}
        
        for m_i in range(m + 1):
            odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
        
        coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                          dtype=np.complex64)
    
        coeffs[0,:,:] = even
        coeffs[1,1:,:] = odd[1:,:]
        
        # change m,n arangement to that used by shtools
        ct = np.swapaxes(coeffs,1,2)
        
        # zero pad to square dim. nmax_grid >= m + 1
        shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
        shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
        x = sh.SHCoeffs.from_array(shs,normalization='4pi',csphase=csp_2)
        
        grid = x.expand(grid=grid_type,extend=False)
        data = np.real(grid.to_array())
        lats = grid.lats()
        lons = grid.lons()
            
        dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
        dlam = (lons[1] - lons[0]) / 180 * np.pi
        
        dim_x = np.shape(data)[1]
                   
        grad_y, grad_x = np.gradient(data,dphi,dlam)
        x_coeffs = np.zeros(np.shape(grad_y))
        y_coeffs = np.zeros(np.shape(grad_y))
            
        for j in range(dim_x):
            x_coeffs[:,j] = nui / a
            y_coeffs[:,j] = 1 / a # Dave uses nu here instead, manual 1 / a
        
        u = x_coeffs * (grad_x)
        v = y_coeffs * grad_y
    
        u_arr[:,:,lvl_i] = u_arr[:,:,lvl_i] + u 
        v_arr[:,:,lvl_i] = v_arr[:,:,lvl_i] + v
        
        # if lvl_i == 60:
        #     plt.figure()
        #     plt.title('u div')
        #     plt.pcolormesh(u)
        #     plt.colorbar()
            
        #     plt.figure()
        #     plt.title('v div')
        #     plt.pcolormesh(v)
        #     plt.colorbar()
            

    # set dimensions to python configuration of [zdim,ydim,xdim] from y, u, z
    u_arr = np.swapaxes(u_arr,0,2) # z, u, y
    u_arr = np.swapaxes(u_arr,1,2) # z, y, u
    v_arr = np.swapaxes(v_arr,0,2)
    v_arr = np.swapaxes(v_arr,1,2)
    p_field = np.swapaxes(p_field,0,2)
    p_field = np.swapaxes(p_field,1,2)
    
    if interp:
        # interpolate to input pressure levels (Pa)
        u_arr = interplevel(u_arr,p_field,levels)
        v_arr = interplevel(v_arr,p_field,levels)
            
    return u_arr, v_arr, p_field, ps


def read_uv_band(file,
                 folder,
                 day_begin,
                 day_final,
                 grid_dim,
                 levels=None,
                 interp=False,
                 grid_type='DH'):
    """
    Construct horizontal wind field from divergence and vorticity spectral
    coefficients.
    
    INPUT: - record number (corresponding to some model time)
           - dimension of output grid (e.g. 64 by 64 lat-lon)
           - levels: levels (km currently) onto which output is interpolated
           - grid_type: 'DH' for [grid_dim, grid_dim] output
                        'DH2' for [grid_dim, 2 * grid_dim] output size
                        
    OUTPUT: - zonal wind field (u) on [alt,lat,lon] grid
            - meridional wind field (v) on [alt,lat,lon] grid
    """
    
    full_u = np.zeros((day_final * 8,np.shape(levels)[0],grid_dim,grid_dim))
    full_v = np.zeros((day_final * 8,np.shape(levels)[0],grid_dim,grid_dim))

    a = 6.37e6 # radius earth [m]
    
    f = FortranFile(folder + file,'r')
    
    # read in grid dimensions / metadata
    md = f.read_record(dtype=np.int32)
    (k1, k2, m, n, l) = (md[0], md[1], md[2], md[3], md[4])
    
    slv = f.read_record(dtype=np.float32)
    lat = f.read_record(dtype=np.float32)
    zstdlvl = f.read_record(dtype=np.float32)
       
    nmax_grid = int(grid_dim / 2) # max n + 1, controls output grid resolution
    
    if grid_type == 'DH':
        grid_multiplier = 1
    else:
        grid_multiplier = 2
    
    for t_steps in range(day_final * 8):
        if t_steps == 0:
            # y, x, z dim
            u_arr = np.zeros((nmax_grid * 2,nmax_grid * 2 * grid_multiplier, l))
            v_arr = np.zeros((nmax_grid * 2,nmax_grid * 2 * grid_multiplier, l))
            
            for k in range(day_begin * 8 + 1):
                # read in spherical harmonics
                jvor = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                jdiv = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                jpot = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                pssp = f.read_record(dtype=np.complex64).reshape(m+1,n+1,order='F')
                geob = f.read_record(dtype=np.float32).reshape(k1,k2,order='F')
            
            # surface pressure field
            even = pssp[:,:]      # p_{m,n} coefficients
            odd = np.copy(even)   # p_{-m,n}
        
            for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
            coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                              dtype=np.complex64)
        
            coeffs[0,:,:] = even
            coeffs[1,1:,:] = odd[1:,:]
            
            # change m,n arangement to that used by shtools
            ct = np.swapaxes(coeffs,1,2)
            
            # zero pad to square dim. nmax_grid >= m + 1
            shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
            shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
            x = sh.SHCoeffs.from_array(shs,normalization='4pi',csphase=csp_2)
            
            grid = x.expand(grid=grid_type)
            ps = np.exp(np.real(grid.to_array())) # exponent due to coordinate system
            p_field = np.zeros((nmax_grid * 2,nmax_grid * 2 * grid_multiplier, l))
            
            for level in range(l):
                # pressure coordinates of (u,v)
                p_field[:,:,level] = ps * slv[level]
                
            # u,v field from vorticity and divergence
            psi = np.copy(jvor)
            xi = np.copy(jdiv)
            
            # Del operator x-y coeffs from numerical recipe gauleg(ymin,ymax,ydim)
            mu = gauleg(-1.,1.,nmax_grid * 2)
            nui = 1 / np.sqrt(1 - mu**2) 
                
            for lvl in range(l):
                for n_sh in range(1,n + 1):
                    # n_sh runs from 1 to (n + 1); jvor(n=0) and jdiv(n=0) are (0,0j)
                    eigenvalue_leg = -n_sh * (n_sh + 1) / (a ** 2)  # Saravanan p. 17
                    psi[:,n_sh,lvl] = 1 / eigenvalue_leg * jvor[:,n_sh,lvl]
                    xi[:,n_sh,lvl] = 1 / eigenvalue_leg * jdiv[:,n_sh,lvl]
                            
            # vorticity term
            for lvl_i in range(l):
                even = psi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                x = sh.SHCoeffs.from_array(shs,normalization='4pi',csphase=csp_2)
                
                grid = x.expand(grid=grid_type)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
        
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
        
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_x))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    x_coeffs[:,j] = nui / a # nui is symmetric 
                    y_coeffs[:,j] = 1 / a   # dave uses np.sqrt(1 - mu**2) / a here instead?
                
                u, v = uv_cross_product(x_coeffs * grad_x,
                                        y_coeffs * grad_y,
                                        lats,
                                        lons)
            
                u_arr[:,:,lvl_i] = u
                v_arr[:,:,lvl_i] = v
                
            # divergence term
            for lvl_i in range(l):
                even = xi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                x = sh.SHCoeffs.from_array(shs,normalization='4pi',csphase=csp_2)
                
                grid = x.expand(grid=grid_type)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
                    
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
                
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_y))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    x_coeffs[:,j] = nui / a
                    y_coeffs[:,j] = 1 / a # Dave uses nu here instead, manual 1 / a
                
                u = x_coeffs * (grad_x)
                v = y_coeffs * grad_y
            
                u_arr[:,:,lvl_i] = u_arr[:,:,lvl_i] + u 
                v_arr[:,:,lvl_i] = v_arr[:,:,lvl_i] + v
        
            # set dimensions to python configuration of [zdim,ydim,xdim] from y, u, z
            u_arr = np.swapaxes(u_arr,0,2) # z, u, y
            u_arr = np.swapaxes(u_arr,1,2) # z, y, u
            v_arr = np.swapaxes(v_arr,0,2)
            v_arr = np.swapaxes(v_arr,1,2)
            p_field = np.swapaxes(p_field,0,2)
            p_field = np.swapaxes(p_field,1,2)
            
            if interp:
                # interpolate to input pressure levels (Pa)
                u_arr = interplevel(u_arr,p_field,levels)
                v_arr = interplevel(v_arr,p_field,levels)
                
            print('Reading record: ', t_steps)

            full_u[t_steps,:,:,:] = u_arr
            full_v[t_steps,:,:,:] = v_arr
            
        else:
            # y, x, z dim
            u_arr = np.zeros((nmax_grid * 2,nmax_grid * 2 * grid_multiplier, l))
            v_arr = np.zeros((nmax_grid * 2,nmax_grid * 2 * grid_multiplier, l))
            
            # read in spherical harmonics
            jvor = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            jdiv = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            jpot = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            pssp = f.read_record(dtype=np.complex64).reshape(m+1,n+1,order='F')
            geob = f.read_record(dtype=np.float32).reshape(k1,k2,order='F')
            
            # surface pressure field
            even = pssp[:,:]      # p_{m,n} coefficients
            odd = np.copy(even)   # p_{-m,n}
        
            for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
            coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                              dtype=np.complex64)
        
            coeffs[0,:,:] = even
            coeffs[1,1:,:] = odd[1:,:]
            
            # change m,n arangement to that used by shtools
            ct = np.swapaxes(coeffs,1,2)
            
            # zero pad to square dim. nmax_grid >= m + 1
            shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
            shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
            x = sh.SHCoeffs.from_array(shs,normalization='4pi',csphase=csp_2)
            
            grid = x.expand(grid=grid_type)
            ps = np.exp(np.real(grid.to_array())) # exponent due to coordinate system
            p_field = np.zeros((nmax_grid * 2,nmax_grid * 2 * grid_multiplier, l))
            
            for level in range(l):
                # pressure coordinates of (u,v)
                p_field[:,:,level] = ps * slv[level]
                
            # u,v field from vorticity and divergence
            psi = np.copy(jvor)
            xi = np.copy(jdiv)
            
            # Del operator x-y coeffs from numerical recipe gauleg(ymin,ymax,ydim)
            mu = gauleg(-1.,1.,nmax_grid * 2)
            nui = 1 / np.sqrt(1 - mu**2) 
                
            for lvl in range(l):
                for n_sh in range(1,n + 1):
                    # n_sh runs from 1 to (n + 1); jvor(n=0) and jdiv(n=0) are (0,0j)
                    eigenvalue_leg = -n_sh * (n_sh + 1) / (a ** 2)  # Saravanan p. 17
                    psi[:,n_sh,lvl] = 1 / eigenvalue_leg * jvor[:,n_sh,lvl]
                    xi[:,n_sh,lvl] = 1 / eigenvalue_leg * jdiv[:,n_sh,lvl]
                            
            # vorticity term
            for lvl_i in range(l):
                even = psi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                x = sh.SHCoeffs.from_array(shs,normalization='4pi',csphase=csp_2)
                
                grid = x.expand(grid=grid_type)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
        
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
        
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_x))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    x_coeffs[:,j] = nui / a # nui is symmetric 
                    y_coeffs[:,j] = 1 / a   # Dave uses np.sqrt(1 - mu**2) / a here instead?
                
                u, v = uv_cross_product(x_coeffs * grad_x,
                                        y_coeffs * grad_y,
                                        lats,
                                        lons)
            
                u_arr[:,:,lvl_i] = u
                v_arr[:,:,lvl_i] = v
                
            # divergence term
            for lvl_i in range(l):
                even = xi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                x = sh.SHCoeffs.from_array(shs,normalization='4pi',csphase=csp_2)
                
                grid = x.expand(grid=grid_type)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
                    
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
                
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_y))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    x_coeffs[:,j] = nui / a
                    y_coeffs[:,j] = 1 / a # Dave uses nu here instead, manual 1 / a
                
                u = x_coeffs * (grad_x)
                v = y_coeffs * grad_y
            
                u_arr[:,:,lvl_i] = u_arr[:,:,lvl_i] + u 
                v_arr[:,:,lvl_i] = v_arr[:,:,lvl_i] + v
        
            # set dimensions to python configuration of [zdim,ydim,xdim] from y, u, z
            u_arr = np.swapaxes(u_arr,0,2) # z, u, y
            u_arr = np.swapaxes(u_arr,1,2) # z, y, u
            v_arr = np.swapaxes(v_arr,0,2)
            v_arr = np.swapaxes(v_arr,1,2)
            p_field = np.swapaxes(p_field,0,2)
            p_field = np.swapaxes(p_field,1,2)
            
            if interp:
                # interpolate to input pressure levels (Pa)
                u_arr = interplevel(u_arr,p_field,levels)
                v_arr = interplevel(v_arr,p_field,levels)
            
            print('Reading record: ', t_steps)
            
            full_u[t_steps,:,:,:] = u_arr
            full_v[t_steps,:,:,:] = v_arr
            
    f.close()
            
    return full_u, full_v, p_field


def read_temperature_plvl(file,
                          folder,
                          rec_n,
                          grid_dim,
                          levels=None,
                          interp=False,
                          grid_type='DH'):
    """
    Construct horizontal temperature field on pressure lvls from 
    JPOT spherical harmonic coefficients.
    
    INPUT: - record number (corresponding to some model time)
           - dimension of output grid (e.g. 64 by 64 lat-lon)
           - grid_type: 'DH' for [grid_dim, grid_dim] output
                        'DH2' for [grid_dim, 2 * grid_dim] output size
                        
    OUTPUT: - temperature field (K) on pressure levels
    """
    
    f = FortranFile(folder + file,'r')
    
    # read in grid dimensions / metadata..irrelevant errors here, m[4] etc. INT
    md = f.read_record(dtype=np.int32)
    (k1, k2, m, n, l) = (md[0], md[1], md[2], md[3], md[4])
    
    slv = f.read_record(dtype=np.float32)
    lat = f.read_record(dtype=np.float32)
    zstdlvl = f.read_record(dtype=np.float32)
       
    nmax_grid = int(grid_dim / 2) # max n + 1, controls output grid resolution
    
    if grid_type == 'DH':
        grid_multiplier = 1
    elif grid_type == 'DH2':
        grid_multiplier = 2

    t_arr = np.zeros((nmax_grid * 2,nmax_grid * 2 * grid_multiplier, l))
    
    for k in range(rec_n):
        # read in spherical harmonics
        jvor = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
        jdiv = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
        jpot = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
        pssp = f.read_record(dtype=np.complex64).reshape(m+1,n+1,order='F')
        geob = f.read_record(dtype=np.float32).reshape(k1,k2,order='F')
    f.close()
    
    # surface pressure field
    even = pssp[:,:] # p_{m,n} coefficients
    odd = np.copy(even)   # p_{-m,n}
    for m_i in range(m + 1):
            odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
        
    coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                      dtype=np.complex64)

    coeffs[0,:,:] = even
    coeffs[1,1:,:] = odd[1:,:]
    
    # change m,n arangement to that used by shtools
    ct = np.swapaxes(coeffs,1,2)
    
    # zero pad to square dim. nmax_grid >= m + 1
    shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
    shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
    x = sh.SHCoeffs.from_array(shs,normalization='4pi',csphase=csp_2)
    
    grid = x.expand(grid=grid_type,extend=False)
    ps = np.exp(np.real(grid.to_array())) # exponent due to coordinate system
    p_field = np.zeros((nmax_grid * 2,nmax_grid * 2 * grid_multiplier, l))
    
    for level in range(l):
        # pressure coordinates of (u,v)
        p_field[:,:,level] = ps * slv[level]
        
    temp = np.copy(jpot) # temperature spherical harmonics    
                
    # temperature to grid routine
    for lvl_i in range(l):
        even = temp[:,:,lvl_i] # p_{m,n} coefficients
        odd = np.copy(even)    # p_{-m,n}
        
        for m_i in range(m + 1):
            odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
        
        coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                          dtype=np.complex64)
    
        coeffs[0,:,:] = even
        coeffs[1,1:,:] = odd[1:,:]
        
        # change m,n arangement to that used by shtools
        ct = np.swapaxes(coeffs,1,2)
        
        # zero pad to square dim. nmax_grid >= m + 1
        shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
        shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
        x = sh.SHCoeffs.from_array(shs,normalization='4pi',csphase=csp_2)
        
        grid = x.expand(grid=grid_type,extend=False)
        data = np.real(grid.to_array())
    
        t_arr[:,:,lvl_i] = data

    # set dimensions to the usual python array: [zdim,ydim,xdim]
    t_arr = np.swapaxes(t_arr,0,2)
    t_arr = np.swapaxes(t_arr,1,2)
    p_field = np.swapaxes(p_field,0,2)
    p_field = np.swapaxes(p_field,1,2)
    
    if interp:
        t_arr = interplevel(t_arr,p_field,levels)
        
    return t_arr


def zonal_mean_sigprim_to_python(folder,
                                 file_name):
    
    """
    Files used by this function must be initialized using
    zonal_mean_fortran_python.f
    
    OUTPUT:
        
       - z_ax: Geometric height coordinate in km
       - y_ax: Latitude coordine in degrees
       - t_ax: Time coordinate in day of year
       - t_field: zonal mean temperature field
       - u_field: zonal mean zonal wind field
       
    """
            
    # read in dimensions of the grid
    f = FortranFile(folder + file_name + '_grid_dim.dat','r')
    
    a = f.read_reals(dtype=np.float32)
    ydim, zdim, tdim = int(a[0]),int(a[1]),int(a[2])
    
    f.close()
    
    # read in grid axis
    f = FortranFile(folder + file_name + '_grid.dat','r')
    
    a = f.read_record(dtype=np.float32)
    y_ax = a[:ydim]
    z_ax = a[ydim:(ydim + zdim)]
    t_ax = a[(ydim + zdim):(ydim + zdim + tdim)]
    
    f.close()
    
    # read in zonal mean wind and temperature fields
    f = FortranFile(folder + file_name + '_fields.dat','r')
    
    u_field = np.empty((tdim,zdim,ydim))
    t_field = np.empty((tdim,zdim,ydim))
    
    for t in range(tdim):
        t_field[t,:,:] = f.read_record(dtype=np.float32).reshape(ydim, zdim, 
                                                                 order='F').T
        u_field[t,:,:] = f.read_record(dtype=np.float32).reshape(ydim, zdim, 
                                                                 order='F').T
    
    f.close()

    return t_field, u_field, y_ax, z_ax, t_ax


def uvt_to_sigprim(u_in_f,
                   v_in_f,
                   t_in_f,
                   dphi,
                   dlam,
                   z_levels,
                   t_steps,
                   wave_numbers,
                   folder,
                   output_name):
    
    """
    ALTITUDE CAN BE OMITTED HERE AND IN FORTRAN. CHANGE?
    
    Creates a file containing spherical harmonics of the horizontal wind
    and temperature fields which can be assimilated by sigprim.
    
    Input files (arrays) on equidistant [lat,lon] grid with spacing [dphi,dlam]
    in radians. Preferably on an even square (DH) or [dim, 2 * dim] (DH2) grid. 
    If not, arrays will have to be interpolated to DH2 first.
    
    Array should have dimensions [time,altitude,latitude,longitude], where
        - time 0 to end
        - altitude high from low (that is, low Pa to high Pa)
        - latitude from -90 to 90
        - longitude from 0 to 360
        
    array_size, determined by max_m and max_n, should preferably
    be kept fixed (and be sufficiently large to cover all future needs?).
    
    - for the winds: calculate curl and transform to spectral space. 
    
    """
        
    max_m = np.float32(45)
    max_n = np.float32(45)
    max_l = np.float32(np.shape(u_in_f)[1]) # <= L1MAX (201) in SIGPRIM
    
    u_field = u_in_f
    v_field = v_in_f
    t_field = t_in_f
    
    tdim, zdim, ydim, xdim = np.shape(u_field)
        
    # write meta-data and z,t grid
    f = FortranFile(folder + output_name,'w')
    
    zlvls = np.zeros(int(max_l),dtype=np.float32)
    zlvls[:zdim] = z_levels
    
    # write grid dims as integers. *2 factor to make rec len match complex len
    rec_int = np.zeros(int(max_m * max_n * max_l * 2),dtype=np.int32)
    rec_int[0], rec_int[1], rec_int[2], rec_int[3] = tdim, max_m, max_n, max_l
    f.write_record(rec_int)
    
    # write time and altitude array as floats
    rec_fl = np.zeros(int(max_m * max_n * max_l * 2),dtype=np.float32)
    rec_fl[:int(max_l)] = np.array(zlvls)
        
    rec_fl[int(max_l):(int(max_l) + tdim)] = np.array(t_steps)
    f.write_record(rec_fl)
        
    # write spherical harmonics of the curl of (u,v) at each level, time step
    for t in range(tdim):
        curl_coeffs = np.zeros((int(max_n),int(max_m),int(max_l)),
                               dtype=np.complex64)
        temp_coeffs = np.zeros((int(max_n),int(max_m),int(max_l)),
                               dtype=np.complex64)
        
        for l in range(zdim):
            # curlz returns even shs of curlz field with coeff dim [N,M].
            curlz_coeff = curlz(u_field[t,l,:,:],
                                v_field[t,l,:,:],
                                dphi,
                                dlam)
                        
            # shs of temperature field, expand assumes lat[0] = 90, flip order
            temp_coeff = sh.expand.SHExpandDHC(t_field[t,l,:,:],csphase=csp_1)
            temp_coeff = temp_coeff[0,:,:] # even shs
            
            # N should not exceed max_n used by the model. Future error msg.
            dmax = np.shape(curlz_coeff)[0]
            
            # minus curlz for the time tendency, see Saravanan p. 10
            curl_coeffs[:dmax,:dmax,l] = -1 * curlz_coeff
            temp_coeffs[:dmax,:dmax,l] = temp_coeff
        
        curl_temp = np.zeros(np.shape(curl_coeffs),dtype=np.complex64)
        temp_temp = np.zeros(np.shape(curl_coeffs),dtype=np.complex64)
        
        # select zonal wave number fields to be assimilated
        for m in wave_numbers:
            curl_temp[:,int(m),:] = curl_coeffs[:,int(m),:]
            temp_temp[:,int(m),:] = temp_coeffs[:,int(m),:]
            
        curl_coeffs = curl_temp
        temp_coeffs = temp_temp
        
        # reshape to [M,N] sh arangement used by sigprim
        curl_coeffs = np.swapaxes(curl_coeffs,0,1)
        temp_coeffs = np.swapaxes(temp_coeffs,0,1)
                
        # write record (beginning now from rec_n = 3)
        print('writing rec (ND): ', t)
        f.write_record(curl_coeffs.reshape(int(max_m * max_n * max_l),
                                           order='F'))
        f.write_record(temp_coeffs.reshape(int(max_m * max_n * max_l),
                                           order='F'))
        
    f.close()
    
    return None


def uvt_to_sigprim_div(u_in_f,
                       v_in_f,
                       t_in_f,
                       dphi,
                       dlam,
                       z_levels,
                       t_steps,
                       wave_numbers,
                       folder,
                       output_name):
    
    """
    ALTITUDE CAN BE OMITTED HERE AND IN FORTRAN. CHANGE?
    
    Creates a file containing spherical harmonics of the horizontal wind
    and temperature fields which can be assimilated by sigprim.
    
    Input files (arrays) on equidistant [lat,lon] grid with spacing [dphi,dlam]
    in radians. Preferably on an even square (DH) or [dim, 2 * dim] (DH2) grid. 
    If not, arrays will have to be interpolated to DH2 first.
    
    Array should have dimensions [time,altitude,latitude,longitude], where
        - time 0 to end
        - altitude high from low (that is, low Pa to high Pa)
        - latitude from -90 to 90
        - longitude from 0 to 360
        
    array_size, determined by max_m and max_n, should preferably
    be kept fixed (and be sufficiently large to cover all future needs?).
    
    - for the winds: calculate curl and transform to spectral space. 
    
    """
        
    max_m = np.float32(45)
    max_n = np.float32(45)
    max_l = np.float32(np.shape(u_in_f)[1]) # <= L1MAX (201) in SIGPRIM
    
    u_field = u_in_f
    v_field = v_in_f
    t_field = t_in_f
    
    tdim, zdim, ydim, xdim = np.shape(u_field)
        
    # write meta-data and z,t grid
    f = FortranFile(folder + output_name,'w')
    
    zlvls = np.zeros(int(max_l),dtype=np.float32)
    zlvls[:zdim] = z_levels
    
    # write grid dims as integers. *2 factor to make rec len match complex len
    rec_int = np.zeros(int(max_m * max_n * max_l * 2),dtype=np.int32)
    rec_int[0], rec_int[1], rec_int[2], rec_int[3] = tdim, max_m, max_n, max_l
    f.write_record(rec_int)
    
    # write time and altitude array as floats
    rec_fl = np.zeros(int(max_m * max_n * max_l * 2),dtype=np.float32)
    rec_fl[:int(max_l)] = np.array(zlvls)
        
    rec_fl[int(max_l):(int(max_l) + tdim)] = np.array(t_steps)
    f.write_record(rec_fl)
        
    # write spherical harmonics of the curl of (u,v) at each level, time step
    for t in range(tdim):
        curl_coeffs = np.zeros((int(max_n),int(max_m),int(max_l)),
                               dtype=np.complex64)
        dive_coeffs = np.zeros((int(max_n),int(max_m),int(max_l)),
                               dtype=np.complex64)
        temp_coeffs = np.zeros((int(max_n),int(max_m),int(max_l)),
                               dtype=np.complex64)
        
        for l in range(zdim):
            # curlz returns even shs of curlz field with coeff dim [N,M].
            curlz_coeff = curlz(u_field[t,l,:,:],
                                v_field[t,l,:,:],
                                dphi,
                                dlam)
            
            dive_coeff = div(u_field[t,l,:,:],
                             v_field[t,l,:,:],
                             dphi,
                             dlam)
                        
            # shs of temperature field, expand assumes lat[0] = 90, flip order
            temp_coeff = sh.expand.SHExpandDHC(t_field[t,l,:,:],csphase=csp_1)
            temp_coeff = temp_coeff[0,:,:] # even shs
            
            # N should not exceed max_n used by the model. Future error msg.
            dmax = np.shape(curlz_coeff)[0]
            
            # minus curlz for the time tendency, see Saravanan p. 10
            curl_coeffs[:dmax,:dmax,l] = -1 * curlz_coeff
            dive_coeffs[:dmax,:dmax,l] = -1 * dive_coeff
            temp_coeffs[:dmax,:dmax,l] = temp_coeff
        
        curl_temp = np.zeros(np.shape(curl_coeffs),dtype=np.complex64)
        dive_temp = np.zeros(np.shape(curl_coeffs),dtype=np.complex64)
        temp_temp = np.zeros(np.shape(curl_coeffs),dtype=np.complex64)
        
        # select zonal wave number fields to be assimilated
        for m in wave_numbers:
            curl_temp[:,int(m),:] = curl_coeffs[:,int(m),:]
            dive_temp[:,int(m),:] = dive_coeffs[:,int(m),:]
            temp_temp[:,int(m),:] = temp_coeffs[:,int(m),:]
            
        curl_coeffs = curl_temp
        dive_coeffs = dive_temp
        temp_coeffs = temp_temp
        
        # reshape to [M,N] sh arangement used by sigprim
        curl_coeffs = np.swapaxes(curl_coeffs,0,1)
        dive_coeffs = np.swapaxes(dive_coeffs,0,1)
        temp_coeffs = np.swapaxes(temp_coeffs,0,1)
                
        # write record (beginning now from rec_n = 3)
        print('writing rec (ND): ', t)
        f.write_record(curl_coeffs.reshape(int(max_m * max_n * max_l),
                                           order='F'))
        f.write_record(dive_coeffs.reshape(int(max_m * max_n * max_l),
                                           order='F'))
        f.write_record(temp_coeffs.reshape(int(max_m * max_n * max_l),
                                           order='F'))
        
    f.close()
    
    return None


def uvt_to_sigprim_zm_input(u_in_f,
                            v_in_f,
                            t_in_f,
                            dphi,
                            dlam,
                            z_levels,
                            t_steps,
                            folder,
                            output_name):
    
    """    
    Creates a file containing spherical harmonics of zonal mean horizontal wind
    and temperature fields that can be assimilated by sigprim.
    
    Input files (arrays) on equidistant [lat,lon] grid with spacing [dphi,dlam]
    in radians on an even square (DH) or [dim, 2 * dim] (DH2) grid. 
    
    Array should have dimensions [time,altitude,latitude,longitude], where
        - time from time start to time end (e.g., DOY -31 for padded files)
        - altitude in pressure coordinates from low pressure to high pressure
          corresponding to sigprim model levels
        - latitude from south to north
        - longitude from 0 to 360
        
    So what happens to meridional winds here?
    """
        
    max_m = np.float32(1)
    max_n = np.float32(45)
    max_l = np.float32(np.shape(u_in_f)[1]) 
       
    tdim, zdim, ydim = np.shape(u_in_f)
    
    full_u = np.empty((zdim,ydim,ydim))
    full_v = np.empty((zdim,ydim,ydim))
    full_t = np.empty((zdim,ydim,ydim))
        
    # write meta-data and z,t grid
    f = FortranFile(folder + output_name,'w')
    
    zlvls = np.zeros(int(max_l), dtype=np.float32)
    zlvls[:zdim] = z_levels
    
    # write grid dims as integers. Factor 2 to make rec len match complex len
    rec_int = np.zeros(int(max_m * max_n * max_l * 2),dtype=np.int32)
    rec_int[0], rec_int[1], rec_int[2], rec_int[3] = tdim, max_m, max_n, max_l
    f.write_record(rec_int)
    
    # write time and altitude array as floats
    rec_fl = np.zeros(int(max_m * max_n * max_l * 2),dtype=np.float32)
    rec_fl[:int(max_l)] = np.array(zlvls)
    rec_fl[int(max_l):(int(max_l) + tdim)] = np.array(t_steps)
    f.write_record(rec_fl)
        
    # write spherical harmonics of the curl of (u,v) at each level, time step
    for t in range(tdim):
        curl_coeffs = np.zeros((int(max_n),int(max_m),int(max_l)),
                               dtype=np.complex64)
        temp_coeffs = np.zeros((int(max_n),int(max_m),int(max_l)),
                               dtype=np.complex64)
        
        # populate square dim lat-lon grid to calculate spherical harmonics
        for y in range(ydim):
            full_u[:,:,y] = u_in_f[t,:,:]
            full_v[:,:,y] = v_in_f[t,:,:]
            full_t[:,:,y] = t_in_f[t,:,:]
                            
        for l in range(zdim):
            # curlz returns even shs of curlz field with coeff dim [N,M].
            curlz_coeff = curlz(full_u[l,:,:],
                                full_v[l,:,:],
                                dphi,
                                dlam)
                        
            # shs of temperature field, expand assumes lat[0] = 90, flip order
            temp_coeff = sh.expand.SHExpandDHC(full_t[l,:,:])
            temp_coeff = temp_coeff[0,:,:] # even shs
            
            # sh.expand always square output, i.e. M=N
            dmax = np.shape(curlz_coeff)[0]
            mmax = np.minimum(int(dmax),int(max_m),dtype='int')
            nmax = np.minimum(int(dmax),int(max_n),dtype='int')
                        
            # minus curlz for the time tendency (?), see Saravanan p. 10
            curl_coeffs[:nmax,:mmax,l] = -1 * curlz_coeff[:nmax,:mmax]
            temp_coeffs[:nmax,:mmax,l] = temp_coeff[:nmax,:mmax]
        
        # reshape from [N,M] to [M,N] sh arangement used by sigprim
        curl_coeffs = np.swapaxes(curl_coeffs,0,1)
        temp_coeffs = np.swapaxes(temp_coeffs,0,1)
                
        # write record (beginning now from rec_n = 3)
        print('writing rec (ND): ', t)
        f.write_record(curl_coeffs.reshape(int(max_m * max_n * max_l),
                                           order='F'))
        f.write_record(temp_coeffs.reshape(int(max_m * max_n * max_l),
                                           order='F'))
        
    f.close()
    
    return None


def uvt_to_sigprim_zm_div(u_in_f,
                          v_in_f,
                          t_in_f,
                          dphi,
                          dlam,
                          z_levels,
                          t_steps,
                          folder,
                          output_name):
    
    """    
    Creates a file containing spherical harmonics of zonal mean horizontal wind
    and temperature fields that can be assimilated by sigprim.
    
    Input files (arrays) on equidistant [lat,lon] grid with spacing [dphi,dlam]
    in radians on an even square (DH) or [dim, 2 * dim] (DH2) grid. 
    
    Array should have dimensions [time,altitude,latitude,longitude], where
        - time from time start to time end (e.g., DOY -31 for padded files)
        - altitude in pressure coordinates from low pressure to high pressure
          corresponding to sigprim model levels
        - latitude from south to north
        - longitude from 0 to 360
        
    """
        
    max_m = np.float32(1)  # only zonal mean fields
    max_n = np.float32(45)
    max_l = np.float32(81) # number of vertical levels used by sigprim
    
    u_field = u_in_f
    v_field = v_in_f
    t_field = t_in_f
    
    tdim, zdim, ydim, xdim = np.shape(u_field)
        
    # write meta-data and z,t grid
    f = FortranFile(folder + output_name,'w')
    
    zlvls = np.zeros(int(max_l), dtype=np.float32)
    zlvls[:zdim] = z_levels
    
    # write grid dims as integers. Factor 2 to make rec len match complex len
    rec_int = np.zeros(int(max_m * max_n * max_l * 2),dtype=np.int32)
    rec_int[0], rec_int[1], rec_int[2], rec_int[3] = tdim, max_m, max_n, max_l
    f.write_record(rec_int)
    
    # write time and altitude array as floats
    rec_fl = np.zeros(int(max_m * max_n * max_l * 2),dtype=np.float32)
    rec_fl[:int(max_l)] = np.array(zlvls)
    rec_fl[int(max_l):(int(max_l) + tdim)] = np.array(t_steps)
    f.write_record(rec_fl)
        
    # write spherical harmonics of the curl of (u,v) at each level, time step
    for t in range(tdim):
        curl_coeffs = np.zeros((int(max_n),int(max_m),int(max_l)),
                               dtype=np.complex64)
        div_coeffs = np.zeros((int(max_n),int(max_m),int(max_l)),
                               dtype=np.complex64)
        temp_coeffs = np.zeros((int(max_n),int(max_m),int(max_l)),
                               dtype=np.complex64)
                
        for l in range(zdim):
            # curlz returns even shs of curlz field with coeff dim [N,M].
            curlz_coeff = curlz(u_field[t,l,:,:],
                                v_field[t,l,:,:],
                                dphi,
                                dlam)
            
            # div returns even shs of div field with coeff dim [N,M].
            div_coeff = div(u_field[t,l,:,:],
                            v_field[t,l,:,:],
                            dphi,
                            dlam)
                        
            # shs of temperature field, expand assumes lat[0] = 90, flip order
            temp_coeff = sh.expand.SHExpandDHC(t_field[t,l,:,:],csphase=csp_1)
            temp_coeff = temp_coeff[0,:,:] # even shs
            
            # sh.expand always square output, i.e. M=N
            dmax = np.shape(curlz_coeff)[0]
            mmax = np.minimum(int(dmax),int(max_m),dtype='int')
            nmax = np.minimum(int(dmax),int(max_n),dtype='int')
                        
            # see Saravanan p. 10
            curl_coeffs[:nmax,:mmax,l] = curlz_coeff[:nmax,:mmax]
            div_coeffs[:nmax,:mmax,l] = div_coeff[:nmax,:mmax]
            temp_coeffs[:nmax,:mmax,l] = temp_coeff[:nmax,:mmax]
        
        # reshape to [M,N] sh arangement used by sigprim
        curl_coeffs = np.swapaxes(curl_coeffs,0,1)
        div_coeffs = np.swapaxes(div_coeffs,0,1)
        temp_coeffs = np.swapaxes(temp_coeffs,0,1)
                
        # write record (beginning now from rec_n = 3)
        print('writing rec (ND): ', t)
        f.write_record(curl_coeffs.reshape(int(max_m * max_n * max_l),
                                           order='F'))
        f.write_record(div_coeffs.reshape(int(max_m * max_n * max_l),
                                           order='F'))
        f.write_record(temp_coeffs.reshape(int(max_m * max_n * max_l),
                                           order='F'))
        
    f.close()
    
    return None


def uv_read_lat(file_n,
                folder,
                output_dim,
                latitude,
                day_start,
                n_days,
                p_lvls):

    """
    
    Latitude between (-90, 90]
    
    """
        
    u_values = np.zeros((n_days * 8,5,output_dim))
    v_values = np.zeros((n_days * 8,5,output_dim))
    
    model_lvls = sigprim_plvls()
    model_lvls_3d =  np.zeros((np.shape(model_lvls)[0],
                               output_dim,
                               output_dim))
    
    # populate 3d model level pressure field
    for x in range(output_dim):
        for y in range(output_dim):
            model_lvls_3d[:,y,x,] = model_lvls
            
    out_lats = -np.linspace(-90,90,output_dim,endpoint=False)[::-1] + 90
    latitude = latitude + 90
    lat_ind = np.argmin(np.abs(out_lats - latitude))      

    for t in range(n_days * 8):
        u, v, p = read_uv(file_n,
                          folder,
                          day_start * 8 + t + 1,
                          output_dim,
                          levels=p_lvls,
                          interp=True,
                          grid_type='DH')
            
        # u_temp = interplevel(u,model_lvls_3d,p_lvls)
        # v_temp = interplevel(v,model_lvls_3d,p_lvls)
        u_values[t,:,:] = u[:,lat_ind,:]
        v_values[t,:,:] = v[:,lat_ind,:]
        
        print('Reading timestep: ', t)
        
    return u_values, v_values


def pressure_to_altitude(p,
                         t):
    
    # zero meters at bottom boundary
    z = np.zeros(np.shape(p))
            
    g_0 = 9.80665   # m s**-2
    r_d = 8.3144598 # J mol**-1 K**-1
    r_e = 6.3781e6  # m
    m = 0.0289644   # kg mol**-1
    
    for k in np.arange(80,0,-1):
        # height adjusted gravitational acceleration assuming spherical earth
        g = g_0 * (r_e / (r_e + z[k,:,:]))**2
        
        # forward difference at bottom boundary
        if k == 80:
            dz = -(p[k-1,:,:] - p[k,:,:]) / (p[k,:,:] * m * g) * r_d * t[k,:,:]
            
        # central difference otherwise
        else:
            dz = (-0.5 * ((p[k-1,:,:] - p[k+1,:,:])) / 
                  (p[k,:,:] * m * g) * r_d * t[k,:,:])
        
        # increment altitude with dz
        z[k-1,:,:] = z[k,:,:] + dz

    return z


def pressure_to_altitude_mm(p,
                            t):
    
    """
    Numerically integrate barometric equation.
    
    - takes into account gravity and molar mass as function of altitude
    - molar mass profile based on global mean from MSISE-00 
    
    inputs: 
    
        p : pressure (Pa)
        t : temperature (K)
    
    """
    
    # zero meters at bottom boundary
    z = np.zeros(np.shape(p))
            
    g_0 = 9.80665   # m s**-2
    r_d = 8.3144598 # J mol**-1 K**-1
    r_e = 6.3781e6  # m
    
    path_f = '/home/wim/Desktop/Projects/Model/Data/'
    molar_mass = np.load(path_f + 'molar_mass.npy')
    mm_heights = np.load(path_f + 'molar_mass_height_m.npy')
    
    l_lvls = np.shape(p)[0]
            
    for k in np.arange(l_lvls - 1,0,-1):
        # height adjusted gravitational acceleration assuming spherical earth
        g = g_0 * (r_e / (r_e + z[k,:,:]))**2
        
        # global mean molar mass estimate
        mean_height = np.mean(z[k,:,:])
        
        func_interp = interpolate.interp1d(mm_heights,molar_mass)
        m = func_interp(mean_height)
        
        # forward difference at bottom boundary
        if k == (l_lvls - 1):
            dz = -(p[k-1,:,:] - p[k,:,:]) / (p[k,:,:] * m * g) * r_d * t[k,:,:]
            
        # central difference otherwise
        else:
            dz = (-0.5 * ((p[k-1,:,:] - p[k+1,:,:])) / 
                  (p[k,:,:] * m * g) * r_d * t[k,:,:])
        
        # increment altitude with dz
        z[k-1,:,:] = z[k,:,:] + dz
        
        # temp[k-1] = np.mean(z[k-1,:,:])
        
    # np.save('/home/wim/Desktop/' + 'sigprim_zlvls_m.npy',temp)
    
    return z

def pressure_to_altitude_topography(p,
                                    t,
                                    geopotential,
                                    model_latitudes):
    
    """
    Numerically integrate barometric equation.
    
    - takes into account gravity and molar mass as function of altitude
    - molar mass profile based on global mean from MSISE-00 
    
    inputs: 
    
        p : pressure (Pa)
        t : temperature (K)
    
    """
    
    # array to hold geometric altitude values
    z = np.zeros(np.shape(p))
            
    g_0 = 9.80665   # m s**-2
    r_d = 8.3144598 # J mol**-1 K**-1
    r_e = 6.3781e6  # m
    
    path_f = '/home/wim/Desktop/Projects/Model/Data/'
    molar_mass = np.load(path_f + 'molar_mass.npy')
    mm_heights = np.load(path_f + 'molar_mass_height_m.npy')
      
    mean_f = geopotential / g_0 # geopotential to altitude (m)
    
    # k1 latitude dim (gaussian), k2 longitude (equidistant)
    k1_surf, k2_surf = np.shape(mean_f)
    
    grid_lon = np.linspace(0,360,k2_surf,endpoint=False)
    grid_lat = model_latitudes # gaussian sigprim grid
            
    ff = interpolate.interp2d(grid_lon,
                              grid_lat,
                              mean_f,
                              kind='linear')
    
    dz_dim, dy_dim, dx_dim = np.shape(p)
    
    sig_lon = np.linspace(0,360,dx_dim,endpoint=False)
    sig_lat = np.linspace(-90,90,dy_dim,endpoint=False)
    
    z[-1,:,:] = ff(sig_lon,
                   sig_lat)
    
    l_lvls = np.shape(p)[0]
            
    for k in np.arange(l_lvls - 1,0,-1):
        # height adjusted gravitational acceleration assuming spherical earth
        g = g_0 * (r_e / (r_e + z[k,:,:]))**2
        
        # molar mass from global mean height-profile
        mean_height = np.mean(z[k,:,:])
        func_interp = interpolate.interp1d(mm_heights,molar_mass)

        if mean_height < 0:
            m = func_interp(0)
        else:
            m = func_interp(mean_height)
            
        # forward difference at bottom boundary
        if k == (l_lvls - 1):
            dz = -(p[k-1,:,:] - p[k,:,:]) / (p[k,:,:] * m * g) * r_d * t[k,:,:]
            
        # central difference otherwise
        else:
            dz = (-0.5 * ((p[k-1,:,:] - p[k+1,:,:])) / 
                  (p[k,:,:] * m * g) * r_d * t[k,:,:])
        
        # increment altitude with dz
        z[k-1,:,:] = z[k,:,:] + dz
            
    return z

# def hydro_ideal(p,
#                 t,
#                 geopotential,
#                 model_latitudes):
    
#     z = np.zeros(np.shape(p))
    
#     r = 287.058     # specific gas constant [J kg-1 K-1]
#     g_0 = 9.80665   # gravitational acceleration [m s**-2]
#     r_e = 6.3781e6  # radius earth [m]
    
#     mean_f = geopotential / g_0 # geopotential to altitude (m)
    
#     # k1 latitude dim (gaussian), k2 longitude (equidistant)
#     k1_surf, k2_surf = np.shape(mean_f)
    
#     grid_lon = np.linspace(0,360,k2_surf,endpoint=False)
#     grid_lat = model_latitudes # gaussian sigprim grid
            
#     ff = interpolate.interp2d(grid_lon,
#                               grid_lat,
#                               mean_f,
#                               kind='linear')
    
#     dz_dim, dy_dim, dx_dim = np.shape(p)
    
#     sig_lon = np.linspace(0,360,dx_dim,endpoint=False)
#     sig_lat = np.linspace(-90,90,dy_dim,endpoint=False)
    
#     z[-1,:,:] = ff(sig_lon,
#                     sig_lat)
    
#     l_lvls = np.shape(p)[0]
    
#     for k in np.arange(l_lvls - 1,0,-1):
#         # height adjusted gravitational acceleration assuming spherical earth
#         g = g_0 * (r_e / (r_e + z[k,:,:]))**2
        
#         # forward difference at bottom boundary
#         if k == (l_lvls - 1):
#             dz = -(p[k-1,:,:] - p[k,:,:]) / p[k,:,:] * r * t[k,:,:] / g
            
#             # dz = -(p[k-1,:,:] - p[k,:,:]) / (p[k,:,:] * m * g) * r_d * t[k,:,:]
            
#         # central difference otherwise
#         else:
#             # dz = -0.5 * ((p[k-1,:,:] - p[k+1,:,:])) / p[k,:,:] * r * t[k,:,:] / g
#             dz = -(p[k-1,:,:] - p[k,:,:]) / p[k,:,:] * r * t[k,:,:] / g
                  
#             # dz = (-0.5 * ((p[k-1,:,:] - p[k+1,:,:])) / 
#             #       (p[k,:,:] * m * g) * r_d * t[k,:,:])
        
#         # increment altitude with dz
#         z[k-1,:,:] = z[k,:,:] + dz
        
#     return z

def hydro_ideal_init(p,
                     t):
    
    z = np.zeros(np.shape(p))
    
    r = 287.058     # specific gas constant [J kg-1 K-1]
    g_0 = 9.80665   # gravitational acceleration [m s**-2]
    r_e = 6.3781e6  # radius earth [m]
                    
    l_lvls = np.shape(p)[0]
    ln_p = np.log(p)
    
    for k in np.arange(l_lvls - 1,0,-1):
        # height adjusted gravitational acceleration assuming spherical earth
        g = g_0 * (r_e / (r_e + z[k]))**2
        
        # forward difference at bottom boundary
        if k == (l_lvls - 1):
            dz = -(ln_p[k-1] - ln_p[k]) * r * t[k] / g
            
            # dz = -(p[k-1,:,:] - p[k,:,:]) / (p[k,:,:] * m * g) * r_d * t[k,:,:]
            
        # central difference otherwise
        else:
            # dz = -0.5 * ((p[k-1,:,:] - p[k+1,:,:])) / p[k,:,:] * r * t[k,:,:] / g
            dz = -(ln_p[k-1] - ln_p[k]) * r * t[k] / g
                  
            # dz = (-0.5 * ((p[k-1,:,:] - p[k+1,:,:])) / 
            #       (p[k,:,:] * m * g) * r_d * t[k,:,:])
        
        # increment altitude with dz
        z[k-1] = z[k] + dz
        
    return z


def hydro_ideal(p,
                t,
                geopotential,
                model_latitudes):
    
    z = np.zeros(np.shape(p))
    
    r = 287.058     # specific gas constant [J kg-1 K-1]
    g_0 = 9.80665   # gravitational acceleration [m s**-2]
    r_e = 6.3781e6  # radius earth [m]
    
    mean_f = geopotential / g_0 # geopotential to altitude (m)
    
    # k1 latitude dim (gaussian), k2 longitude (equidistant)
    k1_surf, k2_surf = np.shape(mean_f)
    
    grid_lon = np.linspace(0,360,k2_surf,endpoint=False)
    grid_lat = model_latitudes # gaussian sigprim grid
            
    ff = interpolate.interp2d(grid_lon,
                              grid_lat,
                              mean_f,
                              kind='linear')
    
    dz_dim, dy_dim, dx_dim = np.shape(p)
    
    sig_lon = np.linspace(0,360,dx_dim,endpoint=False)
    sig_lat = np.linspace(-90,90,dy_dim,endpoint=False)
    
    z[-1,:,:] = ff(sig_lon,
                    sig_lat)
    
    l_lvls = np.shape(p)[0]
    ln_p = np.log(p)
    
    for k in np.arange(l_lvls - 1,0,-1):
        # height adjusted gravitational acceleration assuming spherical earth
        g = g_0 * (r_e / (r_e + z[k,:,:]))**2
        
        # forward difference at bottom boundary
        if k == (l_lvls - 1):
            dz = -(ln_p[k-1,:,:] - ln_p[k,:,:]) * r * t[k,:,:] / g
            
            # dz = -(p[k-1,:,:] - p[k,:,:]) / (p[k,:,:] * m * g) * r_d * t[k,:,:]
            
        # central difference otherwise
        else:
            # dz = -0.5 * ((p[k-1,:,:] - p[k+1,:,:])) / p[k,:,:] * r * t[k,:,:] / g
            dz = -(ln_p[k-1,:,:] - ln_p[k,:,:]) * r * t[k,:,:] / g
                  
            # dz = (-0.5 * ((p[k-1,:,:] - p[k+1,:,:])) / 
            #       (p[k,:,:] * m * g) * r_d * t[k,:,:])
        
        # increment altitude with dz
        z[k-1,:,:] = z[k,:,:] + dz
        
    return z
        
def pressure_to_altitude_scaleheight(p,
                                     scale_height,
                                     geopotential,
                                     model_latitudes):
    
    """
    Numerically integrate barometric equation.
    
    - takes into account gravity and molar mass as function of altitude
    - molar mass profile based on global mean from MSISE-00 
    
    inputs: 
    
        p : pressure (Pa)
        t : temperature (K)
    
    """
    
    # array to hold geometric altitude values
    z = np.zeros(np.shape(p))
            
    g_0 = 9.80665   # m s**-2
    r_d = 8.3144598 # J mol**-1 K**-1
    r_e = 6.3781e6  # m
    
    path_f = '/home/wim/Desktop/Projects/Model/Data/'
    molar_mass = np.load(path_f + 'molar_mass.npy')
    mm_heights = np.load(path_f + 'molar_mass_height_m.npy')
      
    mean_f = geopotential / g_0 # geopotential to altitude (m)
    
    # k1 latitude dim (gaussian), k2 longitude (equidistant)
    k1_surf, k2_surf = np.shape(mean_f)
    
    grid_lon = np.linspace(0,360,k2_surf,endpoint=False)
    grid_lat = model_latitudes # gaussian sigprim grid
            
    ff = interpolate.interp2d(grid_lon,
                              grid_lat,
                              mean_f,
                              kind='linear')
    
    dz_dim, dy_dim, dx_dim = np.shape(p)
    
    sig_lon = np.linspace(0,360,dx_dim,endpoint=False)
    sig_lat = np.linspace(-90,90,dy_dim,endpoint=False)
    
    z[-1,:,:] = ff(sig_lon,
                   sig_lat)
    
    l_lvls = np.shape(p)[0]
            
    for k in np.arange(l_lvls - 1,0,-1):
        
        # increment altitude with dz
        z[k-1,:,:] = z[-1,:,:] - scale_height * np.log(p[k-1,:,:]/p[-1,:,:])
            
    return z
        

def uv_read_altitude(file,
                     folder,
                     day_begin,
                     day_final,
                     grid_dim,
                     interp_levels=[100000]):
    
    """
    
    Uses barometric formula to convert model lvls to altitude
    
    """
        
    full_u = np.zeros((day_final * 8,np.shape(interp_levels)[0],grid_dim,grid_dim))
    full_v = np.zeros((day_final * 8,np.shape(interp_levels)[0],grid_dim,grid_dim))

    a = 6.37e6 # radius earth [m]
    
    f = FortranFile(folder + file,'r')
    
    # read in grid dimensions 
    md = f.read_record(dtype=np.int32)
    (k1, k2, m, n, l) = (md[0], md[1], md[2], md[3], md[4])
    
    # model levels in sigma coordinates
    slv = f.read_record(dtype=np.float32)
    
    # skip next two unused entries
    f.read_record(dtype=np.float32)
    f.read_record(dtype=np.float32)
       
    nmax_grid = int(grid_dim / 2) # max n + 1, controls output grid resolution
        
    for t_steps in range(day_final * 8):
        if t_steps == 0:
            # y, x, z dim arrays to hold read in fields
            u_arr = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
            v_arr = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
            t_arr = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
            
            for k in range(day_begin * 8 + 1):
                # read in spherical harmonic coefficients (shc)
                jvor = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                jdiv = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                jpot = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                pssp = f.read_record(dtype=np.complex64).reshape(m+1,n+1,order='F')
                f.read_record(dtype=np.float32).reshape(k1,k2,order='F') 
            
            # read in surface pressure field and construct odd shcs
            even = pssp[:,:]      
            odd = np.copy(even)   
            
            for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
            coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                              dtype=np.complex64)
            
            # change [m,n] arangement to [n,m] used by shtools
            coeffs[0,:,:] = even
            coeffs[1,1:,:] = odd[1:,:]
            ct = np.swapaxes(coeffs,1,2)
            
            # zero pad to square dim, if necessary. nmax_grid >= m + 1
            shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
            shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
            x = sh.SHCoeffs.from_array(shs,normalization='4pi')
            
            # coordinate system exp transform to get surface pressure in Pa
            grid = x.expand(grid='DH')
            ps = np.exp(np.real(grid.to_array())) 
            p_field = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
            
            for level in range(l):
                # pressure coordinates of (u,v)
                p_field[:,:,level] = ps * slv[level]
                
            # u,v field from vorticity and divergence
            psi = np.copy(jvor)
            xi = np.copy(jdiv)
            
            # Del operator x-y coeffs from numerical recipe gauleg(ymin,ymax,ydim)
            mu = gauleg(-1.,1.,nmax_grid * 2)
            nui = 1 / np.sqrt(1 - mu**2) 
                
            for lvl in range(l):
                for n_sh in range(1,n + 1):
                    # n_sh runs from 1 to (n + 1); jvor(n=0) and jdiv(n=0) are (0,0j)
                    eigenvalue_leg = -n_sh * (n_sh + 1) / (a ** 2) # Saravanan p. 17
                    psi[:,n_sh,lvl] = 1 / eigenvalue_leg * jvor[:,n_sh,lvl]
                    xi[:,n_sh,lvl] = 1 / eigenvalue_leg * jdiv[:,n_sh,lvl]
                            
            # vorticity term
            for lvl_i in range(l):
                even = psi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                x = sh.SHCoeffs.from_array(shs,normalization='4pi')
                
                grid = x.expand(grid='DH')
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
        
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
        
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_x))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    # Saravanan recipe
                    x_coeffs[:,j] = nui / a 
                    y_coeffs[:,j] = 1 / a   
                
                # cross product in spherical coordinates
                u, v = uv_cross_product(x_coeffs * grad_x,
                                        y_coeffs * grad_y,
                                        lats,
                                        lons)
            
                u_arr[:,:,lvl_i] = u
                v_arr[:,:,lvl_i] = v
                
            # divergence term
            for lvl_i in range(l):
                even = xi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                # change m,n arangement to that used by shtools
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                x = sh.SHCoeffs.from_array(shs,normalization='4pi')
                
                grid = x.expand(grid='DH')
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
                    
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
                
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_y))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    x_coeffs[:,j] = nui / a
                    y_coeffs[:,j] = 1 / a 
                
                u = x_coeffs * (grad_x)
                v = y_coeffs * grad_y
            
                u_arr[:,:,lvl_i] = u_arr[:,:,lvl_i] + u 
                v_arr[:,:,lvl_i] = v_arr[:,:,lvl_i] + v
                
            # temperature to grid routine
            for lvl_i in range(l):
                even = jpot[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)    # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                x = sh.SHCoeffs.from_array(shs,normalization='4pi')
                
                grid = x.expand(grid='DH')
                data = np.real(grid.to_array())
            
                t_arr[:,:,lvl_i] = data
                
            # set dimensions to python standard [zdim,ydim,xdim], from y, u, z
            u_arr = np.swapaxes(u_arr,0,2) # z, u, y
            u_arr = np.swapaxes(u_arr,1,2) # z, y, u
            v_arr = np.swapaxes(v_arr,0,2)
            v_arr = np.swapaxes(v_arr,1,2)
            t_arr = np.swapaxes(t_arr,0,2)
            t_arr = np.swapaxes(t_arr,1,2)
            p_field = np.swapaxes(p_field,0,2)
            p_field = np.swapaxes(p_field,1,2)
        
            # calculate geometric height using t,p fields
            z_arr = pressure_to_altitude(p_field,
                                         t_arr)
            
            # interpolate to input levels (Pa)
            u_arr = interplevel(u_arr,z_arr,interp_levels)
            v_arr = interplevel(v_arr,z_arr,interp_levels)
                
            print('Reading record: ', t_steps)
            full_u[t_steps,:,:,:] = u_arr
            full_v[t_steps,:,:,:] = v_arr
            
        else:
            # y, x, z dim arrays to hold read in fields
            u_arr = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
            v_arr = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
            t_arr = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
            
            # read in spherical harmonic coefficients (shc)
            jvor = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            jdiv = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            jpot = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            pssp = f.read_record(dtype=np.complex64).reshape(m+1,n+1,order='F')
            f.read_record(dtype=np.float32).reshape(k1,k2,order='F') # unused
            
            # read in surface pressure field and construct odd shcs
            even = pssp[:,:]      
            odd = np.copy(even)   
            
            for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
            coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                              dtype=np.complex64)
            
            # change [m,n] arangement to [n,m] used by shtools
            coeffs[0,:,:] = even
            coeffs[1,1:,:] = odd[1:,:]
            ct = np.swapaxes(coeffs,1,2)
            
            # zero pad to square dim, if necessary. nmax_grid >= m + 1
            shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
            shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
            x = sh.SHCoeffs.from_array(shs,normalization='4pi')
            
            # coordinate system exp transform to get surface pressure in Pa
            grid = x.expand(grid='DH')
            ps = np.exp(np.real(grid.to_array())) 
            p_field = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
            
            for level in range(l):
                # pressure coordinates of (u,v)
                p_field[:,:,level] = ps * slv[level]
                
        
            # u,v field from vorticity and divergence
            psi = np.copy(jvor)
            xi = np.copy(jdiv)
            
            # Del operator x-y coeffs from numerical recipe gauleg(ymin,ymax,ydim)
            mu = gauleg(-1.,1.,nmax_grid * 2)
            nui = 1 / np.sqrt(1 - mu**2) 
                
            for lvl in range(l):
                for n_sh in range(1,n + 1):
                    # n_sh runs from 1 to (n + 1); jvor(n=0) and jdiv(n=0) are (0,0j)
                    eigenvalue_leg = -n_sh * (n_sh + 1) / (a ** 2) # Saravanan p. 17
                    psi[:,n_sh,lvl] = 1 / eigenvalue_leg * jvor[:,n_sh,lvl]
                    xi[:,n_sh,lvl] = 1 / eigenvalue_leg * jdiv[:,n_sh,lvl]
                            
            # vorticity term
            for lvl_i in range(l):
                even = psi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                x = sh.SHCoeffs.from_array(shs,normalization='4pi')
                
                grid = x.expand(grid='DH')
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
        
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
        
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_x))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    # Saravanan recipe
                    x_coeffs[:,j] = nui / a 
                    y_coeffs[:,j] = 1 / a   
                
                # cross product in spherical coordinates
                u, v = uv_cross_product(x_coeffs * grad_x,
                                        y_coeffs * grad_y,
                                        lats,
                                        lons)
            
                u_arr[:,:,lvl_i] = u
                v_arr[:,:,lvl_i] = v
                
            # divergence term
            for lvl_i in range(l):
                even = xi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                # change m,n arangement to that used by shtools
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                x = sh.SHCoeffs.from_array(shs,normalization='4pi')
                
                grid = x.expand(grid='DH')
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
                    
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
                
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_y))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    x_coeffs[:,j] = nui / a
                    y_coeffs[:,j] = 1 / a 
                
                u = x_coeffs * (grad_x)
                v = y_coeffs * grad_y
            
                u_arr[:,:,lvl_i] = u_arr[:,:,lvl_i] + u 
                v_arr[:,:,lvl_i] = v_arr[:,:,lvl_i] + v
                
            # temperature to grid routine
            for lvl_i in range(l):
                even = jpot[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)    # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                x = sh.SHCoeffs.from_array(shs,normalization='4pi')
                
                grid = x.expand(grid='DH')
                data = np.real(grid.to_array())
            
                t_arr[:,:,lvl_i] = data
                
            # set dimensions to python standard [zdim,ydim,xdim], from y, u, z
            u_arr = np.swapaxes(u_arr,0,2) # z, u, y
            u_arr = np.swapaxes(u_arr,1,2) # z, y, u
            v_arr = np.swapaxes(v_arr,0,2)
            v_arr = np.swapaxes(v_arr,1,2)
            t_arr = np.swapaxes(t_arr,0,2)
            t_arr = np.swapaxes(t_arr,1,2)
            p_field = np.swapaxes(p_field,0,2)
            p_field = np.swapaxes(p_field,1,2)
                            
            # calculate geometric height using t,p fields
            z_arr = pressure_to_altitude(p_field,
                                         t_arr)
            
            u_arr = interplevel(u_arr,z_arr,interp_levels)
            v_arr = interplevel(v_arr,z_arr,interp_levels)
                
            print('Reading record: ', t_steps)
            full_u[t_steps,:,:,:] = u_arr
            full_v[t_steps,:,:,:] = v_arr
                    
    f.close()
        
    return full_u, full_v


def uvt_read_altitude(file,
                      folder,
                      day_begin,
                      day_final,
                      grid_dim,
                      model_dt,
                      interp_levels=[100000]):
    """
    Uses hydrostatic equilibrium to convert model lvls to altitude
    
    """
    
    a = 6.37e6 # radius earth [m]    
    model_tstep = int(24 / model_dt)         
    
    # data arrays (square lat-lon) with size based on input values
    full_u = np.zeros((day_final * model_tstep,np.shape(interp_levels)[0],
                       grid_dim,grid_dim))
    full_v = np.zeros((day_final * model_tstep,np.shape(interp_levels)[0],
                       grid_dim,grid_dim))
    full_t = np.zeros((day_final * model_tstep,np.shape(interp_levels)[0],
                       grid_dim,grid_dim))
    
    f = FortranFile(folder + file,'r')
    
    # read in grid dimensions 
    md = f.read_record(dtype=np.int32)
    (k1, k2, m, n, l) = (md[0], md[1], md[2], md[3], md[4])
            
    # model levels in sigma coordinates
    slv = f.read_record(dtype=np.float32)
        
    # latitudes and zsdtlv (unused)
    lat_sigprim = f.read_record(dtype=np.float32)
    f.read_record(dtype=np.float32)
       
    nmax_grid = int(grid_dim / 2) # max n + 1, controls output grid resolution
        
    for t_steps in range(day_final * model_tstep):
        if t_steps == 0:
            # y, x, z dim arrays to hold read in fields
            u_arr = np.zeros((grid_dim,grid_dim,l))
            v_arr = np.zeros((grid_dim,grid_dim,l))
            t_arr = np.zeros((grid_dim,grid_dim,l))
                        
            for k in range(day_begin * model_tstep + 1):
                jvor = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                jdiv = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                jpot = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                pssp = f.read_record(dtype=np.complex64).reshape(m+1,n+1,order='F')
                geop = f.read_record(dtype=np.float32).reshape(k1,k2,order='F') 
            
            # read in surface pressure field and construct odd shcs
            even = pssp[:,:]      
            odd = np.copy(even)   
            
            for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
            coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                              dtype=np.complex64)
            
            # change [m,n] arangement to [n,m] used by shtools
            coeffs[0,:,:] = even
            coeffs[1,1:,:] = odd[1:,:]
            ct = np.swapaxes(coeffs,1,2)
            
            # zero pad to square dim, if necessary. nmax_grid >= m + 1
            shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
            shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
            
            # Condon-Shortley phase is used (-1) or not (1)
            x = sh.SHCoeffs.from_array(shs,                                       
                                       normalization='4pi',
                                       csphase=1)
                        
            # coordinate system exp transform to get surface pressure in Pa
            grid = x.expand(grid='DH',extend=False)
            ps = np.exp(np.real(grid.to_array()))
            p_field = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
                        
            for level in range(l):
                # pressure coordinates of (u,v)
                p_field[:,:,level] = ps * slv[level]
                
            # u,v field from vorticity and divergence
            psi = np.copy(jvor)
            xi = np.copy(jdiv)
            
            # Del operator x-y coeffs from numerical recipe gauleg(ymin,ymax,ydim)
            mu = gauleg(-1.,1.,nmax_grid * 2)
            nui = 1 / np.sqrt(1 - mu**2) 
                
            for lvl in range(l):
                for n_sh in range(1,n + 1):
                    # n_sh runs from 1 to (n + 1); jvor(n=0) and jdiv(n=0) are (0,0j)
                    eigenvalue_leg = -n_sh * (n_sh + 1) / (a ** 2) # Saravanan p. 17
                    psi[:,n_sh,lvl] = 1 / eigenvalue_leg * jvor[:,n_sh,lvl]
                    xi[:,n_sh,lvl] = 1 / eigenvalue_leg * jdiv[:,n_sh,lvl]
                            
            # vorticity term
            for lvl_i in range(l):
                even = psi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH',extend=False)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
        
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
        
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_x))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    # Saravanan recipe
                    x_coeffs[:,j] = nui / a 
                    y_coeffs[:,j] = 1 / a   
                
                # cross product in spherical coordinates
                u, v = uv_cross_product(x_coeffs * grad_x,
                                        y_coeffs * grad_y,
                                        lats,
                                        lons)
            
                u_arr[:,:,lvl_i] = u
                v_arr[:,:,lvl_i] = v 
                                
            # divergence term
            for lvl_i in range(l):
                even = xi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                # change m,n arangement to that used by shtools
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH',extend=False)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
                    
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
                
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_y))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    x_coeffs[:,j] = nui / a
                    y_coeffs[:,j] = 1 / a 
                
                u = x_coeffs * (grad_x)
                v = y_coeffs * grad_y
                            
                u_arr[:,:,lvl_i] = u_arr[:,:,lvl_i] + u 
                v_arr[:,:,lvl_i] = v_arr[:,:,lvl_i] + v
                                
            # temperature to grid routine
            for lvl_i in range(l):
                even = jpot[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)    # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH',extend=False)
                data = np.real(grid.to_array())
            
                t_arr[:,:,lvl_i] = data
                
            # set dimensions to python standard [zdim,ydim,xdim], from y, u, z
            u_arr = np.swapaxes(u_arr,0,2) # z, u, y
            u_arr = np.swapaxes(u_arr,1,2) # z, y, u
            v_arr = np.swapaxes(v_arr,0,2)
            v_arr = np.swapaxes(v_arr,1,2)
            t_arr = np.swapaxes(t_arr,0,2)
            t_arr = np.swapaxes(t_arr,1,2)
            p_field = np.swapaxes(p_field,0,2)
            p_field = np.swapaxes(p_field,1,2)
            
            # set geopotential dimensions for python standard
            geop = np.swapaxes(geop[:,::-1],0,1)
            
            # calculate geometric height using t,p fields
            z_arr = pressure_to_altitude_topography(p_field,
                                                    t_arr,
                                                    geop,
                                                    lat_sigprim)
            
            # # calculate geometric height using t,p fields
            # z_arr = hydro_ideal(p_field,
            #                     t_arr,
            #                     geop,
            #                     lat_sigprim)
                                    
            # interpolate to input levels (Pa)
            u_arr = interplevel(u_arr,z_arr,interp_levels)
            v_arr = interplevel(v_arr,z_arr,interp_levels)
            t_arr = interplevel(t_arr,z_arr,interp_levels)
                
            print('Reading record: ', t_steps)
            full_u[t_steps,:,:,:] = u_arr
            full_v[t_steps,:,:,:] = v_arr
            full_t[t_steps,:,:,:] = t_arr
            
        else:
            # y, x, z dim arrays to hold read in fields
            u_arr = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
            v_arr = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
            t_arr = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
            
            # read in spherical harmonic coefficients (shc)
            jvor = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            jdiv = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            jpot = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            pssp = f.read_record(dtype=np.complex64).reshape(m+1,n+1,order='F')
            geop = f.read_record(dtype=np.float32).reshape(k1,k2,order='F') 
            
            # read in surface pressure field and construct odd shcs
            even = pssp[:,:]      
            odd = np.copy(even)   
            
            for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
            coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                              dtype=np.complex64)
            
            # change [m,n] arangement to [n,m] used by shtools
            coeffs[0,:,:] = even
            coeffs[1,1:,:] = odd[1:,:]
            ct = np.swapaxes(coeffs,1,2)
            
            # zero pad to square dim, if necessary. nmax_grid >= m + 1
            shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
            shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
            
            # Condon-Shortley phase is used (-1) or not (1)
            x = sh.SHCoeffs.from_array(shs,
                                       normalization='4pi',
                                       csphase=1)
            
            # coordinate system exp transform to get surface pressure in Pa
            grid = x.expand(grid='DH',extend=False)
            ps = np.exp(np.real(grid.to_array())) 
            p_field = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
            
            for level in range(l):
                # pressure coordinates of (u,v)
                p_field[:,:,level] = ps * slv[level]
                                
            # u,v field from vorticity and divergence
            psi = np.copy(jvor)
            xi = np.copy(jdiv)
            
            # Del operator x-y coeffs from numerical recipe gauleg(ymin,ymax,ydim)
            mu = gauleg(-1.,1.,nmax_grid * 2)
            nui = 1 / np.sqrt(1 - mu**2) 
                
            for lvl in range(l):
                for n_sh in range(1,n + 1):
                    # n_sh runs from 1 to (n + 1); jvor(n=0) and jdiv(n=0) are (0,0j)
                    eigenvalue_leg = -n_sh * (n_sh + 1) / (a ** 2) # Saravanan p. 17
                    psi[:,n_sh,lvl] = 1 / eigenvalue_leg * jvor[:,n_sh,lvl]
                    xi[:,n_sh,lvl] = 1 / eigenvalue_leg * jdiv[:,n_sh,lvl]
                            
            # vorticity term
            for lvl_i in range(l):
                even = psi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH',extend=False)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
        
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
        
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_x))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    # Saravanan recipe
                    x_coeffs[:,j] = nui / a 
                    y_coeffs[:,j] = 1 / a   
                
                # cross product in spherical coordinates
                u, v = uv_cross_product(x_coeffs * grad_x,
                                        y_coeffs * grad_y,
                                        lats,
                                        lons)
            
                u_arr[:,:,lvl_i] = u
                v_arr[:,:,lvl_i] = v
                            
            # divergence term
            for lvl_i in range(l):
                even = xi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                # change m,n arangement to that used by shtools
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH',extend=False)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
                    
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
                
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_y))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    x_coeffs[:,j] = nui / a
                    y_coeffs[:,j] = 1 / a 
                
                u = x_coeffs * (grad_x)
                v = y_coeffs * grad_y
            
                u_arr[:,:,lvl_i] = u_arr[:,:,lvl_i] + u 
                v_arr[:,:,lvl_i] = v_arr[:,:,lvl_i] + v
                
            # temperature to grid routine
            for lvl_i in range(l):
                even = jpot[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)    # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH',extend=False)
                data = np.real(grid.to_array())
            
                t_arr[:,:,lvl_i] = data
                
            # set dimensions to python standard [zdim,ydim,xdim], from y, u, z
            u_arr = np.swapaxes(u_arr,0,2) # z, u, y
            u_arr = np.swapaxes(u_arr,1,2) # z, y, u
            v_arr = np.swapaxes(v_arr,0,2)
            v_arr = np.swapaxes(v_arr,1,2)
            t_arr = np.swapaxes(t_arr,0,2)
            t_arr = np.swapaxes(t_arr,1,2)
            p_field = np.swapaxes(p_field,0,2)
            p_field = np.swapaxes(p_field,1,2)
            
            # set geopotential dimensions for python standard
            geop = np.swapaxes(geop[:,::-1],0,1)
                        
            # calculate geometric height using t,p fields
            z_arr = pressure_to_altitude_topography(p_field,
                                                    t_arr,
                                                    geop,
                                                    lat_sigprim)
            
            # # calculate geometric height using t,p fields
            # z_arr = hydro_ideal(p_field,
            #                     t_arr,
            #                     geop,
            #                     lat_sigprim)
            
            u_arr = interplevel(u_arr,z_arr,interp_levels)
            v_arr = interplevel(v_arr,z_arr,interp_levels)
            t_arr = interplevel(t_arr,z_arr,interp_levels)
                
            print('Reading record: ', t_steps)
            full_u[t_steps,:,:,:] = u_arr
            full_v[t_steps,:,:,:] = v_arr
            full_t[t_steps,:,:,:] = t_arr
                    
    f.close()
        
    return full_u, full_v, full_t


def uvt_read_pressure(file,
                      folder,
                      day_begin,
                      day_final,
                      grid_dim,
                      model_dt,
                      interp_levels=[100000]):
    """
    Uses hydrostatic equilibrium to convert model lvls to altitude
    
    """
    
    a = 6.37e6 # radius earth [m]    
    model_tstep = int(24 / model_dt)         
    
    # data arrays (square lat-lon) with size based on input values
    full_u = np.zeros((day_final * model_tstep,np.shape(interp_levels)[0],
                       grid_dim,grid_dim))
    full_v = np.zeros((day_final * model_tstep,np.shape(interp_levels)[0],
                       grid_dim,grid_dim))
    full_t = np.zeros((day_final * model_tstep,np.shape(interp_levels)[0],
                       grid_dim,grid_dim))
    
    f = FortranFile(folder + file,'r')
    
    # read in grid dimensions 
    md = f.read_record(dtype=np.int32)
    (k1, k2, m, n, l) = (md[0], md[1], md[2], md[3], md[4])
            
    # model levels in sigma coordinates
    slv = f.read_record(dtype=np.float32)
        
    # latitudes and zsdtlv (unused)
    lat_sigprim = f.read_record(dtype=np.float32)
    f.read_record(dtype=np.float32)
       
    nmax_grid = int(grid_dim / 2) # max n + 1, controls output grid resolution
        
    for t_steps in range(day_final * model_tstep):
        if t_steps == 0:
            # y, x, z dim arrays to hold read in fields
            u_arr = np.zeros((grid_dim,grid_dim,l))
            v_arr = np.zeros((grid_dim,grid_dim,l))
            t_arr = np.zeros((grid_dim,grid_dim,l))
                        
            for k in range(day_begin * model_tstep + 1):
                jvor = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                jdiv = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                jpot = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                pssp = f.read_record(dtype=np.complex64).reshape(m+1,n+1,order='F')
                geop = f.read_record(dtype=np.float32).reshape(k1,k2,order='F') 
            
            # read in surface pressure field and construct odd shcs
            even = pssp[:,:]      
            odd = np.copy(even)   
            
            for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
            coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                              dtype=np.complex64)
            
            # change [m,n] arangement to [n,m] used by shtools
            coeffs[0,:,:] = even
            coeffs[1,1:,:] = odd[1:,:]
            ct = np.swapaxes(coeffs,1,2)
            
            # zero pad to square dim, if necessary. nmax_grid >= m + 1
            shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
            shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
            
            # Condon-Shortley phase is used (-1) or not (1)
            x = sh.SHCoeffs.from_array(shs,                                       
                                       normalization='4pi',
                                       csphase=1)
                        
            # coordinate system exp transform to get surface pressure in Pa
            grid = x.expand(grid='DH',extend=False)
            ps = np.exp(np.real(grid.to_array()))
            p_field = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
                        
            for level in range(l):
                # pressure coordinates of (u,v)
                p_field[:,:,level] = ps * slv[level]
                
            # u,v field from vorticity and divergence
            psi = np.copy(jvor)
            xi = np.copy(jdiv)
            
            # Del operator x-y coeffs from numerical recipe gauleg(ymin,ymax,ydim)
            mu = gauleg(-1.,1.,nmax_grid * 2)
            nui = 1 / np.sqrt(1 - mu**2) 
                
            for lvl in range(l):
                for n_sh in range(1,n + 1):
                    # n_sh runs from 1 to (n + 1); jvor(n=0) and jdiv(n=0) are (0,0j)
                    eigenvalue_leg = -n_sh * (n_sh + 1) / (a ** 2) # Saravanan p. 17
                    psi[:,n_sh,lvl] = 1 / eigenvalue_leg * jvor[:,n_sh,lvl]
                    xi[:,n_sh,lvl] = 1 / eigenvalue_leg * jdiv[:,n_sh,lvl]
                            
            # vorticity term
            for lvl_i in range(l):
                even = psi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH',extend=False)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
        
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
        
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_x))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    # Saravanan recipe
                    x_coeffs[:,j] = nui / a 
                    y_coeffs[:,j] = 1 / a   
                
                # cross product in spherical coordinates
                u, v = uv_cross_product(x_coeffs * grad_x,
                                        y_coeffs * grad_y,
                                        lats,
                                        lons)
            
                u_arr[:,:,lvl_i] = u
                v_arr[:,:,lvl_i] = v 
                                
            # divergence term
            for lvl_i in range(l):
                even = xi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                # change m,n arangement to that used by shtools
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH',extend=False)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
                    
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
                
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_y))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    x_coeffs[:,j] = nui / a
                    y_coeffs[:,j] = 1 / a 
                
                u = x_coeffs * (grad_x)
                v = y_coeffs * grad_y
                            
                u_arr[:,:,lvl_i] = u_arr[:,:,lvl_i] + u 
                v_arr[:,:,lvl_i] = v_arr[:,:,lvl_i] + v
                                
            # temperature to grid routine
            for lvl_i in range(l):
                even = jpot[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)    # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH',extend=False)
                data = np.real(grid.to_array())
            
                t_arr[:,:,lvl_i] = data
                
            # set dimensions to python standard [zdim,ydim,xdim], from y, u, z
            u_arr = np.swapaxes(u_arr,0,2) # z, u, y
            u_arr = np.swapaxes(u_arr,1,2) # z, y, u
            v_arr = np.swapaxes(v_arr,0,2)
            v_arr = np.swapaxes(v_arr,1,2)
            t_arr = np.swapaxes(t_arr,0,2)
            t_arr = np.swapaxes(t_arr,1,2)
            p_field = np.swapaxes(p_field,0,2)
            p_field = np.swapaxes(p_field,1,2)
            
            # set geopotential dimensions for python standard
            geop = np.swapaxes(geop[:,::-1],0,1)
            
            # calculate geometric height using t,p fields
            z_arr = p_field
                                    
            # interpolate to input levels (Pa)
            u_arr = interplevel(u_arr,z_arr,interp_levels)
            v_arr = interplevel(v_arr,z_arr,interp_levels)
            t_arr = interplevel(t_arr,z_arr,interp_levels)
                
            print('Reading record: ', t_steps)
            full_u[t_steps,:,:,:] = u_arr
            full_v[t_steps,:,:,:] = v_arr
            full_t[t_steps,:,:,:] = t_arr
            
        else:
            # y, x, z dim arrays to hold read in fields
            u_arr = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
            v_arr = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
            t_arr = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
            
            # read in spherical harmonic coefficients (shc)
            jvor = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            jdiv = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            jpot = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            pssp = f.read_record(dtype=np.complex64).reshape(m+1,n+1,order='F')
            geop = f.read_record(dtype=np.float32).reshape(k1,k2,order='F') 
            
            # read in surface pressure field and construct odd shcs
            even = pssp[:,:]      
            odd = np.copy(even)   
            
            for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
            coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                              dtype=np.complex64)
            
            # change [m,n] arangement to [n,m] used by shtools
            coeffs[0,:,:] = even
            coeffs[1,1:,:] = odd[1:,:]
            ct = np.swapaxes(coeffs,1,2)
            
            # zero pad to square dim, if necessary. nmax_grid >= m + 1
            shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
            shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
            
            # Condon-Shortley phase is used (-1) or not (1)
            x = sh.SHCoeffs.from_array(shs,
                                       normalization='4pi',
                                       csphase=1)
            
            # coordinate system exp transform to get surface pressure in Pa
            grid = x.expand(grid='DH',extend=False)
            ps = np.exp(np.real(grid.to_array())) 
            p_field = np.zeros((nmax_grid * 2,nmax_grid * 2, l))
            
            for level in range(l):
                # pressure coordinates of (u,v)
                p_field[:,:,level] = ps * slv[level]
                                
            # u,v field from vorticity and divergence
            psi = np.copy(jvor)
            xi = np.copy(jdiv)
            
            # Del operator x-y coeffs from numerical recipe gauleg(ymin,ymax,ydim)
            mu = gauleg(-1.,1.,nmax_grid * 2)
            nui = 1 / np.sqrt(1 - mu**2) 
                
            for lvl in range(l):
                for n_sh in range(1,n + 1):
                    # n_sh runs from 1 to (n + 1); jvor(n=0) and jdiv(n=0) are (0,0j)
                    eigenvalue_leg = -n_sh * (n_sh + 1) / (a ** 2) # Saravanan p. 17
                    psi[:,n_sh,lvl] = 1 / eigenvalue_leg * jvor[:,n_sh,lvl]
                    xi[:,n_sh,lvl] = 1 / eigenvalue_leg * jdiv[:,n_sh,lvl]
                            
            # vorticity term
            for lvl_i in range(l):
                even = psi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH',extend=False)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
        
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
        
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_x))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    # Saravanan recipe
                    x_coeffs[:,j] = nui / a 
                    y_coeffs[:,j] = 1 / a   
                
                # cross product in spherical coordinates
                u, v = uv_cross_product(x_coeffs * grad_x,
                                        y_coeffs * grad_y,
                                        lats,
                                        lons)
            
                u_arr[:,:,lvl_i] = u
                v_arr[:,:,lvl_i] = v
                            
            # divergence term
            for lvl_i in range(l):
                even = xi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                # change m,n arangement to that used by shtools
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH',extend=False)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
                    
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
                
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_y))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    x_coeffs[:,j] = nui / a
                    y_coeffs[:,j] = 1 / a 
                
                u = x_coeffs * (grad_x)
                v = y_coeffs * grad_y
            
                u_arr[:,:,lvl_i] = u_arr[:,:,lvl_i] + u 
                v_arr[:,:,lvl_i] = v_arr[:,:,lvl_i] + v
                
            # temperature to grid routine
            for lvl_i in range(l):
                even = jpot[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)    # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH',extend=False)
                data = np.real(grid.to_array())
            
                t_arr[:,:,lvl_i] = data
                
            # set dimensions to python standard [zdim,ydim,xdim], from y, u, z
            u_arr = np.swapaxes(u_arr,0,2) # z, u, y
            u_arr = np.swapaxes(u_arr,1,2) # z, y, u
            v_arr = np.swapaxes(v_arr,0,2)
            v_arr = np.swapaxes(v_arr,1,2)
            t_arr = np.swapaxes(t_arr,0,2)
            t_arr = np.swapaxes(t_arr,1,2)
            p_field = np.swapaxes(p_field,0,2)
            p_field = np.swapaxes(p_field,1,2)
            
            # set geopotential dimensions for python standard
            geop = np.swapaxes(geop[:,::-1],0,1)
                        
            # calculate geometric height using t,p fields
            z_arr = p_field
            
            u_arr = interplevel(u_arr,z_arr,interp_levels)
            v_arr = interplevel(v_arr,z_arr,interp_levels)
            t_arr = interplevel(t_arr,z_arr,interp_levels)
                
            print('Reading record: ', t_steps)
            full_u[t_steps,:,:,:] = u_arr
            full_v[t_steps,:,:,:] = v_arr
            full_t[t_steps,:,:,:] = t_arr
                    
    f.close()
        
    return full_u, full_v, full_t

def uvt_read_altitude_metradar(file,
                              folder,
                              day_begin,
                              day_final,
                              grid_dim,
                              interp_levels=[100000]):
    """
    Uses hydrostatic equilibrium to convert model lvls to altitude
    
    """
    
    a = 6.37e6 # radius earth [m]        
    
    # data arrays (square lat-lon) with size based on input values
    full_u = np.zeros((day_final * 8,np.shape(interp_levels)[0],
                       grid_dim,grid_dim * 2))
    full_v = np.zeros((day_final * 8,np.shape(interp_levels)[0],
                       grid_dim,grid_dim * 2))
    full_t = np.zeros((day_final * 8,np.shape(interp_levels)[0],
                       grid_dim,grid_dim * 2))
    
    f = FortranFile(folder + file,'r')
    
    # read in grid dimensions 
    md = f.read_record(dtype=np.int32)
    (k1, k2, m, n, l) = (md[0], md[1], md[2], md[3], md[4])
            
    # model levels in sigma coordinates
    slv = f.read_record(dtype=np.float32)
        
    # latitudes and zsdtlv (unused)
    lat_sigprim = f.read_record(dtype=np.float32)
    f.read_record(dtype=np.float32)
       
    nmax_grid = int(grid_dim / 2) # max n + 1, controls output grid resolution
        
    for t_steps in range(day_final * 8):
        if t_steps == 0:
            # y, x, z dim arrays to hold read in fields
            u_arr = np.zeros((grid_dim,grid_dim * 2,l))
            v_arr = np.zeros((grid_dim,grid_dim * 2,l))
            t_arr = np.zeros((grid_dim,grid_dim * 2,l))
                        
            for k in range(day_begin * 8 + 1):
                jvor = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                jdiv = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                jpot = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                pssp = f.read_record(dtype=np.complex64).reshape(m+1,n+1,order='F')
                geop = f.read_record(dtype=np.float32).reshape(k1,k2,order='F') 
            
            # read in surface pressure field and construct odd shcs
            even = pssp[:,:]      
            odd = np.copy(even)   
            
            for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
            coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                              dtype=np.complex64)
            
            # change [m,n] arangement to [n,m] used by shtools
            coeffs[0,:,:] = even
            coeffs[1,1:,:] = odd[1:,:]
            ct = np.swapaxes(coeffs,1,2)
            
            # zero pad to square dim, if necessary. nmax_grid >= m + 1
            shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
            shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
            
            # Condon-Shortley phase is used (-1) or not (1)
            x = sh.SHCoeffs.from_array(shs,                                       
                                       normalization='4pi',
                                       csphase=1)
                        
            # coordinate system exp transform to get surface pressure in Pa
            grid = x.expand(grid='DH2',extend=False)
            ps = np.exp(np.real(grid.to_array()))
            p_field = np.zeros((nmax_grid * 2,nmax_grid * 2  * 2, l))
                        
            for level in range(l):
                # pressure coordinates of (u,v)
                p_field[:,:,level] = ps * slv[level]
                
            # u,v field from vorticity and divergence
            psi = np.copy(jvor)
            xi = np.copy(jdiv)
            
            # Del operator x-y coeffs from numerical recipe gauleg(ymin,ymax,ydim)
            mu = gauleg(-1.,1.,nmax_grid * 2)
            nui = 1 / np.sqrt(1 - mu**2) 
                
            for lvl in range(l):
                for n_sh in range(1,n + 1):
                    # n_sh runs from 1 to (n + 1); jvor(n=0) and jdiv(n=0) are (0,0j)
                    eigenvalue_leg = -n_sh * (n_sh + 1) / (a ** 2) # Saravanan p. 17
                    psi[:,n_sh,lvl] = 1 / eigenvalue_leg * jvor[:,n_sh,lvl]
                    xi[:,n_sh,lvl] = 1 / eigenvalue_leg * jdiv[:,n_sh,lvl]
                            
            # vorticity term
            for lvl_i in range(l):
                even = psi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH2',extend=False)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
        
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
        
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_x))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    # Saravanan recipe
                    x_coeffs[:,j] = nui / a 
                    y_coeffs[:,j] = 1 / a   
                
                # cross product in spherical coordinates
                u, v = uv_cross_product(x_coeffs * grad_x,
                                        y_coeffs * grad_y,
                                        lats,
                                        lons)
            
                u_arr[:,:,lvl_i] = u
                v_arr[:,:,lvl_i] = v 
                                
            # divergence term
            for lvl_i in range(l):
                even = xi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                # change m,n arangement to that used by shtools
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH2',extend=False)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
                    
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
                
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_y))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    x_coeffs[:,j] = nui / a
                    y_coeffs[:,j] = 1 / a 
                
                u = x_coeffs * (grad_x)
                v = y_coeffs * grad_y
                            
                u_arr[:,:,lvl_i] = u_arr[:,:,lvl_i] + u 
                v_arr[:,:,lvl_i] = v_arr[:,:,lvl_i] + v
                                
            # temperature to grid routine
            for lvl_i in range(l):
                even = jpot[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)    # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH2',extend=False)
                data = np.real(grid.to_array())
            
                t_arr[:,:,lvl_i] = data
                
            # set dimensions to python standard [zdim,ydim,xdim], from y, u, z
            u_arr = np.swapaxes(u_arr,0,2) # z, u, y
            u_arr = np.swapaxes(u_arr,1,2) # z, y, u
            v_arr = np.swapaxes(v_arr,0,2)
            v_arr = np.swapaxes(v_arr,1,2)
            t_arr = np.swapaxes(t_arr,0,2)
            t_arr = np.swapaxes(t_arr,1,2)
            p_field = np.swapaxes(p_field,0,2)
            p_field = np.swapaxes(p_field,1,2)
            
            # set geopotential dimensions for python standard
            geop = np.swapaxes(geop[:,::-1],0,1)
            
            # calculate geometric height using t,p fields
            z_arr = hydro_ideal(p_field,
                                t_arr,
                                geop,
                                lat_sigprim)
            
            print(np.mean(z_arr,axis=(1,2))[0] / 1000)
            print(np.mean(p_field,axis=(1,2))[0] / 100)
            
            print(np.mean(z_arr,axis=(1,2))[30] / 1000)
            print(np.mean(p_field,axis=(1,2))[30] / 100)
            
            print(np.mean(z_arr,axis=(1,2))[38] / 1000)
            print(np.mean(p_field,axis=(1,2))[38] / 100)
            
                                    
            # interpolate to input levels (Pa)
            u_arr = interplevel(u_arr,z_arr,interp_levels)
            v_arr = interplevel(v_arr,z_arr,interp_levels)
            t_arr = interplevel(t_arr,z_arr,interp_levels)
                
            print('Reading record: ', t_steps)
            full_u[t_steps,:,:,:] = u_arr
            full_v[t_steps,:,:,:] = v_arr
            full_t[t_steps,:,:,:] = t_arr
            
        else:
            # y, x, z dim arrays to hold read in fields
            u_arr = np.zeros((nmax_grid * 2,nmax_grid * 2 * 2, l))
            v_arr = np.zeros((nmax_grid * 2,nmax_grid * 2 * 2, l))
            t_arr = np.zeros((nmax_grid * 2,nmax_grid * 2 * 2, l))
            
            # read in spherical harmonic coefficients (shc)
            jvor = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            jdiv = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            jpot = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            pssp = f.read_record(dtype=np.complex64).reshape(m+1,n+1,order='F')
            geop = f.read_record(dtype=np.float32).reshape(k1,k2,order='F') 
            
            # read in surface pressure field and construct odd shcs
            even = pssp[:,:]      
            odd = np.copy(even)   
            
            for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
            coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                              dtype=np.complex64)
            
            # change [m,n] arangement to [n,m] used by shtools
            coeffs[0,:,:] = even
            coeffs[1,1:,:] = odd[1:,:]
            ct = np.swapaxes(coeffs,1,2)
            
            # zero pad to square dim, if necessary. nmax_grid >= m + 1
            shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
            shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
            
            # Condon-Shortley phase is used (-1) or not (1)
            x = sh.SHCoeffs.from_array(shs,
                                       normalization='4pi',
                                       csphase=1)
            
            # coordinate system exp transform to get surface pressure in Pa
            grid = x.expand(grid='DH2',extend=False)
            ps = np.exp(np.real(grid.to_array())) 
            p_field = np.zeros((nmax_grid * 2,nmax_grid * 2 * 2, l))
            
            for level in range(l):
                # pressure coordinates of (u,v)
                p_field[:,:,level] = ps * slv[level]
                                
            # u,v field from vorticity and divergence
            psi = np.copy(jvor)
            xi = np.copy(jdiv)
            
            # Del operator x-y coeffs from numerical recipe gauleg(ymin,ymax,ydim)
            mu = gauleg(-1.,1.,nmax_grid * 2)
            nui = 1 / np.sqrt(1 - mu**2) 
                
            for lvl in range(l):
                for n_sh in range(1,n + 1):
                    # n_sh runs from 1 to (n + 1); jvor(n=0) and jdiv(n=0) are (0,0j)
                    eigenvalue_leg = -n_sh * (n_sh + 1) / (a ** 2) # Saravanan p. 17
                    psi[:,n_sh,lvl] = 1 / eigenvalue_leg * jvor[:,n_sh,lvl]
                    xi[:,n_sh,lvl] = 1 / eigenvalue_leg * jdiv[:,n_sh,lvl]
                            
            # vorticity term
            for lvl_i in range(l):
                even = psi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH2',extend=False)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
        
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
        
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_x))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    # Saravanan recipe
                    x_coeffs[:,j] = nui / a 
                    y_coeffs[:,j] = 1 / a   
                
                # cross product in spherical coordinates
                u, v = uv_cross_product(x_coeffs * grad_x,
                                        y_coeffs * grad_y,
                                        lats,
                                        lons)
            
                u_arr[:,:,lvl_i] = u
                v_arr[:,:,lvl_i] = v
                            
            # divergence term
            for lvl_i in range(l):
                even = xi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                # change m,n arangement to that used by shtools
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH2',extend=False)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
                    
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
                
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_y))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    x_coeffs[:,j] = nui / a
                    y_coeffs[:,j] = 1 / a 
                
                u = x_coeffs * (grad_x)
                v = y_coeffs * grad_y
            
                u_arr[:,:,lvl_i] = u_arr[:,:,lvl_i] + u 
                v_arr[:,:,lvl_i] = v_arr[:,:,lvl_i] + v
                
            # temperature to grid routine
            for lvl_i in range(l):
                even = jpot[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)    # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH2',extend=False)
                data = np.real(grid.to_array())
            
                t_arr[:,:,lvl_i] = data
                
            # set dimensions to python standard [zdim,ydim,xdim], from y, u, z
            u_arr = np.swapaxes(u_arr,0,2) # z, u, y
            u_arr = np.swapaxes(u_arr,1,2) # z, y, u
            v_arr = np.swapaxes(v_arr,0,2)
            v_arr = np.swapaxes(v_arr,1,2)
            t_arr = np.swapaxes(t_arr,0,2)
            t_arr = np.swapaxes(t_arr,1,2)
            p_field = np.swapaxes(p_field,0,2)
            p_field = np.swapaxes(p_field,1,2)
            
            # set geopotential dimensions for python standard
            geop = np.swapaxes(geop[:,::-1],0,1)
                        
            # calculate geometric height using t,p fields
            z_arr = hydro_ideal(p_field,
                                t_arr,
                                geop,
                                lat_sigprim)
                        
            u_arr = interplevel(u_arr,z_arr,interp_levels)
            v_arr = interplevel(v_arr,z_arr,interp_levels)
            t_arr = interplevel(t_arr,z_arr,interp_levels)
                
            print('Reading record: ', t_steps)
            full_u[t_steps,:,:,:] = u_arr
            full_v[t_steps,:,:,:] = v_arr
            full_t[t_steps,:,:,:] = t_arr
                    
    f.close()
        
    return full_u, full_v, full_t

def uvt_read_altitude_scaleheight(file,
                                  folder,
                                  day_begin,
                                  day_final,
                                  grid_dim,
                                  interp_levels=[100000],
                                  H=7000):
    """
    Uses hydrostatic equilibrium to convert model lvls to altitude
    
    """
    
    a = 6.37e6 # radius earth [m]        
    
    # data arrays (square lat-lon) with size based on input values
    full_u = np.zeros((day_final * 8,np.shape(interp_levels)[0],
                       grid_dim,grid_dim * 2))
    full_v = np.zeros((day_final * 8,np.shape(interp_levels)[0],
                       grid_dim,grid_dim * 2))
    full_t = np.zeros((day_final * 8,np.shape(interp_levels)[0],
                       grid_dim,grid_dim * 2))
    
    f = FortranFile(folder + file,'r')
    
    # read in grid dimensions 
    md = f.read_record(dtype=np.int32)
    (k1, k2, m, n, l) = (md[0], md[1], md[2], md[3], md[4])
            
    # model levels in sigma coordinates
    slv = f.read_record(dtype=np.float32)
        
    # latitudes and zsdtlv (unused)
    lat_sigprim = f.read_record(dtype=np.float32)
    f.read_record(dtype=np.float32)
       
    nmax_grid = int(grid_dim / 2) # max n + 1, controls output grid resolution
        
    for t_steps in range(day_final * 8):
        if t_steps == 0:
            # y, x, z dim arrays to hold read in fields
            u_arr = np.zeros((grid_dim,grid_dim * 2,l))
            v_arr = np.zeros((grid_dim,grid_dim * 2,l))
            t_arr = np.zeros((grid_dim,grid_dim * 2,l))
                        
            for k in range(day_begin * 8 + 1):
                jvor = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                jdiv = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                jpot = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
                pssp = f.read_record(dtype=np.complex64).reshape(m+1,n+1,order='F')
                geop = f.read_record(dtype=np.float32).reshape(k1,k2,order='F') 
            
            # read in surface pressure field and construct odd shcs
            even = pssp[:,:]      
            odd = np.copy(even)   
            
            for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
            coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                              dtype=np.complex64)
            
            # change [m,n] arangement to [n,m] used by shtools
            coeffs[0,:,:] = even
            coeffs[1,1:,:] = odd[1:,:]
            ct = np.swapaxes(coeffs,1,2)
            
            # zero pad to square dim, if necessary. nmax_grid >= m + 1
            shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
            shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
            
            # Condon-Shortley phase is used (-1) or not (1)
            x = sh.SHCoeffs.from_array(shs,                                       
                                       normalization='4pi',
                                       csphase=1)
                        
            # coordinate system exp transform to get surface pressure in Pa
            grid = x.expand(grid='DH2',extend=False)
            ps = np.exp(np.real(grid.to_array()))
            p_field = np.zeros((nmax_grid * 2,nmax_grid * 2  * 2, l))
                        
            for level in range(l):
                # pressure coordinates of (u,v)
                p_field[:,:,level] = ps * slv[level]
                
            # u,v field from vorticity and divergence
            psi = np.copy(jvor)
            xi = np.copy(jdiv)
            
            # Del operator x-y coeffs from numerical recipe gauleg(ymin,ymax,ydim)
            mu = gauleg(-1.,1.,nmax_grid * 2)
            nui = 1 / np.sqrt(1 - mu**2) 
                
            for lvl in range(l):
                for n_sh in range(1,n + 1):
                    # n_sh runs from 1 to (n + 1); jvor(n=0) and jdiv(n=0) are (0,0j)
                    eigenvalue_leg = -n_sh * (n_sh + 1) / (a ** 2) # Saravanan p. 17
                    psi[:,n_sh,lvl] = 1 / eigenvalue_leg * jvor[:,n_sh,lvl]
                    xi[:,n_sh,lvl] = 1 / eigenvalue_leg * jdiv[:,n_sh,lvl]
                            
            # vorticity term
            for lvl_i in range(l):
                even = psi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH2',extend=False)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
        
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
        
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_x))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    # Saravanan recipe
                    x_coeffs[:,j] = nui / a 
                    y_coeffs[:,j] = 1 / a   
                
                # cross product in spherical coordinates
                u, v = uv_cross_product(x_coeffs * grad_x,
                                        y_coeffs * grad_y,
                                        lats,
                                        lons)
            
                u_arr[:,:,lvl_i] = u
                v_arr[:,:,lvl_i] = v 
                                
            # divergence term
            for lvl_i in range(l):
                even = xi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                # change m,n arangement to that used by shtools
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH2',extend=False)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
                    
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
                
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_y))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    x_coeffs[:,j] = nui / a
                    y_coeffs[:,j] = 1 / a 
                
                u = x_coeffs * (grad_x)
                v = y_coeffs * grad_y
                            
                u_arr[:,:,lvl_i] = u_arr[:,:,lvl_i] + u 
                v_arr[:,:,lvl_i] = v_arr[:,:,lvl_i] + v
                                
            # temperature to grid routine
            for lvl_i in range(l):
                even = jpot[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)    # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH2',extend=False)
                data = np.real(grid.to_array())
            
                t_arr[:,:,lvl_i] = data
                
            # set dimensions to python standard [zdim,ydim,xdim], from y, u, z
            u_arr = np.swapaxes(u_arr,0,2) # z, u, y
            u_arr = np.swapaxes(u_arr,1,2) # z, y, u
            v_arr = np.swapaxes(v_arr,0,2)
            v_arr = np.swapaxes(v_arr,1,2)
            t_arr = np.swapaxes(t_arr,0,2)
            t_arr = np.swapaxes(t_arr,1,2)
            p_field = np.swapaxes(p_field,0,2)
            p_field = np.swapaxes(p_field,1,2)
            
            # set geopotential dimensions for python standard
            geop = np.swapaxes(geop[:,::-1],0,1)
            
            # # # calculate geometric height using t,p fields
            # z_arr = pressure_to_altitude_topography(p_field,
            #                                         t_arr,
            #                                         geop,
            #                                         lat_sigprim)
            
            print('height km: ', (np.mean(z_arr,axis=(1,2)) / 1000)[8])
            # print('hPa   : ', (np.mean(p_field,axis=(1,2))/100)[38])
            
            # z_arr = pressure_to_altitude_scaleheight(p_field,
            #                                          H,
            #                                          geop,
            #                                          lat_sigprim)
            
            # print('height km: ', np.mean(z_arr,axis=(1,2)) / 1000)
            # print('hPa   : ', np.mean(p_field,axis=(1,2))/100)
            
            # print('height km: ', (np.mean(z_arr,axis=(1,2)) / 1000)[30])
            # print('hPa   : ', (np.mean(p_field,axis=(1,2))/100)[30])
            
            z_arr = hydro_ideal(p_field,
                                t_arr,
                                geop,
                                lat_sigprim)
            
            # print('height km: ', np.mean(z_arr,axis=(1,2)) / 1000)
            # print('hPa   : ', np.mean(p_field,axis=(1,2))/100)
            
            # print('height km: ', (np.mean(z_arr,axis=(1,2)) / 1000)[30])
            # print('hPa   : ', (np.mean(p_field,axis=(1,2))/100)[30])
                        
            # interpolate to input levels (Pa)
            u_arr = interplevel(u_arr,z_arr,interp_levels)
            v_arr = interplevel(v_arr,z_arr,interp_levels)
            t_arr = interplevel(t_arr,z_arr,interp_levels)
                
            print('Reading record: ', t_steps)
            full_u[t_steps,:,:,:] = u_arr
            full_v[t_steps,:,:,:] = v_arr
            full_t[t_steps,:,:,:] = t_arr
            
        else:
            # y, x, z dim arrays to hold read in fields
            u_arr = np.zeros((nmax_grid * 2,nmax_grid * 2 * 2, l))
            v_arr = np.zeros((nmax_grid * 2,nmax_grid * 2 * 2, l))
            t_arr = np.zeros((nmax_grid * 2,nmax_grid * 2 * 2, l))
            
            # read in spherical harmonic coefficients (shc)
            jvor = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            jdiv = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            jpot = f.read_record(dtype=np.complex64).reshape(m+1,n+1,l,order='F')
            pssp = f.read_record(dtype=np.complex64).reshape(m+1,n+1,order='F')
            geop = f.read_record(dtype=np.float32).reshape(k1,k2,order='F') 
            
            # read in surface pressure field and construct odd shcs
            even = pssp[:,:]      
            odd = np.copy(even)   
            
            for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
            coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                              dtype=np.complex64)
            
            # change [m,n] arangement to [n,m] used by shtools
            coeffs[0,:,:] = even
            coeffs[1,1:,:] = odd[1:,:]
            ct = np.swapaxes(coeffs,1,2)
            
            # zero pad to square dim, if necessary. nmax_grid >= m + 1
            shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
            shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
            
            # Condon-Shortley phase is used (-1) or not (1)
            x = sh.SHCoeffs.from_array(shs,
                                       normalization='4pi',
                                       csphase=1)
            
            # coordinate system exp transform to get surface pressure in Pa
            grid = x.expand(grid='DH2',extend=False)
            ps = np.exp(np.real(grid.to_array())) 
            p_field = np.zeros((nmax_grid * 2,nmax_grid * 2 * 2, l))
            
            for level in range(l):
                # pressure coordinates of (u,v)
                p_field[:,:,level] = ps * slv[level]
                                
            # u,v field from vorticity and divergence
            psi = np.copy(jvor)
            xi = np.copy(jdiv)
            
            # Del operator x-y coeffs from numerical recipe gauleg(ymin,ymax,ydim)
            mu = gauleg(-1.,1.,nmax_grid * 2)
            nui = 1 / np.sqrt(1 - mu**2) 
                
            for lvl in range(l):
                for n_sh in range(1,n + 1):
                    # n_sh runs from 1 to (n + 1); jvor(n=0) and jdiv(n=0) are (0,0j)
                    eigenvalue_leg = -n_sh * (n_sh + 1) / (a ** 2) # Saravanan p. 17
                    psi[:,n_sh,lvl] = 1 / eigenvalue_leg * jvor[:,n_sh,lvl]
                    xi[:,n_sh,lvl] = 1 / eigenvalue_leg * jdiv[:,n_sh,lvl]
                            
            # vorticity term
            for lvl_i in range(l):
                even = psi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH2',extend=False)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
        
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
        
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_x))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    # Saravanan recipe
                    x_coeffs[:,j] = nui / a 
                    y_coeffs[:,j] = 1 / a   
                
                # cross product in spherical coordinates
                u, v = uv_cross_product(x_coeffs * grad_x,
                                        y_coeffs * grad_y,
                                        lats,
                                        lons)
            
                u_arr[:,:,lvl_i] = u
                v_arr[:,:,lvl_i] = v
                            
            # divergence term
            for lvl_i in range(l):
                even = xi[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)   # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                # change m,n arangement to that used by shtools
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH2',extend=False)
                data = np.real(grid.to_array())
                lats = grid.lats()
                lons = grid.lons()
                    
                dphi = latfac * (lats[1] - lats[0]) / 180 * np.pi
                dlam = (lons[1] - lons[0]) / 180 * np.pi
                
                dim_x = np.shape(data)[1]
                           
                grad_y, grad_x = np.gradient(data,dphi,dlam)
                x_coeffs = np.zeros(np.shape(grad_y))
                y_coeffs = np.zeros(np.shape(grad_y))
                    
                for j in range(dim_x):
                    x_coeffs[:,j] = nui / a
                    y_coeffs[:,j] = 1 / a 
                
                u = x_coeffs * (grad_x)
                v = y_coeffs * grad_y
            
                u_arr[:,:,lvl_i] = u_arr[:,:,lvl_i] + u 
                v_arr[:,:,lvl_i] = v_arr[:,:,lvl_i] + v
                
            # temperature to grid routine
            for lvl_i in range(l):
                even = jpot[:,:,lvl_i] # p_{m,n} coefficients
                odd = np.copy(even)    # p_{-m,n}
                
                for m_i in range(m + 1):
                    odd[m_i,:] = (-1)**m_i * np.conj(even[m_i,:])
                
                coeffs = np.zeros((2,np.shape(even)[0],np.shape(even)[1]),
                                  dtype=np.complex64)
            
                coeffs[0,:,:] = even
                coeffs[1,1:,:] = odd[1:,:]
                
                # change m,n arangement to that used by shtools
                ct = np.swapaxes(coeffs,1,2)
                
                # zero pad to square dim. nmax_grid >= m + 1
                shs = np.zeros((2,nmax_grid,nmax_grid),dtype=np.complex64)
                shs[:,:nmax_grid,:(m + 1)] = ct[:,:nmax_grid,:]
                
                # Condon-Shortley phase is used (-1) or not (1)
                x = sh.SHCoeffs.from_array(shs,
                                           normalization='4pi',
                                           csphase=1)
                
                grid = x.expand(grid='DH2',extend=False)
                data = np.real(grid.to_array())
            
                t_arr[:,:,lvl_i] = data
                
            # set dimensions to python standard [zdim,ydim,xdim], from y, u, z
            u_arr = np.swapaxes(u_arr,0,2) # z, u, y
            u_arr = np.swapaxes(u_arr,1,2) # z, y, u
            v_arr = np.swapaxes(v_arr,0,2)
            v_arr = np.swapaxes(v_arr,1,2)
            t_arr = np.swapaxes(t_arr,0,2)
            t_arr = np.swapaxes(t_arr,1,2)
            p_field = np.swapaxes(p_field,0,2)
            p_field = np.swapaxes(p_field,1,2)
            
            # set geopotential dimensions for python standard
            geop = np.swapaxes(geop[:,::-1],0,1)
                        
            # calculate geometric height using t,p fields
            z_arr = pressure_to_altitude_topography(p_field,
                                                    t_arr,
                                                    geop,
                                                    lat_sigprim)
            
            # z_arr = hydro_ideal(p_field,
            #                     t_arr,
            #                     geop,
            #                     lat_sigprim)
                        
            u_arr = interplevel(u_arr,z_arr,interp_levels)
            v_arr = interplevel(v_arr,z_arr,interp_levels)
            t_arr = interplevel(t_arr,z_arr,interp_levels)
                
            print('Reading record: ', t_steps)
            full_u[t_steps,:,:,:] = u_arr
            full_v[t_steps,:,:,:] = v_arr
            full_t[t_steps,:,:,:] = t_arr
                    
    f.close()
        
    return full_u, full_v, full_t


"""
###############################################################################
# INTERPOLATION/GRID FUNCTIONS                                                #
###############################################################################

"""

def zm_interp_lat_grid(input_zonal_wind,
                       input_temperature,
                       input_y_grid,
                       output_y_grid):
    """
    Interpolate zonal mean fields to different latitude grid (radians).    
    input and output dimensions: [time,height,latitude]
    """
    
    int_array_u = np.zeros((np.shape(input_zonal_wind)[0],
                               np.shape(input_zonal_wind)[1],
                               np.shape(output_y_grid)[0]))
    int_array_t = np.zeros((np.shape(input_zonal_wind)[0],
                               np.shape(input_zonal_wind)[1],
                               np.shape(output_y_grid)[0]))

    t_steps = np.shape(input_zonal_wind)[0]
    z_steps = np.shape(input_zonal_wind)[1]

    for t_ind in range(t_steps):
        for z_ind in range(z_steps):
            int_array_u[t_ind,z_ind,:] = np.interp(output_y_grid,
                                                   input_y_grid,
                                                   input_zonal_wind[t_ind,
                                                                    z_ind,:])
            int_array_t[t_ind,z_ind,:] = np.interp(output_y_grid,
                                                   input_y_grid,
                                                   input_temperature[t_ind,
                                                                     z_ind,:])
            
    return int_array_u, int_array_t


def interp_3D_alt_grid(u_array,
                       v_array,
                       t_array,
                       p_array,
                       new_zlvls):
    """
    Interpolate [time,alt,lat,lon] u,v,t grids to altitude levels
    set by new_zlvls.
    """
    
    tdim, zdim, ydim, xdim = np.shape(u_array)
    
    u_new = np.zeros((tdim,np.shape(new_zlvls)[0],ydim,xdim))
    v_new = np.zeros((tdim,np.shape(new_zlvls)[0],ydim,xdim))
    t_new = np.zeros((tdim,np.shape(new_zlvls)[0],ydim,xdim))
                
    for t in range(tdim):
        # function assumes [time, altitude, latitude, longitude] dims
        u_new[t,:,:,:] = interplevel(u_array[t,:,:,:],p_array[t,:,:,:],
             new_zlvls)
        v_new[t,:,:,:] = interplevel(v_array[t,:,:,:],p_array[t,:,:,:],
             new_zlvls)
        t_new[t,:,:,:] = interplevel(t_array[t,:,:,:],p_array[t,:,:,:],
             new_zlvls)

    return u_new, v_new, t_new


def sigprim_plvls():
    """     
    Function call returns list of Sigprim model -- initialized -- levels (Pa).
    
    """
    
    lvls = np.array([3.57821e-006,7.14165e-006,1.39804e-005,2.65079e-005,
                 4.79696e-005,8.19642e-005,0.000132236,0.000203615,0.000303671,
                 0.000444207,0.000642736,0.000924545,0.00132580,0.00221539,
                 0.00316592,0.00452402,0.00646411,0.00923494,0.0131908,
                 0.0188355,0.0268835,0.0383447,0.0546378,0.0777402,0.110373,
                 0.156212,0.226033,0.323055,0.461724,0.659914,0.943176,
                 1.34802,1.92665,2.75364,3.93562,5.62495,8.03939,11.4902,
                 16.4223,23.4714,30.0000,40.0000,50.0000,70.0000,100.000,
                 200.000,300.000,400.000,500.000,700.000,1000.00,2000.00,
                 3000.00,4000.00,5000.00,7000.00,10000.0,15000.0,20000.0,
                 25000.0,30000.0,35000.0,40000.0,45000.0,50000.0,55000.0,
                 60000.0,65000.0,70000.0,72500.0,75000.0,77500.0,80000.0,
                 82500.0,85000.0,87500.0,90000.0,92500.0,95000.0,
                 97500.0,100000.])
    
    return lvls


def sigprim_plvls_era5():
    """     
    Function call returns list of Sigprim model - initialized - ERA_Pgrid lvls
    
    """
    
    lvls = np.array([3.57821e-006,7.14165e-006,1.39804e-005,2.65079e-005,
                 4.79696e-005,8.19642e-005,0.000132236,0.000203615,0.000303671,
                 0.000444207,0.000642736,0.000924545,0.00132580,0.00221539,
                 0.00316592,0.00452402,0.00646411,0.00923494,0.0131908,
                 0.0188355,0.0268835,0.0383447,0.0546378,0.0777402,0.110373,
                 0.156212,0.226033,0.323055,0.461724,0.659914,1.000,2.5500,
                 3.88000,5.75000,8.29000,11.68000,16.11000,21.8000,28.9900,
                 37.9300,48.9200,62.24000,78.2100,97.1600,119.420,175.3100,
                 292.9800,398.2900,528.5100,686.540,981.770,1996.810,3017.7600,
                 4053.3100,4959.5200,7111.870,9841.640,14790.580,19736.79,
                 24836.340,29651.55,35198.87,39851.16,45058.58,50750.21,
                 54803.120,61066.46,65244.24,69330.43,73253.25,75134.26,
                 78705.28,80386.22,81993.02,84976.68,87649.57,90016.69,92091.93,
                 94702.40,97918.52,101204.94])
    
    return lvls


def sigprim_plvls_era5_hr():
    """     
    Function call returns list of Sigprim model - initialized - ERA_Pgrid lvls
    
    """
    
    lvls = np.load('/home/wim/Desktop/Projects/Model/ERA5/ERA_Pgrid_hr.npy')
    
    return lvls


"""
###############################################################################
# WACCMX-DART / NAVGEM-HA / ERA5 FUNCTIONS                                    #
###############################################################################

"""

def WACCM_file_creator(var_i,
                       y_step,
                       x_step,
                       folder):
    """
    Function returns WACCMX-DART variable var_i.
    Variables: (U, V, T, Z) with corresponding var_i (0, 1, 2, 3)     
    
    Folder (path to files) input e.g. '/media/wim/NieuwVolume/WACCMX-DART'
    """
        
    # file names of WACCMX-DART data files
    periods = ['2009011000-2009011923','2009012000-2009012923',
               '2009013000-2009013123','2009020100-2009020923',
               '2009021000-2009021923','2009022000-2009022823',
               '2009030100-2009030923','2009031000-2009031923']
    
    varibs = ['U','V','T','Z']
    var_codes = ['U','V','T','GPH']
    
    for index, item in enumerate(periods):
        
        file_n = 'WACCMX+DART_' + varibs[var_i] + '_' + item
        
        # time, pressure, latitude, longitude
        data = Dataset(folder + file_n + '.nc',mode='r')
        data_parameter = data[var_codes[var_i]][:,:,::y_step,::x_step]
        
        time = data['YYYYMMDDHH'][:]
        lats = data['LATITUDE'][::y_step]
        lons = data['LONGITUDE'][::x_step]
        plvls = data['PRESSURE'][:]
               
        # concatenate the files together w.r.t. time
        if index == 0:
            conc_file = data_parameter
            conc_time = time
        else:
            conc_file = np.concatenate((conc_file,data_parameter),axis=0)
            conc_time = np.concatenate((conc_time,time),axis=0)
            
    return conc_file, conc_time, lats, lons, plvls


def WACCMX_DART_hrly_to_daily_mean(field,
                                   day_0):
    """
    output day_0 = 11 for 2009 WACCMX-DART
    """
    
    # mean days are gonna start coming after 11th of Jan 2009
    
    tdim, zdim, ydim, xdim = np.shape(field)
    t_steps = int(tdim / 24 - 2)
    
    mean = np.empty((t_steps,zdim,ydim,xdim))
    time_arr = np.empty(t_steps,dtype=np.float32)
    
    for t in range(t_steps):
        mean[t,:,:,:] = np.mean(field[(int(24*t)):(int(24*t+48)),:,:,:],axis=0)
        time_arr[t] = day_0 + t

    return mean, time_arr


"""
###############################################################################
# HWM14 + NRLMSISE-00 WIND AND TEMPERATURE ON ALTITUDE/PRESSURE GRID          #
###############################################################################

"""

def hwm_msise_uvtp(date,
                   alt,
                   lats,
                   lons,
                   folder):
    
    """
    Function reads f107a, f107 and ap from file based on DOY. Calculates
    pressure from MSISE number densities, temperatures, and the ideal gas law. 
    
    Input:
        - date (datetime object; year, month, day)
        - altitude (km)
        - latitude array (degrees)  (e.g. np.arange(-90,91))
        - longitude array (degrees) (e.g. np.arange(0,360))
        - folder; path to .txt files containing geophysical indices 
    
    Output:
        - lat/lon NRLMSISE-00 temperature and pressure 
        - lat/lon HWM14 zonal and meridional winds
        
    f107a: float
		The observed f107a (81-day running mean of f107) centred at date.
	f107: float
		The observed f107 value on the previous day.
	ap: float
		The ap value at date.
        
    NAVGEM-HA top = 0.1435314 Pa (~ 90 km)
    DOPE top      = 3.57821e-6 Pa
    
    n_levels DOPE between NAVGEM-HA and DOPE top: 32
    """
    
    doy = (date - datetime(2008,1,1,0)).days
    
    # read in f10.7 from file based on input date
    if alt <= 80:
        f107a, f107, ap = (150,150,4)
    else:
        # f10.7 and (derived) f10.7A
        date_f = float(date.strftime('%Y%m%d'))
        file = np.loadtxt(folder + 'f107_08_16.txt')
        dates, vals = file[:,0], file[:,1]
                
        ind = 0
        found = True
        while found:
            if dates[ind] == date_f:
                found = False
            else:
                ind += 1
                
        f107 = vals[int(ind - 1)] # day before today / doy
        f107a = np.mean(vals[int(ind - 40):int(ind + 41)]) # 81 day mean
                
        # ap index
        file_ap = open(folder + 'ap_08_16.txt','r')

        for m in range(doy):
            line = file_ap.readline()
            ap = float(line[53:55])

        file_ap.close()
            
    # index 5 = density, index -1 = temperature
    msise = msise_flat(date,
                       alt,
                       lats[:,None],
                       lons[None,:],
                       f107a,
                       f107,
                       ap)
        
    r = 8.314               # J / K / mol
    avg = 6.022140e23       # avocado's number
    
    # number density m**-3
    he = 1e6 * msise[:,:,0] 
    o = 1e6 * msise[:,:,1]  
    n2 = 1e6 * msise[:,:,2] 
    o2 = 1e6 * msise[:,:,3] 
    ar = 1e6 * msise[:,:,4] 
    h = 1e6 * msise[:,:,6]  
    n = 1e6 * msise[:,:,7] 
    ox = 1e6 * msise[:,:,8] 
    
    n_density = he + o + n2 + o2 + ar + h + n + ox
    
    temp = msise[:,:,-1]    # [K]
    
    p = n_density / avg * r * temp # pressure using ideal gas law
        
    # horizontal winds fro HWM14. Does not use f107 and f107a
    winds = HWM142D(alt=alt, glatlim=[lats[0], lats[-1]],
                    glatstp=(lats[1] - lats[0]), glonlim=[lons[0], lons[-1]], 
                    glonstp=(lons[1] - lons[0]), option=6,verbose=False,
                    day=int(date.strftime("%j")), 
                    ut=int(date.strftime('%H')),
                    year=int(date.strftime('%Y')),
                    ap=[-1,ap])
    
    u = winds.Uwind
    v = winds.Vwind
    
    return u, v, temp, p


def hwm_msise_uvtp_climatology(date,
                               alt,
                               lats,
                               lons,
                               folder):
    """
    Function reads f107a, f107 and ap from file based on DOY. Calculates
    pressure from MSISE number densities, temperatures, and the ideal gas law. 
    
    Input:
        - date (datetime object; year, month, day)
        - altitude (km)
        - latitude array (degrees)  (e.g. np.arange(-90,91))
        - longitude array (degrees) (e.g. np.arange(0,360))
        - folder; path to .txt files containing geophysical indices 
    
    Output:
        - lat/lon NRLMSISE-00 temperature and pressure 
        - lat/lon HWM14 zonal and meridional winds
        
    f107a: float
		The observed f107a (81-day running mean of f107) centred at date.
	f107: float
		The observed f107 value on the previous day.
	ap: float
		The ap value at date.
        
    NAVGEM-HA top = 0.1435314 Pa (~ 90 km)
    DOPE top      = 3.57821e-6 Pa
    
    n_levels DOPE between NAVGEM-HA and DOPE top: 32
    
    """
        
    # standard geomagnetic parameters
    f107a, f107, ap = (150,150,4)
               
    # index 5 = density, index -1 = temperature
    msise = msise_flat(date,
                       alt,
                       lats[:,None],
                       lons[None,:],
                       f107a,
                       f107,
                       ap)
    
    r = 8.314               # J / K / mol
    avg = 6.022140e23       # avocado's number
    
    # number density m**-3
    he = 1e6 * msise[:,:,0] 
    o = 1e6 * msise[:,:,1]  
    n2 = 1e6 * msise[:,:,2] 
    o2 = 1e6 * msise[:,:,3] 
    ar = 1e6 * msise[:,:,4] 
    h = 1e6 * msise[:,:,6]  
    n = 1e6 * msise[:,:,7] 
    ox = 1e6 * msise[:,:,8] 
    
    n_density = he + o + n2 + o2 + ar + h + n + ox
    
    temp = msise[:,:,-1] # [K]
    
    p = n_density / avg * r * temp # pressure using ideal gas law
    
    # horizontal winds fro HWM14. Does not use f107 and f107a
    winds = HWM142D(alt=alt, glatlim=[lats[0], lats[-1]],
                    glatstp=(lats[1] - lats[0]), glonlim=[lons[0], lons[-1]], 
                    glonstp=(lons[1] - lons[0]), option=6,verbose=False,
                    day=int(date.strftime("%j")), 
                    ut=int(date.strftime('%H')),
                    year=int(date.strftime('%Y')),
                    ap=[-1,ap])
    
    u = winds.Uwind
    v = winds.Vwind
    
    return u, v, temp, p


def era5_model_levels():
    
    """
    Returns pressure levels (Pa) of ERA5 model levels from 2014 and 2015
    heating fields; soon to be obsolete
    
    """
    
    lvls = [0.01,0.0388,0.0829,0.1611,0.2899,0.4892,0.7821,1.1942,1.7531,
            2.4875,3.4270,4.6010,6.0388,7.7686,9.8177,12.2123,14.9770,18.1348,
            21.7076,26.7735,30.1776,35.1111,40.5321,46.4498,54.5299,71.1187,
            93.3527,121.0526,155.3448,197.3679,248.3634,309.6684,351.9887,
            398.8516,450.5858,507.5021,548.0312,589.6797,631.6194,673.0352,
            713.1631,751.3426,787.0528,819.9302,835.2358,849.7668,863.5190,
            876.4957,888.7066,900.1669,920.9193,938.9532,954.5059,967.8315,
            979.1852,988.8133,996.9452,1003.7906]
    
    lvls = np.array(lvls) * 1e2
    
    return lvls


def merra2_model_levels():
    
    """
    Returns pressure levels (Pa) of MERRA2 model levels 
    
    """
    
    lvls = [0.0100, 0.0200, 0.0327, 0.0476, 0.0660, 0.0893, 0.1197, 0.1595,
            0.2113, 0.2785, 0.3650, 0.4758, 0.6168, 0.7951, 1.0194, 1.3005,
            1.6508, 2.0850, 2.6202, 3.2764, 4.0766, 5.0468, 6.2168, 7.6198,
            9.2929, 11.2769, 13.6434, 16.4571, 19.7916, 23.7304, 28.3678,
            33.8100, 40.1754, 47.6439, 56.3879, 66.6034, 78.5123, 92.3657,
            108.6630, 127.8370, 150.3930, 176.9300, 208.1520, 244.8750,
            288.0830, 337.5000, 375.0000, 412.5000, 450.0000, 487.5000,
            525.0000, 562.5000, 600.0000, 637.5000, 675.0000, 700.0000,
            725.0000, 750.0000, 775.0000, 800.0000, 820.0000, 835.0000,
            850.0000, 865.0000, 880.0000, 895.0000, 910.0000, 925.0000,
            940.0000, 955.0000, 970.0000, 985.0000]
    
    lvls = np.array(lvls) * 1e2
    
    return lvls


def merra2_pressure_levels():
    
    """
    Returns pressure levels (Pa) of MERRA2 pressure levels 
    
    """
    
    lvls = [0.0100, 0.0200, 0.0327, 0.0476, 0.0660, 0.0893, 0.1197, 0.1595,
            0.2113, 0.2785, 0.3650, 0.4758, 0.6168, 0.7951, 1.0194, 1.3005,
            1.6508, 2.0850, 2.6202, 3.2764, 4.0766, 5.0468, 6.2168, 7.6198,
            9.2929, 11.2769, 13.6434, 16.4571, 19.7916, 23.7304, 28.3678,
            33.8100, 40.1754, 47.6439, 56.3879, 66.6034, 78.5123, 92.3657,
            108.6630, 127.8370, 150.3930, 176.9300, 208.1520, 244.8750,
            288.0830, 337.5000, 375.0000, 412.5000, 450.0000, 487.5000,
            525.0000, 562.5000, 600.0000, 637.5000, 675.0000, 700.0000,
            725.0000, 750.0000, 775.0000, 800.0000, 820.0000, 835.0000,
            850.0000, 865.0000, 880.0000, 895.0000, 910.0000, 925.0000,
            940.0000, 955.0000, 970.0000, 985.0000]
    
    lvls = np.array(lvls) * 1e2
    
    return lvls


def merra2_plevels():
    """
    Returns MERRA2 pressure levels in Pa
    
    """

    merra_lvls = [100000,97500,95000,92500,90000,87500,85000,82500,80000,
              77500,75000,72500,70000,65000,60000,55000,50000,45000,40000,
              35000,30000,25000,20000,15000,10000,7000,5000,4000,3000,2000,
              1000,700,500,400,300,200,100,70,50,40,30,10]
    
    merra_lvls = np.array(merra_lvls)
    
    return merra_lvls


def NAVGEM_HA_plvls():
    """
    Returns NAVGEM-HA pressure levels in Pa
    
    """
    nav_lvls = [0.0178022, 0.040599, 0.0797394, 0.143531, 0.244708, 0.401782,
                0.640952, 0.99877, 1.52532, 2.2883, 3.3788, 4.91932, 7.07405,
                10.0608, 14.1641, 19.7495, 27.2815, 37.3509, 50.7206, 68.3903,
                91.6937, 122.436, 163.072, 216.958, 288.624, 384.057, 511.047,
                679.719, 903.143, 1197.6, 1582.16, 2077.65, 2705.44, 3485.69, 
                4436.77, 5576.92, 6926.21, 8506.79]
    
    return nav_lvls

# new internal sigprim functions 
def uvt_read_altitude_grid(file,
                           folder,
                           day_begin,
                           day_final,
                           model_dt,
                           interp_levels=[100000]):

    model_tstep = int(24 / model_dt)         
    
    f = FortranFile(folder + file,'r')

    # read in grid dimensions / metadata gunk
    md = f.read_record(dtype=np.int32)
    (k1, k2, m, n, l) = (md[0], md[1], md[2], md[3], md[4])
    
    slv = f.read_record(dtype=np.float32)
    lat = f.read_record(dtype=np.float32)
    lon = np.linspace(0,360,k1,endpoint=False)
    zstdlvl = f.read_record(dtype=np.float32)
    
    n_steps = day_final * model_tstep
    
    full_u = np.zeros((n_steps,np.shape(interp_levels)[0],k2,k1))
    full_v = np.zeros((n_steps,np.shape(interp_levels)[0],k2,k1))
    full_t = np.zeros((n_steps,np.shape(interp_levels)[0],k2,k1))
    
    for k in range(day_begin * model_tstep):
        # skip the first day_begin number of days (spin-up)
        u = f.read_record(dtype=np.float32).reshape(k1,k2,l,order='F')
        v = f.read_record(dtype=np.float32).reshape(k1,k2,l,order='F')
        t = f.read_record(dtype=np.float32).reshape(k1,k2,l,order='F')
        surf_p = f.read_record(dtype=np.float32).reshape(k1,k2,order='F')
        geob = f.read_record(dtype=np.float32).reshape(k1,k2,order='F')
        
    for t_steps in range(n_steps):
        # read in and interpolate 
        u = f.read_record(dtype=np.float32).reshape(k1,k2,l,order='F')
        v = f.read_record(dtype=np.float32).reshape(k1,k2,l,order='F')
        t = f.read_record(dtype=np.float32).reshape(k1,k2,l,order='F')
        surf_p = f.read_record(dtype=np.float32).reshape(k1,k2,order='F')
        geob = f.read_record(dtype=np.float32).reshape(k1,k2,order='F')
                    
        # change from [x,y,z] to python standard [z,y,x]
        u = np.swapaxes(u,0,2)
        v = np.swapaxes(v,0,2)
        t = np.swapaxes(t,0,2)
        surf_p = np.exp(np.swapaxes(surf_p,0,1)) # log pressure sigma
        geob = np.swapaxes(geob,0,1)
            
        dy, dx = np.shape(surf_p)
        p = np.zeros((l,dy,dx))
                    
        for lev in range(l):
            p[lev,:,:] = surf_p * slv[lev]
                    
        # calculate geometric height using t,p fields
        z_arr = pressure_to_altitude_topography_grid(p,
                                                     t,
                                                     geob)
        
        # plt.figure()
        # plt.plot(np.mean(z_arr,axis=(1,2)) / 1000)
                                            
        # interpolate to input levels (Pa)
        u_arr = interplevel(u,z_arr,interp_levels)
        v_arr = interplevel(v,z_arr,interp_levels)
        t_arr = interplevel(t,z_arr,interp_levels)
                            
        print('Reading record: ', t_steps)
        full_u[t_steps,:,:,:] = u_arr
        full_v[t_steps,:,:,:] = v_arr
        full_t[t_steps,:,:,:] = t_arr
        
    f.close()
        
    return full_u, full_v, full_t, lat, lon



def pressure_to_altitude_topography_grid(p,
                                         t,
                                         geopotential):
    
    """
    Numerically integrate barometric equation.
    
    - takes into account gravity and molar mass as function of altitude
    - molar mass profile based on global mean from MSISE-00 
    
    inputs: 
    
        p : pressure (Pa)
        t : temperature (K)
    
    """
            
    g_0 = 9.80665   # m s**-2
    r_d = 8.3144598 # J mol**-1 K**-1
    r_e = 6.3781e6  # m
    
    # array to hold geometric altitude values
    z = np.zeros(np.shape(p))
    z[-1,:,:] = geopotential / g_0
        
    path_f = '/home/wim/Desktop/Projects/Model/Data/'
    molar_mass = np.load(path_f + 'molar_mass.npy')
    mm_heights = np.load(path_f + 'molar_mass_height_m.npy')
          
    dz_dim, dy_dim, dx_dim = np.shape(p)
                    
    for k in np.arange(dz_dim - 1,0,-1):
        # height adjusted gravitational acceleration assuming spherical earth
        g = g_0 * (r_e / (r_e + z[k,:,:])) ** 2
        
        # molar mass from global mean height-profile
        mean_height = np.mean(z[k,:,:])
        func_interp = interpolate.interp1d(mm_heights,molar_mass)

        if mean_height < 0:
            m = func_interp(0)
        else:
            m = func_interp(mean_height)
            
        # forward difference at bottom boundary
        if k == (dz_dim - 1):
            dz = -(p[k-1,:,:] - p[k,:,:]) / (p[k,:,:] * m * g) * r_d * t[k,:,:]
            
        # central difference otherwise
        else:
            dz = (-0.5 * ((p[k-1,:,:] - p[k+1,:,:])) / 
                  (p[k,:,:] * m * g) * r_d * t[k,:,:])
        
        # increment altitude with dz
        z[k-1,:,:] = z[k,:,:] + dz
            
    return z


def read_uv_grid(file,
                 folder,
                 rec_n):     
    
    f = FortranFile(folder + file,'r')

    # read in grid dimensions / metadata gunk
    md = f.read_record(dtype=np.int32)
    (k1, k2, m, n, l) = (md[0], md[1], md[2], md[3], md[4])
    
    slv = f.read_record(dtype=np.float32)
    lat = f.read_record(dtype=np.float32)
    lon = np.linspace(0,360,k1,endpoint=False)
    zstdlvl = f.read_record(dtype=np.float32)
    
    # np.save('/home/wim/Desktop/Projects/Model/vertical_grids/sigprim_lvls_ERAn.npy',slv)
        
    for k in range(rec_n):
        # skip the first day_begin number of days (spin-up)
        u = f.read_record(dtype=np.float32).reshape(k1,k2,l,order='F')
        v = f.read_record(dtype=np.float32).reshape(k1,k2,l,order='F')
        t = f.read_record(dtype=np.float32).reshape(k1,k2,l,order='F')
        surf_p = f.read_record(dtype=np.float32).reshape(k1,k2,order='F')
        geob = f.read_record(dtype=np.float32).reshape(k1,k2,order='F')
        
    # change from [x,y,z] to python standard [z,y,x]
    u = np.swapaxes(u,0,2)
    v = np.swapaxes(v,0,2)
    t = np.swapaxes(t,0,2)
    surf_p = np.exp(np.swapaxes(surf_p,0,1))
    geob = np.swapaxes(geob,0,1)
            
    f.close()
        
    return u, v, t, lat, lon, surf_p, geob

#%%
