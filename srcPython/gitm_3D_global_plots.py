#!/usr/bin/env python
#-----------------------------------------------------------------------------
# gitm_3D_global_plots
#
# Author: Angeline G. Burrell, UMichigan, Jan 2013
#
# Comments: Routine to make 3D color plots of thermospheric and ionospheric
#           quantities that are either at a single altitude or are constant
#           with altitude as functions of geographic latitude and longitude.
#           Map contours are also availabel if you have the mpl_toolkit
#           baseline installed (available on fink, macports, or online at
#           http://matplotlib.org/basemap).
#
# Includes: plot_single_3D_image          - plots a single rectangular or polar
#                                           alt slice
#           plot_single_nsglobal_3D_image - plots northern and southern polar
#                                           altitude slice
#           plot_global_3D_snapshot       - plots northern and southern polar
#                                           and a midlatitude rectangular
#                                           altitude slice
#           plot_mult_3D_slices           - plot multiple polar or rectangular
#                                           altitude slices
#           -----------------------------------------------------------------
#           plot_rectangular_3D_global    - plot a rectangular geographic
#                                           contour with or without a map
#           plot_polar_3D_global          - plot a polar geographic contour
#                                           with or without a map
#           -----------------------------------------------------------------
#           add_colorbar                  - add a colorbar to a contour plot
#           center_polar_cap              - center radial coordinates for a
#                                           polar plot
#           find_zlimits                  - find the upper and lower limits
#                                           for a list of GITM data arrays
#----------------------------------------------------------------------------

'''
Plot data from a 2D GITM file (or 3D GITM file at a single altitude) for
different geographic configurations
'''

# Import modules
from mpl_toolkits.basemap import Basemap
import sys
import string 
import math
import numpy as np
from spacepy.pybats import gitm
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.ticker import ScalarFormatter,FormatStrFormatter,MultipleLocator

def plot_single_3D_image(plot_type, zkey, gData, title=None, figname=None,
                         aindex=-1, nlat=90, slat=-90, linc=6, earth=False,
                         *args, **kwargs):
    '''
    Creates a rectangular or polar map projection plot for a specified latitude
    range.
    Input: plot_type = key to determine plot type (rectangular, polar)
           zkey      = key for z variable (ie 'Vertical TEC')
           gData     = gitm bin structure
           title     = plot title
           figname   = file name to save figure as (default is none)
           aindex    = altitude index (default -1 if it is a 2D parameter)
           nlat      = northern latitude limit (degrees North, default 90)
           slat      = southern latitude limit (degrees North, defalut -90)
           linc      = number of latitude tick incriments (default 6)
           earth     = include continent outlines for Earth (default False)
    '''

    # Initialize the z variable limits
    zmin, zmax = find_zlimits([gData], zkey, aindex, 6)

    # Initialize the new figure

    gf = True
    pf = False

    if(string.lower(plot_type)=="polar" and not earth):
        pf = True

    f  = plt.figure()
    ax = f.add_subplot(111, polar=pf)

    if(string.lower(plot_type)=="rectangular"):
        plot_rectangular_3D_global(ax, zkey, gData, zmin, zmax, 6, aindex, nlat,
                                   slat, linc, True, "r", title, "t", True,
                                   True, earth)
    elif(string.lower(plot_type)=="polar"):
        plot_polar_3D_global(ax, 1, zkey, gData, zmin, zmax, 6, aindex, nlat,
                             slat, linc, True, "r", title, "t", True, True,
                             earth)
    else:
        print "ERROR: unknown input type [", plot_type, "]\n"
        gf = False

    if gf:
        # Draw to screen.
        if plt.isinteractive():
            plt.draw() #In interactive mode, you just "draw".
        else:
            # W/o interactive mode, "show" stops the user from typing more 
            # at the terminal until plots are drawn.
            plt.show()

        # Save output file

        if figname is not None:
            plt.savefig(figname)

def plot_single_nsglobal_3D_image(zkey, gData, title=None, figname=None,
                                  aindex=-1, plat=90, elat=0, linc=3,
                                  earth=False, *args, **kwargs):
    '''
    Creates a figure with two polar map projections for the northern and 
    southern ends of a specified latitude range.
    Input: zkey      = key for z variable (ie 'Vertical TEC')
           gData     = gitm bin structure
           title     = plot title
           figname   = file name to save figure as (default is none)
           aindex    = altitude index (default -1 if it is a 2D parameter)
           plat      = polar latitude limit (degrees North, default +/-90)
           elat      = equatorial latitude limit (degrees North, defalut 0)
           linc      = number of latitude tick incriments (default 6)
           earth     = include Earth continent outlines (default False)
    '''

    # Initialize the z variable limits
    zmin, zmax = find_zlimits([gData], zkey, aindex, 6)

    # Initialize the new figure

    f  = plt.figure()
    pf = True
    nc = False
    sc = True

    if earth:
        pf = False
        nc = True
        sc = False

    if title:
        f.suptitle(title, size="medium")

    # Northern Plot
    axn = f.add_subplot(121, polar=pf)
    plot_polar_3D_global(axn, 2, zkey, gData, zmin, zmax, 6, aindex, plat,
                         elat, linc, False, "r", "North", "t", True, nc, earth)
    psn = list(axn.get_position().bounds)

    # Southern Plot
    plat = -plat
    if(elat != 0.0):
        elat = -elat

    axs = f.add_subplot(122, polar=pf)
    con = plot_polar_3D_global(axs, 2, zkey, gData, zmin, zmax, 6, aindex,
                               plat, elat, linc, True, "r", "South", "t",
                               True, sc, earth)
    pss = list(axs.get_position().bounds)

    # Adjust plot sizes
    pwidth = 0.5 * (psn[2] + pss[2])
    psn[0] = psn[0] + (pwidth - psn[2]) * 0.5
    psn[2] = pwidth
    pss[0] = pss[0] - (pwidth - pss[2])
    pss[2] = pwidth

    axs.set_position(pss)
    axn.set_position(psn)    

    # Draw to screen.
    if plt.isinteractive():
        plt.draw() #In interactive mode, you just "draw".
    else:
        # W/o interactive mode, "show" stops the user from typing more 
        # at the terminal until plots are drawn.
        plt.show()

    # Save output file

    if figname is not None:
        plt.savefig(figname)


def plot_global_3D_snapshot(zkey, gData, title=None, figname=None,
                            aindex=-1, earth=False, *args, **kwargs):
    '''
    Creates a map projection plot for the entire globe, seperating the polar
    and central latitude regions.
    Input: zkey      = key for z variable (ie 'Vertical TEC')
           gData     = gitm bin structure
           title     = plot title
           figname   = file name to save figure as (default is none)
           aindex    = altitude index (default -1 if it is a 2D parameter)
           nlat      = northern latitude limit (degrees North, default 90)
           slat      = southern latitude limit (degrees North, defalut 90)
           earth     = include Earth continent outlines (default False)
    '''

    # Initialize the z variable limits
    zmin, zmax = find_zlimits([gData], zkey, aindex, 6)

    # Initialize the new figure, starting with the mid- and low-latitudes

    f   = plt.figure()
    axl = f.add_subplot(212)
    plot_rectangular_3D_global(axl, zkey, gData, zmin, zmax, 6, aindex, 45.0,
                               -45.0, 6, False, "r", None, "t", True, True,
                               earth)
    psl = list(axl.get_position().bounds)

    # Add the North pole

    axn = f.add_subplot(221, polar=True)
    plot_polar_3D_global(axn, 1, zkey, gData, zmin, zmax, 6, aindex, 90, 45, 3,
                         False, "r", None, "t", True, False, earth)

    # Add the South pole

    axs     = f.add_subplot(222, polar=True)
    con, cb = plot_polar_3D_global(axs, 2, zkey, gData, zmin, zmax, 6, aindex,
                                   -90, -45, 3, False, "r", None, "t", True,
                                   True, earth)
    pss = list(axs.get_position().bounds)

    # Add a colorbar for the entire plot

    orient = 'vertical'
    if(cloc == 'l' or cloc == 'r'):
        orient = 'horizontal'

    cb = add_colorbar(con, zmin, zmax, 6, orient, gData[zkey].attrs['scale'],
                      gData[zkey].attrs['name'], gData[zkey].attrs['units'])
    bp = list(cb.ax.get_position().bounds)

    bp[0] = .9
    bp[1] = psl[1]
    bp[3] = pss[3] + .5
    bp[2] = bp[2] / 4.0 

    pss[0] = pss[0] - .08
    psl[0] = psl[0] - .02
    axs.set_position(pss)
    axl.set_position(psl)
    cb.ax.set_position(bp)

    if title:
        f.suptitle(title, size="medium")

    # Draw to screen.
    if plt.isinteractive():
        plt.draw() #In interactive mode, you just "draw".
    else:
        # W/o interactive mode, "show" stops the user from typing more 
        # at the terminal until plots are drawn.
        plt.show()

    # Save output file
        
    if figname is not None:
        plt.savefig(figname)

def plot_mult_3D_slices(plot_type, zkey, gData, aindex, title=None,
                        figname=None, nlat=90, slat=-90, linc=6, earth=False,
                        *args, **kwargs):
    '''
    Creates a rectangular or polar map projection plot for a specified latitude
    range.
    Input: plot_type = key to determine plot type (rectangular, polar)
           zkey      = key for z variable (ie 'Vertical TEC')
           gData     = gitm bin structure
           aindex    = list of altitude indices
           title     = plot title
           figname   = file name to save figure as (default is none)
           nlat      = northern latitude limit (degrees North, default 90)
           slat      = southern latitude limit (degrees North, defalut -90)
           linc      = number of latitude tick incriments (default 6)
           earth     = include Earth continent outlines (default False)
    '''

    # Initialize the z variable limits
    zmin, zmax = find_zlimits([gData], zkey, -2, 6)

    # Initialize the new figure

    pf = False
    sn = len(aindex)

    if(sn < 1):
        print "plot_mult_3D_slices ERROR: no altitude slices specified"

    if(string.lower(plot_type)=="polar" and not earth):
        pf = True

    f  = plt.figure()
    ax = list()
    tl = " "

    # Adjust the figure height to accomadate the number of subplots

    if(sn > 2):
        fheight = f.get_figheight()
        f.set_figheight(fheight * 0.5 * sn)

    if(sn > 1 and string.lower(plot_type) == "polar"):
        fwidth = f.get_figwidth()
        f.set_figwidth(fwidth * 0.5)

    for snum in reversed(range(0, sn)):
        cl   = False
        xl   = False
        yl   = False
        fnum = (sn * 100) + 11 + snum
        ax   = f.add_subplot(fnum, polar=pf)

        if(sn == snum + 1):
            xl = True

        if(math.floor(sn * 0.5) == snum):
            yl = True

        if(snum == 0):
            cl = True

        if(string.lower(plot_type)=="rectangular"):

            con = plot_rectangular_3D_global(ax, zkey, gData, zmin, zmax, 6,
                                             aindex[snum], nlat, slat, linc,
                                             False, "t", tl, "r", xl, yl, earth)
        elif(string.lower(plot_type)=="polar"):
            # Because the polar plots take up half the width as the rectangular
            # plots, decrease the level of latitude incriments by half
            # and adjust the y label output pad

            
            con = plot_polar_3D_global(ax, 2, zkey, gData, zmin, zmax, 6,
                                       aindex[snum], nlat, slat, linc/2, False,
                                       "t", tl, "l", xl, yl, earth)
        else:
            print "ERROR: unknown input type [", plot_type, "]\n"
            sys.exit(0)

        if(cl == False):
            cpr = list(con.ax.get_position().bounds)

        else:
            # Add and adjust colorbar

            cbar = add_colorbar(con, zmin, zmax, 6, "horizontal",
                                gData[zkey].attrs['scale'],
                                gData[zkey].attrs['name'],
                                gData[zkey].attrs['units'])

            bp  = list(cbar.ax.get_position().bounds)
            cp  = list(con.ax.get_position().bounds)

            cp[1] = bp[1]
            cp[3] = cpr[3]
            bp[1] = cp[1] + cp[3] + 0.025

            cbar.ax.set_position(bp)
            con.ax.set_position(cp)

    # Draw to screen.
    if plt.isinteractive():
        plt.draw() #In interactive mode, you just "draw".
    else:
        # W/o interactive mode, "show" stops the user from typing more 
        # at the terminal until plots are drawn.
        plt.show()

    # Save output file

    if figname is not None:
        plt.savefig(figname)

# End plot_mult_3D_slices

def plot_rectangular_3D_global(ax, zkey, gData, zmin, zmax, zinc=6, aindex=-1,
                               nlat=90, slat=-90, linc=6, cb=True, cloc="r",
                               title=None, tloc="t", xl=True, yl=True,
                               earth=False, *args, **kwargs):
    '''
    Creates a rectangular map projection plot for a specified latitude range.
    Input: ax        = axis handle
           zkey      = key for z variable (ie 'Vertical TEC')
           gData     = gitm bin structure
           zmin       = minimum value for z variable
           zmax       = maximum value for z variable
           zinc       = number of tick incriments for z variable (default 6)
           aindex    = altitude index (default -1, for 2D parameter)
           nlat      = northern latitude limit (degrees North, default 90)
           slat      = southern latitude limit (degrees North, defalut 90)
           linc      = number of tick incriments in latitude (default 6)
           cb        = Add a colorbar (default is True)
           cloc      = Specify the colorbar location (t=top, r=right, l=left,
                       b=bottom, default is right)
           title     = plot title (default is None)
           tloc      = Specify the title location (t=top, r=right, l=left,
                       b=bottom, default is top)
           xl        = Include x (longitude) label (default is True)
           yl        = Include y (latitude) label (default is True)
           earth     = Include Earth continent outlines (default is False)
    '''

    if(nlat == slat):
        print "plot_rectangular_3D_global ERROR: no latitude range"
        sys.exit(0)

    df = False

    if(aindex is -1):
        df     = True
        aindex = 0

    # Set latitude range values
    yrange = nlat - slat
    ywidth = yrange / linc

    # If desired, map the Earth using the Equidistant Cylindrical Projection
    if earth:
        m = Basemap(projection='cyl',llcrnrlat=slat,urcrnrlat=nlat,
                    llcrnrlon=0,urcrnrlon=360,resolution='i')
        m.drawcoastlines(linewidth=0.5)
        m.drawparallels(np.arange(slat, nlat+1, ywidth))
        m.drawmeridians(np.arange(0,361,60))

    # Set the contour
    v    = np.linspace(zmin, zmax, 70, endpoint=True)
    con  = ax.contourf(gData['dLon'][:,:,aindex], gData['dLat'][:,:,aindex],
                       gData[zkey][:,:,aindex], v, cmap=get_cmap('Spectral_r'),
                       vmin=zmin, vmax=zmax)

    # Configure axis
    ytics  = MultipleLocator(ywidth)
    ax.yaxis.set_major_locator(ytics)
    if yl:
        ax.set_ylabel('Latitude ($degrees$)')
    plt.ylim(slat, nlat)

    xtics = MultipleLocator(60)
    ax.xaxis.set_major_locator(xtics)
    if xl:
        ax.set_xlabel('Longitude ($degrees$)')
    plt.xlim(0., 360.)

    # Set the title
    if title:
        rot  = 'horizontal'
        yloc = 1.1
        xloc = .5

        if(tloc == "l" or tloc == "r"):
            xloc = -.1
            yloc = .5
            rot  = 'vertical'

            if(tloc == "r"):
                xloc = 1.1

        if(tloc == "b"):
            yloc = -.1
            
        if df:
            ax.set_title(r'%s' % (title), size='medium', rotation=rot, y=yloc,
                         x=xloc)
        elif title == " ":
            ax.set_title(r'%5.2f $km$' % (gData['Altitude'][0,0,aindex]/1000.0),
                         size='medium', rotation=rot, y=yloc, x=xloc)
        else:
            ax.set_title(r'%s slice at %5.2f $km$' % (title,
                                          gData['Altitude'][0,0,aindex]/1000.0),
                         size='medium', rotation=rot, y=yloc, x=xloc)

    # Add a colorbar
    if cb:
        orient = 'vertical'

        if(cloc == 't' or cloc == 'b'):
            orient = 'horizontal'

        cbar = add_colorbar(con, zmin, zmax, zinc, orient,
                            gData[zkey].attrs['scale'],
                            gData[zkey].attrs['name'],
                            gData[zkey].attrs['units'])

        if(cloc == 'l' or cloc == 't'):
            bp = list(cbar.ax.get_position().bounds)
            cp = list(con.ax.get_position().bounds)

            if(cloc == 't'):
                cp[1] = bp[1]
                bp[1] = cp[1] + cp[3] + 0.025
            else:
                bp[0] = 0.125
                cp[0] = bp[0] + 0.1 + bp[2]

            con.ax.set_position(cp)
            cbar.ax.set_position(bp)

    return con

#End plot_rectangular_3D_global


def plot_polar_3D_global(ax, nsub, zkey, gData, zmin, zmax, zinc=6, aindex=-1,
                         center_lat=90, edge_lat=0, linc=6, cb=True, cloc="r",
                         title = None, tloc="t", tl = True, rl = True,
                         earth = False, *args, **kwargs):
    '''
    Creates a single polar projection, with the latitude center and range
    determined by the input.
    Input: ax         = axis handle
           nsub       = number of subplots in this row
           zkey       = key for z variable (ie 'Vertical TEC')
           gData      = gitm bin structure
           zmin       = minimum value for z variable
           zmax       = maximum value for z variable
           zinc       = number of tick incriments for z variable (default 6)
           aindex     = altitude index (default -1, for 2D parameter)
           center_lat = upper (center) latitude limit (degrees N, default 90)
           edge_lat   = lower (edge) latitude limit (degrees N, default 0)
           linc       = number of tick incriments in latitude (default 6)
           cb         = Add a colorbar (default is True)
           cloc       = Colorbar location (t=top, r=right, l=left, b=bottom, 
                        default is right)
           title      = plot title (default is none)
           tloc       = title location (t=top, r=right, l=left, b=bottom,
                        default is top)
           tl         = include theta (longitude) label (default is True)
           rl         = include r (latitude) label (default is True)
           earth      = include Earth continent outlines (default is False)
    '''

    df = False

    if(aindex == -1):
        df     = True
        aindex = 0

    # Assign the Longitude, Latitude, and Z data structures
    r     = gData['dLat'][:,:,aindex]
    z     = gData[zkey][:,:,aindex]

    # Set range values
    rrange  = center_lat - edge_lat
    rwidth  = rrange / linc
    v       = np.linspace(zmin, zmax, 70, endpoint=True)

    if earth:
        theta = gData['dLon'][:,:,aindex]
        # If desired, map the Earth using a Polar Azimuthal Equidistant
        # Projection

        if center_lat < 0.0:
            m    = Basemap(projection='spaeqd', lon_0=270, lat_0=center_lat,
                           boundinglat=edge_lat, round=True, resolution='i')
            lats = list(np.arange(edge_lat, center_lat-1, rwidth))
            lon  = 157.5
        else:
            m    = Basemap(projection='npaeqd', lon_0=270, lat_0=center_lat,
                           boundinglat=edge_lat, round=True, resolution='i')
            lats = list(np.arange(edge_lat, center_lat+1, rwidth))
            lats.reverse()
            lon  = 22.5

        m.drawcoastlines(linewidth=0.5)
        m.drawmapboundary(linewidth=0.5)
        llab = m.drawparallels(lats)

        m.drawmeridians(np.arange(0,360,45), labels=[1,1,0,0], labelstyle="+/-")
        con  = m.contourf(theta, r, z, v, cmap=get_cmap('Spectral_r'),
                          vmin=zmin, vmax=zmax, latlon=True) 
        lpad = 20
        x, y = m(list(lon for i in lats), lats)

        for i, label in enumerate(lats):
            ax.text(x[i], y[i], "%.0f$^\circ$" % (label))

    else:
        # Set the contour
        theta = gData['Longitude'][:,:,aindex]
        con   = ax.contourf(theta, center_polar_cap(center_lat, edge_lat, r),
                            z, v, cmap=get_cmap('Spectral_r'), vmin=zmin,
                            vmax=zmax)
        lpad    = 0
        rtics   = [edge_lat + rwidth * (x) for x in range(linc+1)]
        rlabels = map(str, (map(int, rtics)))

        if(center_lat > edge_lat):
            rlabels.reverse()

        ax.set_rmax(edge_lat)
        ax.set_rmin(center_lat)
        ax.set_rticks(rtics)
    
        if(min(rtics) >= 0.0):
            ax.set_rgrids(rtics[1:], labels=rlabels[1:])

             
    # Configure axis.
    if tl:
        ax.set_xlabel('Longitude', labelpad=lpad)

    if rl:
        lpad         = 200 + 50 * (nsub-1)
        label_string = ""

        if earth:
            if cb or nsub > 1:
                lpad = 1.5 * lpad + 50 * (nsub - 1)
            else:
                lpad *= 2
        else:
            lpad = lpad - 10 * (nsub - 1)
            label_string = " ($degrees$)"

        ax.set_ylabel('Latitude%s' % label_string, labelpad=-lpad/nsub)
           
    # Set the title
    if title:
        rot  = 'horizontal'
        yloc = 1.1
        xloc = 0.5

        if tloc == "b":
            yloc = -.1
        elif tloc != "t":
            rot  = 'vertical'
            yloc = 0.5
            xloc = -.2

            if earth:
                xloc = -.3

            if tloc == "r":
                xloc = 1.1
                
                if earth:
                    xloc = 1.2

        if df:
            ax.set_title(r'%s' % (title), size='medium', y=yloc, x=xloc,
                         rotation = rot)
        elif title == " ":
            ax.set_title(r'%5.2f $km$'%(gData['Altitude'][0,0,aindex]/1000.0),
                         size='medium', y=yloc, x=xloc, rotation=rot)
        else:
            ax.set_title(r'%s slice at %5.2f $km$' % (title,
                                          gData['Altitude'][0,0,aindex]/1000.0),
                         size='medium', y=yloc, x=xloc, rotation=rot)

    # Add a colorbar
    if cb:
        orient = 'vertical'

        if(cloc == 't' or cloc == 'b'):
            orient = 'horizontal'

        cbar = add_colorbar(con, zmin, zmax, zinc, orient,
                            gData[zkey].attrs['scale'],
                            gData[zkey].attrs['name'],
                            gData[zkey].attrs['units'])

        if(cloc == 'l' or cloc == 't'):
            bp = list(cbar.ax.get_position().bounds)
            cp = list(con.ax.get_position().bounds)

            if(cloc == 't'):
                cp[1] = bp[1]
                bp[1] = cp[1] + cp[3] + 0.025
            else:
                bp[0] = 0.125
                cp[0] = bp[0] + 0.1 + bp[2]

            con.ax.set_position(cp)
            cbar.ax.set_position(bp)

    return con

def add_colorbar(contour_handle, zmin, zmax, zinc, orient, scale, name, units):
    '''
    Add a colorbar

    Input: contour_handle = handle to contour plot
           zmin           = minimum z value
           zmax           = maximum z value
           zinc           = z tick incriment (recommend 6)
           orient         = orientation of the colorbar (horizontal or vertical)
           scale          = linear or exponential?
           name           = z variable name
           units          = z variable units
    '''

    w  = np.linspace(zmin, zmax, zinc, endpoint=True)
    cb = plt.colorbar(contour_handle, ticks=w, pad=.15, orientation=orient,
                      fraction=.06)
    if(scale is "exponetial"):
        cb.formatter=FormatStrFormatter('%7.2E')
    cb.set_label(r'%s ($%s$)' % (name, units))
    cb.update_ticks()
    return cb

def center_polar_cap(rcenter, redge, r):
    '''
    Adjust the radial axis in a polar plot so that it is centered about
    the northern or southern pole
    '''

    if(rcenter > redge):
        return rcenter - r
    else:
        return r

def find_zlimits(gDataList, zkey, aindex=-1, zinc=6, *args, **kwargs):
    '''
    Establish the appropriate z-axis limits for a list of GitmBin files

    Input: gDataList = A list of GitmBin data structures
           zkey      = key for the desired z value
           aindex    = altitude index (default -1 for 2D measurement,
                       use -2 for no index)
           zinc      = number of tick incriments (default is 6)
    '''

    hold_min = []
    hold_max = []

    for gData in gDataList:
        if(aindex > -2):
            flat = gData[zkey][:,:,aindex].reshape(-1)
        else:
            flat = gData[zkey][:,:,:].reshape(-1)

        hold_min.append(min(flat))
        hold_max.append(max(flat))

    zmin = min(hold_min)
    zmax = max(hold_max)
    zran = round((zmax-zmin)/zinc)

    if(zran != 0.0):
        zmin = math.floor(zmin / zran) * zran
        zmax = math.ceil(zmax / zran) * zran

    return zmin, zmax


#End

