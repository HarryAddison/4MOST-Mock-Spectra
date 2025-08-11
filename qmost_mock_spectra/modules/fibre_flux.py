import numpy as np
import math
from scipy import integrate
import scipy.ndimage


def effective_fibre_mag(sep, gal_ddlr, sersic_index, galmag, pix_size_arcsec, seeing):

    '''This code calculates the effective galaxy magnitude in a 4MOST fibre. I'd touch this
    as little as possible.
    
    sep : host-transient separation in arcseconds.
    gal_ddlr : normalised directional light radius of the galaxy (included in SNANA pop sims)
    sersic index : float, set to 0.5 for my purposes
    galmag : input galaxy magnitude, LSST r-band
    pix_size_arcsec : pixel scale, I find snsep/100 works well without being too slow
    seeing : seeing in arcsecond, I use 0.8 for everything'''

    dlr = gal_ddlr * sep
    gmag = galmag
    seperation = sep
    
    #first get total galaxy flux in erg/s/cm2 (constant here is the zero-point flux in AB system in the LSST r-band filter bandpass)
    flux_gal_total = 4.69542e-6 * (np.e ** (-gmag / 2.5))

    #then calculate effective intensity at centre of fibre (seperation)
    bn = (2 * sersic_index) - (1/3) + (0.009876 * sersic_index)  #from Prugneil and Simien

    const = 2 * sersic_index * (np.e ** bn) * (bn ** -(2 * sersic_index)) * math.factorial(int(((2 * sersic_index) - 1)))
    eff_intensity = flux_gal_total / (const * np.pi * (dlr ** 2))


    #now create the pixel array from integer number of pixels in display range
    pixel_length = pix_size_arcsec
    pixel_no = int((6 * seperation) / pixel_length) + 1 #plus one because always rounds down

    int_array = []
    distance_array = []
    for i in range(pixel_no):
        sub_array = []
        sub_dist = []
        for j in range(pixel_no):
            sub_array.append(0)
            sub_dist.append(np.sqrt((i - (pixel_no/2))**2 + (j - (pixel_no/2))**2) * pixel_length)
        
        int_array.append(sub_array)
        distance_array.append(sub_dist)

    #define the intensities at each pixel
    for h in range(pixel_no):
        for g in range(pixel_no):
            try:
                int_array[h][g] = eff_intensity * np.e ** (-bn * (((distance_array[h][g] / dlr) ** (1/sersic_index)) - 1))
            except ZeroDivisionError:
                int_array[h][g] = 0.0

    #define min and max intensities used in iamge (not max and min overall) and other useful numbers, objects
    max_int = max(max(int_array))
    min_int = min(min(int_array))

    #define conv radius using seeing (FWHM -> sigma)
    conv_radius = (seeing / np.sqrt(8 * np.log(2))) / pixel_length
    new_int_array = scipy.ndimage.gaussian_filter(int_array, sigma=conv_radius, mode='constant', cval=min_int)

    pixel_sep = seperation / pixel_length
    pixel_fibre_size = 0.725 / pixel_length
    true_array_int = np.asarray(new_int_array)

    #now try and get the pixels that are in the fibre. pixel coords are to the bottom left corner, but centres are at integer coords
    #remember circle of fibre is N/2 pixels up and to right 
    #also remember that the pixel intenisties are logged for display purposes and must have np.e raised to their power to get actual intensity values!!!!!
    good_coords = []
    fibre_int_pix = 0
    for h in range(pixel_no):
        for g in range(pixel_no):
            pixel_centre_x = h
            pixel_centre_y = g
            pixel_distance = (pixel_centre_x - ((pixel_no/2) + pixel_sep)) ** 2 + (pixel_centre_y - (pixel_no/2)) ** 2
            pixel_int = true_array_int[h][g]
            if pixel_distance <= pixel_fibre_size ** 2:
                # print((np.e ** pixel_int) * pixel_length * pixel_length, 'this is the pixel flux for pixel', h, g)
                fibre_int_pix = fibre_int_pix + (pixel_int * pixel_length * pixel_length)
                good_coords.append([h,g])
            else:
                continue

    #change the color off the intensities inside the fibre for gra hics purposes (note that int array seems to be in configuration y,x so that's odd)
    #should give coords as [0],[1] so don't know why that isn't working
    for i in range(len(good_coords)):
        coords = good_coords[i]
        int_array[coords[1]][coords[0]] = min_int


    ratio_pix = fibre_int_pix / flux_gal_total

    if ratio_pix > 1:
        ratio_pix = 1
    else:
        ratio_pix = ratio_pix

    #now calculate effective galaxy mag
    eff_mag = gmag - 2.5 * np.log10(ratio_pix)

    return eff_mag


def point_convolute(seeing, sne_mag):
    '''Essentially the same as the above function, but just for transient. Since
    they're point sources it's basically a constant adjustment based on seeing.'''
    print(sne_mag)
    #first turn the seeing value into a sigma for the gaussian
    FWHM = seeing
    sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))

    #generate a radial guassian (gaussian multiplied by 2 pi x)
    gaussian = lambda x:2 * np.pi * x * np.e ** (-(x**2)/(2 * (sigma ** 2)))

    #integrate to infinity to get the normalisation
    normalisation = integrate.quad(gaussian, 0, np.inf)

    #integrate again to fibre radius with normalisation constant dividing gaussian
    gaussian_norm = lambda x:2 * np.pi * x * np.e ** (-(x**2)/(2 * (sigma ** 2))) / normalisation[0]
    fraction = integrate.quad(gaussian_norm ,0 ,0.725)[0] #only want the first bit here

    print(sne_mag)
    print(2.5*np.log(fraction))
    new_mag = sne_mag - 2.5*np.log(fraction)

    return new_mag
