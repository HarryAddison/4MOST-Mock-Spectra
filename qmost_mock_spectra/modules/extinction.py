from dustmaps.sfd import SFDQuery
from astropy.coordinates import SkyCoord


def get_mwebv(ra, dec):
    '''
    Obtain the Milky way line of sight extinction at the
    given celestial frame coordinates.

    Inputs:
        > "ra" = Right ascension (celestial frame) in degrees.
        > "dec" = Declination (celestial frame) in degrees.

    Output:
        > "mwebv" = Milky Way line of sight extinction (B-V).
    '''

    # Store coordinates in astropy SkyCoord object.
    coord = SkyCoord(ra, dec, frame='icrs')

    # Query the dust map at the given coords
    mwebv = float(SFDQuery()(coord))

    return mwebv
