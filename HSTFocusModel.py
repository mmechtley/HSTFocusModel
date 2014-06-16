"""
Interfaces with the HST Focus Model hosted at the Space Telescope Science
Institute (see http://www.stsci.edu/hst/observatory/focus/FocusModel).

STScI provides it as a web service with a form interface, but no programmatic
interface, so it is time-consuming to make many queries. This module builds
HTTP POST requests which are submitted to the STScI CGI script to generate the
output data on their server, and then retrives those output data.

Overview from the STScI website:

For a given time, this model estimates the amount of defocus at a particular
camera. The model is a function of telemetered temperatures and secular terms.
The temperature dependent part of the model is taken from Di Nino et al. 2008
while the long-term secular trend is a more recent determination given in
Niemi et al., 2010. In addition to the temperature terms and secular double
exponential, the model includes zero point offsets characterizing the focus
offsets between cameras and channels.
"""
from __future__ import division
import httplib
import urllib
import re
from datetime import date
from StringIO import StringIO
from numpy import floor, genfromtxt, linspace, trapz
from scipy.interpolate import UnivariateSpline

__author__ = 'Matt Mechtley'
__copyright__ = '2012, Creative Commons Attribution-ShareAlike 3.0'

_server = 'focustool.stsci.edu:80'
_request_url = '/cgi-bin/control3.py'

_txt_table_file_fmt = '/images/focusdata{Year}.{Date}_{Start}-{Stop}.txt'
_png_plot_file_fmt = '/images/focusplot{Year}.{Date}_{Start}-{Stop}.png'

# Formatted in the way numpy.genfromtxt expects
_model_output_columns = 'JulianDate, Month, Day, Year, Time, Model'

# Field widths in the text download format. Sadly the output is not delimiter-
# separated
_output_field_widths = (12, 4, 3, 5, 9, 8)


class HTTPResponseError(Exception):
    def __init__(self, url, response):
        self.url = url
        self.response = response
        self.message = 'Bad response from server\n{} {}'.format(
            response.status, response.reason)


def get_model_data(year, date, start, stop, camera='UVIS1', format='TXT'):
    """
    Gets plaintext-formatted table or png-formatted image of model focus data
    for the specified time range.

    Uses a two-step HTTP request process. The first request causes the server
    to generate the output files, and the second request retrieves the table or
    image.

    Note: Make sure you request sensible times for a given camera, or you may
    get nonsense results. e.g. WFPC2/PC wasn't on the telescope in 2010, but
    the focus model will happily generate model values for that period. Valid
    dates for each camera are given at http://focustool.stsci.edu

    :param year: Year of observation (string or integer)
    :param date: Date of observation (string in MM/DD format)
    :param start: Start of time range (string in 24-hr HH:MM format)
    :param stop: End of time range (string in 24-hr HH:MM format)
    :param camera: One of UVIS1, UVIS2, WFC1, WFC2, HRC, PC. Default is UVIS1.
    :param format: Output format, one of TXT, PNG, or BOTH. Default is TXT.
    :return: model focus data, as plaintext file or png image as requested, or
             2-tuple of text, image if format is BOTH
    """
    form_controls = {'Output': 'Model',
                     'Year': year,
                     'Camera': camera,
                     'Date': date,
                     'Start': start,
                     'Stop': stop}
    ## These are the form controls passed via POST on the website
    conn = httplib.HTTPConnection(_server)
    conn.request('POST', _request_url,
                 urllib.urlencode(form_controls),
                 {'Content-type': 'application/x-www-form-urlencoded'})
    response = conn.getresponse()
    conn.close()
    if response.status != httplib.OK:
        raise HTTPResponseError(_request_url, response)

    ## Build URLs for the generated data files (table and/or image)
    filename_params = {'Year': year,
                       'Date': date.replace('/', '.'),
                       'Start': start.replace(':', ''),
                       'Stop': stop.replace(':', '')}
    txt_table_url = _txt_table_file_fmt.format(**filename_params)
    png_plot_url = _png_plot_file_fmt.format(**filename_params)

    ## Retrieve generated data files
    png_data, txt_data = None, None
    if format in ('TXT', 'BOTH'):
        conn = httplib.HTTPConnection(_server)
        conn.request('GET', txt_table_url, headers={'Accept': 'text/plain'})
        response = conn.getresponse()
        if response.status != httplib.OK:
            raise HTTPResponseError(txt_table_url, response)
        txt_data = response.read()
        conn.close()
    if format in ('PNG', 'BOTH'):
        conn = httplib.HTTPConnection(_server)
        conn.request('GET', png_plot_url, headers={'Accept': 'image/png'})
        response = conn.getresponse()
        if response.status != httplib.OK:
            raise HTTPResponseError(png_plot_url, response)
        png_data = response.read()
        conn.close()

    if format == 'PNG':
        return png_data
    elif format == 'TXT':
        return txt_data
    else:
        return txt_data, png_data


def mean_focus(expstart, expend, camera='UVIS1', spline_order=3,
               not_found_value=None, with_var=False):
    """
    Gets the mean focus over a given observation period. Exposure start and end
    times can be specified as Modified Julian Date float (like the FITS header
    EXPSTART and EXPEND keywords) or a UTC time string in YYYY-MM-DD HH:MM:SS
    format.
    :param expstart: Start time of exposure.
    :param expend: End time of exposure.
    :param camera: One of UVIS1, UVIS2, WFC1, WFC2, HRC, PC. Default is UVIS1.
    :param spline_order: Degree of the spline used to interpolate the model
     data points (passed as k= to scipy.interpolate.UnivariateSpline). Use 1 for
     linear interpolation. Default is 3.
    :param not_found_value: Value to return if the Focus Model does not have
     data for the given time interval. Default value (None) means raise
     HTTPResponseError
     :param with_var: Also include variance in a returned 2-tuple
    :return: Continuous (integral) mean focus between expstart and expend
    """
    # Convert date/time strings to MJD
    try:
        startnums = [int(num) for num in re.split(':|-|/| ', expstart)]
        endnums = [int(num) for num in re.split(':|-|/| ', expend)]
        expstart = _date_time_to_mjd(*startnums)
        expend = _date_time_to_mjd(*endnums)
    except TypeError:
        pass
    # Pad input exposure start and end time, to make sure we get at least one
    # data point before and after. Then split up into year, date, times
    ten_mins = 10 / (24 * 60)
    expstart_pad = float(expstart) - ten_mins
    expend_pad = float(expend) + ten_mins
    start_yr, start_date, start_time = _mjd_to_year_date_time(expstart_pad)
    stop_yr, stop_date, stop_time = _mjd_to_year_date_time(expend_pad)
    # Chop off seconds
    start_time = start_time.rsplit(':', 1)[0]
    stop_time = stop_time.rsplit(':', 1)[0]

    if start_date != stop_date:
        intervals = [(start_yr, start_date, start_time, '23:59'),
                     (stop_yr, stop_date, '00:00', stop_time)]
    else:
        intervals = [(start_yr, start_date, start_time, stop_time)]

    try:
        txt_focus = ''
        # Get text table of focus data for each interval
        for year, date, start, stop in intervals:
            txt_interval = get_model_data(year, date, start, stop, camera,
                                          format='TXT')
            col_names, txt_interval = txt_interval.split('\n', 1)
            txt_focus += txt_interval
        # convert to numpy array
        focus_data = genfromtxt(StringIO(txt_focus), skiprows=0, dtype=None,
                                names=_model_output_columns,
                                delimiter=_output_field_widths)
        # Create interpolating spline
        spline = UnivariateSpline(focus_data['JulianDate'], focus_data['Model'],
                                  k=spline_order)
        # Return the continuous (integral) mean
        mean_foc = spline.integral(expstart, expend) / (expend - expstart)
        # Calculate signal variance (see e.g. Wikipedia article for RMS)
        if with_var:
            xvals = linspace(expstart, expend, focus_data.size*2)
            var_foc = trapz(spline(xvals)**2, xvals) / (expend - expstart)
            var_foc -= mean_foc**2

    except HTTPResponseError, err:
        if err.response.status == httplib.NOT_FOUND \
                or not_found_value is not None:
            mean_foc = not_found_value
            var_foc = not_found_value
        else:
            raise err

    if with_var:
        return mean_foc, var_foc
    else:
        return mean_foc


def add_mean_focus_to_header(filename, ext=0, with_var=False, focus_key='FOCUS',
                             var_key='FOCUSVAR', **kwargs):
    """
    Calculates the mean focus for the given fits file (based on EXPSTART and
    EXPEND header keywords), and saves the calculated focus into the header.
    :param filename: Filename to calculate focus for
    :param ext: FITS extension that contains EXPSTART and EXPEND keywords
    :param with_var: Should variance be calculated and saved as well?
    :param focus_key: Keyword in which to save focus value
    :param var_key: Keyword in which to save focus variance value
    :param kwargs: Additional keyword arguments, passed to mean_focus()
    """
    try:
        import pyfits
    except ImportError, err:
        print 'pyfits module is required for reading and writing fits headers'
        pyfits = None
        raise err
    expstart = pyfits.getval(filename, 'EXPSTART', ext=ext)
    expend = pyfits.getval(filename, 'EXPEND', ext=ext)
    focus = mean_focus(expstart, expend, with_var=with_var, **kwargs)
    if not hasattr(focus, 'count'):
        focus = (focus, )
    comments = ('Estimated mean focus (HST Focus Model)',
                'Estimated focus variance (HST Focus Model)')
    for key, val, comment in zip((focus_key, var_key), focus, comments):
        pyfits.setval(filename, key, value=val, comment=comment)


def _mjd_to_year_date_time(mjd):
    """
    Convert Modified Julian Date number to (year, date, time) tuple
    :param mjd: Modified Julian Date (e.g. from fits header)
    :return: 3-tuple of (YYYY, MM/DD, HH:MM:SS) (int, string, string)
    """
    # Adapted from libastro/XEphem included in pyephem, which erroneously
    # calls its dates Modified Julian Date, when in fact they are Dublin
    # Julian Date (Dec 31 1899 12:00 zero point).
    # https://github.com/brandon-rhodes/pyephem/blob/master/libastro-3.7.5/mjd.c

    _mjd_to_dublin = -15019.5  # Conversion from MJD to Dublin JD (Wikipedia)
    _days_per_year = 365.25  # In JD convention, anyway

    # For Gregorian calendar days.
    _gregorian_change_date = -115860.0  # Date of Julian-Gregorian change in DJD
    # 14.99835726 leap days skipped up to Dec 31 1899 12:00
    _leapskips_epoch = 14.99835726
    _days_per_century = 36524.25

    days1900 = mjd + _mjd_to_dublin + 0.5  # +12h moves ref pt to 1/1/1900 00:00
    day_int = floor(days1900)
    day_frac = days1900 - day_int

    if day_frac == 1:
        day_frac = 0
        day_int += 1

    # Is the date a Gregorian Calendar date?
    if day_int > _gregorian_change_date:
        leapskips = floor((day_int / _days_per_century) + _leapskips_epoch)
        # Add skipped leap days, but subtract off for 400-divisible years
        day_int += 1 + leapskips - floor(leapskips / 4.0)

    ## TODO: Give variables semantic names
    b = floor((day_int / _days_per_year) + .802601)
    ce = day_int - floor((_days_per_year * b) + .750001) + 416
    g = floor(ce / 30.6001)
    month = int(g - 1)
    day = int(ce - floor(30.6001 * g) + day_frac)
    year = int(b + 1899)

    if g > 13.5:
        month = int(g - 13)
    if month < 2.5:
        year = int(b + 1900)
    if year < 1:
        year -= 1

    hour = int(24 * day_frac)
    minute = int((day_frac - hour / 24) * 24 * 60)
    second = int((day_frac - hour / 24 - minute / (24 * 60)) * 24 * 3600)

    return (year, '{0:02d}/{1:02d}'.format(month, day),
            '{0:02d}:{1:02d}:{2:02d}'.format(hour, minute, second))


def _date_time_to_mjd(year, month, day, hours, minutes, seconds):
    """
    Convert a calendar date to a Modified Julian Day Number
    """
    # Adapted from pyephem/libastro, which uses Dublin Julian Date (different
    # zero point from MJD). So at the end, we add a conversion factor
    _last_gregorian_date = date(1582, 10, 14)
    _dublin_to_modified = 15019.5
    m = month
    y = year + 1 if year < 0 else year
    if month < 3:
        m += 12
        y -= 1

    # Handle differently before and after Julian-Gregorian changeover
    if date(year, month, day) < _last_gregorian_date:
        b = 0
    else:
        n_centuries = y // 100
        b = 2 - n_centuries + n_centuries // 4

    # Handle differently before and after epoch zero point
    if y < 0:
        c = ((365.25 * y) - 0.75) - 694025
    else:
        c = (365.25 * y) - 694025

    d = 30.6001 * (m + 1)

    dublin_jd = b + c + d + day - 0.5

    dublin_jd += hours / 24 + minutes / 1440 + seconds / 86400

    return dublin_jd + _dublin_to_modified
