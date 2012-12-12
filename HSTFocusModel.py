"""
Interfaces with the HST Focus Model hosted at the Space Telescope Science
Institute (see http://www.stsci.edu/hst/observatory/focus/FocusModel).

Because the model queries 5-minute interval temperature telemetry from 2003
onwards (much data), STScI provides it as a web service. The STScI website
provides a web form interface, but no programmatic interface, so it is
time-consuming to make many queries. This module builds HTTP POST requests
which are submitted to the STScI CGI script to generate the output data on
their server, and then retrives those output data.

Overview from the STScI website:

For a given time, this model estimates the amount of defocus at a particular
camera. The model is a function of telemetered temperatures and secular terms.
The temperature dependent part of the model is taken from Di Nino et al. 2008
while the long-term secular trend is a more recent determination given in
Niemi et al., 2010. In addition to the temperature terms and secular double
exponential, the model includes zero point offsets characterizing the focus
offsets between cameras and channels.
"""
__author__ = 'Matt Mechtley'
__copyright__  = '2012, Creative Commons Attribution-ShareAlike 3.0'
import httplib, urllib

_focus_tool_server = 'focustool.stsci.edu:80'
_focus_request_url = '/cgi-bin/control3.py'

_txt_table_file_fmt = '/images/focusdata{Year}.{Date}_{Start}-{Stop}.txt'
_png_plot_file_fmt = '/images/focusplot{Year}.{Date}_{Start}-{Stop}.png'

class HTTPResponseError(Exception):
	def __init__(self,response):
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
	form_controls = {'Output':'Model',
	                 'Year':year,
	                 'Camera':camera,
	                 'Date':date,
	                 'Start':start,
	                 'Stop':stop}
	## These are the form controls passed via POST on the website
	conn = httplib.HTTPConnection(_focus_tool_server)
	conn.request('POST', _focus_request_url,
	             urllib.urlencode(form_controls),
	             {'Content-type':'application/x-www-form-urlencoded'})
	response = conn.getresponse()
	if response.status != httplib.OK:
		raise HTTPResponseError(response)

	## Build URLs for the generated data files (table and/or image)
	filename_params = {'Year':year,
	                   'Date':date.replace('/', '.'),
	                   'Start':start.replace(':', ''),
	                   'Stop':stop.replace(':', '')}
	txt_table_url = _txt_table_file_fmt.format(**filename_params)
	png_plot_url = _png_plot_file_fmt.format(**filename_params)

	## Retrieve generated data files
	png_data, txt_data = None, None
	if format in ('TXT', 'BOTH'):
		conn.request('GET', txt_table_url, headers={'Accept': 'text/plain'})
		response = conn.getresponse()
		if response.status != httplib.OK:
			raise HTTPResponseError(response)
		txt_data = response.read()
	if format in ('PNG', 'BOTH'):
		conn.request('GET', png_plot_url, headers={'Accept': 'image/png'})
		response = conn.getresponse()
		if response.status != httplib.OK:
			raise HTTPResponseError(response)
		png_data = response.read()

	conn.close()

	if format == 'PNG':
		return png_data
	elif format == 'TXT':
		return txt_data
	else:
		return txt_data, png_data

# Can also run from the command line. First argument should be date in
# YYYY/MM/DD format, second argument is time range in HH:MM-HH:MM 24-hour
# format
if __name__ == '__main__':
	import sys
	year,date = sys.argv[1].split('/',1)
	start,stop = sys.argv[2].split('-')
	print get_model_data(year,date,start,stop)
