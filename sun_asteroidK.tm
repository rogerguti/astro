# The meta kernel file contains entries pointing to the following SPICE kernels, which the user needs to download.
#   https://naif.jpl.nasa.gov/pub/naif/CASSINI/kernels/spk/https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp
#   https://ssd.jpl.nasa.gov/x/spk.html
#   with data: 
#   Object: 2162173
#   SPK start date (TDB): 2010-Jan-01 00:01
#   SPK stop date (TDB): 2025-Jan-01 00:01
#   SPK file format: Binary

\begindata
KERNELS_TO_LOAD=(
'naif0012.tls.pc'
'de430.bsp',
'2162173.bsp')
\begintext
