# pyrafspec

pyrafspec is a wrapper to use pyraf(iraf) to extract spectrum observed by E9G10 of BFOSC or YFOSC

install iraf could look at https://iraf-community.github.io/install.html

install pyraf (https://pypi.org/project/pyraf/)
- pip install pyraf

install laspec (https://github.com/hypergravity/laspec)
- pip install laspec

install lightkurve
- pip install lightkurve==1.11.3

install pyrafspec
- $git clone https://github.com/lidihei/pyrafspec.git
- $cd pyrafspec
- $python setup.py install

# pyfile/bfoscE9G10/E9G10.py is a example of extracting spectrum observed by Xinglong 216cm BFOSC E9+G10 (pyfile/bfoscE9G10)
# pyfile/bfoscE9G10/HRS240.py is a example of extracting spectrum observed by Lijiang 240cm HRS

# wavelength calibrate
- keywords of iraf
-- k --> down to next order 
-- j --> up to previous order or (check rms)
-- h --> leave rms
-- d --> delete points
-- m --> mark points
-- f --> fit curve
-- w&e -> e --> select zone (xaixs invert: right_bottom --> left_top; yaxis_invert: left_top-->right_bottom)
-- w&a --> restore image
