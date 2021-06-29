# winter_drp
Requirements and installations - 
Sextractor -
1. Download the source code from the repo https://github.com/astromatic/sextractor
2. Follow the instructions here https://sextractor.readthedocs.io/en/latest/Installing.html
3. It will ask you to install all sorts of weird libraries that may or may not be compatible with your device. The configure that worked best for me is without the atlas directories but by using the openblas files instead
./configure --enable-openblas --with-atlas-libdir=/usr/local/opt --with-atlas-incdir=/usr/local/opt/openblas/include
