# winter_drp
Requirements and installations - 
Sextractor -
1. Download the source code from the repo https://github.com/astromatic/sextractor
2. Follow the instructions here https://sextractor.readthedocs.io/en/latest/Installing.html
3. It will ask you to install all sorts of weird libraries that may or may not be compatible with your device. The configure that worked best for me is without the atlas directories but by using the openblas files instead
./configure --enable-openblas --with-atlas-libdir=/usr/local/opt --with-atlas-incdir=/usr/local/opt/openblas/include

Instructions for running the focusloop:
1. Download all images to a data directory
2. Download the log file (should be in the format of log_20210625.csv)
3. Run python fix_headers.py --d <data_dir_path> --l <logfilename>
4. Run python focusLoop.py --d <data_dir_path> --plot
This will go through all files in the data directory to check for files that are "focus" files and do the analysis on them, and plot the results in a file named "focusloop.pdf'

Instructions for .env
1. Copy the `.env.example` file to the root of the project and update the environment variables.
2. Name this file `.env`
