CREATE TABLE IF NOT EXISTS exposures ( expid SERIAL PRIMARY KEY, obsdate INT, obsHistID INT, itid INT, night INT,
        fieldID INT, filter VARCHAR(5), progID INT, visitTime FLOAT, visitExpTime FLOAT, expMJD FLOAT,
        subprogram VARCHAR(20), airmass FLOAT, shutopen REAL, shutclsd REAL, tempture REAL, windspd REAL, Dewpoint REAL,
        Humidity REAL, Pressure REAL, Moonra REAL, Moondec REAL,
        Moonillf REAL, Moonphas REAL, Moonalt REAL, Sunaz REAL, Sunalt REAL, Detsoft VARCHAR(50), Detfirm
        VARCHAR(50), ra REAL, dec REAL, altitude REAL, azimuth REAL, process_flag INT
);