CREATE TABLE IF NOT EXISTS detfiles ( detid SERIAL PRIMARY KEY,
nightid INT FOREIGN KEY REFERENCES nights(nightid),
    expid INT FOREIGN KEY REFERENCES exposures(expid), fid INT FOREIGN KEY REFERENCES filters(fid), fieldid FOREIGN KEY
    REFERENCES fields(fieldid), filename VARCHAR(500), filesize REAL,  checksum VARCHAR(50), CREATED REAL);