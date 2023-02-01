CREATE TABLE IF NOT EXISTS detfiles (
    detid SERIAL PRIMARY KEY,
    nightid INT,
    expid INT,
    fid INT,
    fieldid INT,
    filename VARCHAR(500),
    filesize REAL,
    checksum VARCHAR(50),
    CREATED REAL,
    CONSTRAINT fk_nights_det
            FOREIGN KEY(nightid)
                REFERENCES nights(nightid),

    CONSTRAINT fk_exps_det
            FOREIGN KEY(expid)
                REFERENCES exposures(expid),

    CONSTRAINT fk_filters_det
            FOREIGN KEY(fid)
                REFERENCES filters(fid),

    CONSTRAINT fk_fields_det
            FOREIGN KEY(fieldid)
                REFERENCES fields(fieldid)

);
