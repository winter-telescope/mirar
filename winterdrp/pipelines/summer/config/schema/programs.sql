CREATE TABLE IF NOT EXISTS programs (
    progid SERIAL PRIMARY KEY,
    progname VARCHAR(20),
    piname VARCHAR(20),
    startdate REAL,
    enddate REAL,
    time_hours REAL,
    basepriority REAL
);
